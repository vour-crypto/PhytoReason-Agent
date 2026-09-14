"""MassBank, MoNA and PubChem adapters for Phase B2.

Network responses are kept as evidence.  A failed remote call never becomes a
successful empty result; callers receive ``partial`` with an explicit warning.
"""
from __future__ import annotations

import hashlib
import json
import threading
import time
from typing import Any, Iterable

import requests

from phyto_reason.models.evidence import Evidence, EvidenceType
from phyto_reason.models.tool_result import ToolResult
from phyto_reason.tools.base_tool import BaseTool, ToolParameter
from phyto_reason.tools.tool_registry import register_tool
from phyto_reason.utils.cache import knowledge_cache

MASSBANK_URL = "https://massbank.eu/MassBank-api/records"
MONA_URL = "https://mona.fiehnlab.ucdavis.edu/rest/spectra/search"
PUBCHEM_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound"


def _error_message(exc: Exception) -> str:
    text = str(exc)
    if "10013" in text or "PermissionError" in text:
        return f"疑似本机安全软件拦截，建议检查白名单: {text}"
    if isinstance(exc, (requests.Timeout, requests.ReadTimeout)) or "timed out" in text.lower():
        return f"网络较慢或不可达: {text}"
    return text


def parse_massbank_records(payload: Any) -> list[dict[str, Any]]:
    """Parse the observed MassBank shape, including a truncated JSON prefix."""
    if isinstance(payload, dict):
        if isinstance(payload.get("records"), list):
            return [x for x in payload["records"] if isinstance(x, dict)]
        return [payload]
    if isinstance(payload, list):
        return [x for x in payload if isinstance(x, dict)]
    if isinstance(payload, str):
        try:
            value = json.loads(payload)
            return parse_massbank_records(value)
        except json.JSONDecodeError:
            # The supplied sample is a stream cut mid-record. Recover complete
            # leading objects without pretending the complete response exists.
            decoder = json.JSONDecoder()
            pos = 0
            while pos < len(payload) and payload[pos].isspace():
                pos += 1
            if pos < len(payload) and payload[pos] == "[":
                pos += 1
            records: list[dict[str, Any]] = []
            while pos < len(payload):
                while pos < len(payload) and payload[pos] in " \r\n,":
                    pos += 1
                try:
                    value, end = decoder.raw_decode(payload, pos)
                except json.JSONDecodeError:
                    break
                if isinstance(value, dict):
                    records.append(value)
                pos = end
            return records
    return []


def _massbank_peak_values(record: dict[str, Any]) -> list[dict[str, float]]:
    values = (((record.get("peak") or {}).get("peak") or {}).get("values") or [])
    return [
        {"mz": float(v["mz"]), "intensity": float(v.get("intensity", 0)), "rel": float(v.get("rel", 0))}
        for v in values if isinstance(v, dict) and "mz" in v
    ]


def _massbank_links(record: dict[str, Any]) -> dict[str, list[str]]:
    links: dict[str, list[str]] = {}
    for item in ((record.get("compound") or {}).get("link") or []):
        if not isinstance(item, dict):
            continue
        db = str(item.get("database", "")).upper()
        identifier = str(item.get("identifier", ""))
        if db and identifier:
            links.setdefault(db, []).append(identifier)
    return links


def _massbank_record_summary(record: dict[str, Any]) -> dict[str, Any]:
    compound = record.get("compound") or {}
    acquisition = record.get("acquisition") or {}
    ms = record.get("mass_spectrometry") or {}
    focused = {str(x.get("subtag", "")): x.get("value") for x in (ms.get("focused_ion") or []) if isinstance(x, dict)}
    return {
        "accession": record.get("accession", ""),
        "title": record.get("title", ""),
        "compound": {
            "names": compound.get("names", []), "formula": compound.get("formula"),
            "mass": compound.get("mass"), "links": _massbank_links(record),
        },
        "ion_mode": ((acquisition.get("mass_spectrometry") or {}).get("ion_mode")),
        "ms_type": ((acquisition.get("mass_spectrometry") or {}).get("ms_type")),
        "precursor_mz": focused.get("PRECURSOR_M/Z"),
        "precursor_type": focused.get("PRECURSOR_TYPE"),
        "peaks": _massbank_peak_values(record),
    }


def _read_stream_json(response: requests.Response, max_bytes: int = 64 * 1024 * 1024) -> tuple[Any, int]:
    chunks: list[bytes] = []
    total = 0
    for chunk in response.iter_content(chunk_size=64 * 1024):
        if not chunk:
            continue
        total += len(chunk)
        if total > max_bytes:
            raise ValueError(f"response exceeds {max_bytes} bytes")
        chunks.append(chunk)
    raw = b"".join(chunks).decode(response.encoding or "utf-8", errors="replace")
    return raw, total


class _RemoteTool(BaseTool):
    supported_data_types = ["spectrum", "metabolite"]

    def _evidence(self, source: str, query: str, count: int, records: list[dict[str, Any]], warning: str = "") -> ToolResult:
        evidence = Evidence(
            evidence_type=EvidenceType.ANNOTATION, source=source,
            score=min(count / 10.0, 1.0), confidence=min(count / 10.0, 1.0),
            description=f"{source}: {count} records for {query}",
            metadata={"query": query, "records": records[:10], "actual_source": source},
        )
        return ToolResult(status="partial" if warning else "success", evidence_list=[evidence],
                          warnings=[warning] if warning else [],
                          metadata={"source": source, "query": query, "count": count, "records": records})


@register_tool
class MassBankSpectrumSearchTool(_RemoteTool):
    tool_name = "massbank_spectrum_search"
    version = "1.0.0"
    description = "MassBank3 化合物/质量/峰列表检索（流式读取）"
    parameters = [
        ToolParameter("compound_name", "string", "Compound name", required=False),
        ToolParameter("precursor_mz", "number", "Precursor m/z", required=False),
        ToolParameter("tolerance_ppm", "number", "Mass tolerance", default=10),
        ToolParameter("ion_mode", "string", "POSITIVE or NEGATIVE", required=False),
        ToolParameter("ms_type", "string", "MS1/MS2", required=False),
        ToolParameter("peaks", "array", "m/z,intensity pairs", required=False),
        ToolParameter("limit", "integer", "Requested record limit", default=10),
    ]

    def validate_input(self, **kwargs) -> list[str]:
        return [] if any(kwargs.get(k) not in (None, "", []) for k in ("compound_name", "precursor_mz", "peaks")) else ["需要 compound_name、precursor_mz 或 peaks"]

    def run(self, **kwargs) -> ToolResult:
        query = str(kwargs.get("compound_name") or kwargs.get("precursor_mz") or "peaks")
        limit = min(max(int(kwargs.get("limit", 10)), 1), 10)
        params: dict[str, Any] = {"limit": limit}
        if kwargs.get("compound_name"): params["compound_name"] = kwargs["compound_name"]
        if kwargs.get("precursor_mz") is not None:
            params["exact_mass"] = kwargs["precursor_mz"]
            params["mass_tolerance"] = float(kwargs.get("tolerance_ppm", 10))
        if kwargs.get("ion_mode"): params["ion_mode"] = str(kwargs["ion_mode"]).upper()
        if kwargs.get("ms_type"): params["ms_type"] = kwargs["ms_type"]
        if kwargs.get("peaks"):
            params["peak_list"] = ",".join(f"{float(p[0])};{float(p[1])}" if isinstance(p, (list, tuple)) else str(p) for p in kwargs["peaks"])
        last: Exception | None = None
        for attempt in range(3):
            try:
                response = requests.get(MASSBANK_URL, params=params, timeout=(30, 120), stream=True)
                response.raise_for_status()
                raw, nbytes = _read_stream_json(response)
                parsed = parse_massbank_records(raw)
                # 服务端不尊重 limit/page（2026-09-13 实测 limit=1/5、page=0/1
                # 返回字节级一致的全量体）；分页只能在客户端完成。
                seen: set[str] = set()
                records: list[dict[str, Any]] = []
                duplicates = 0
                for record in (_massbank_record_summary(x) for x in parsed):
                    accession = str(record.get("accession", ""))
                    if accession:
                        if accession in seen:
                            duplicates += 1
                            continue
                        seen.add(accession)
                    records.append(record)
                total_count = len(records)
                warning = ""
                if total_count > limit:
                    warning = (
                        f"MassBank 服务端不执行 limit（实测返回全量 {total_count} 条）；"
                        f"已在客户端截取前 {limit} 条"
                    )
                    records = records[:limit]
                result = self._evidence("MassBank", query, len(records), records, warning)
                result.metadata.update({
                    "requested_limit": limit,
                    "raw_bytes": nbytes,
                    "total_count": total_count,
                    "has_more": total_count > limit,
                    "returned": len(records),
                    "duplicates_dropped": duplicates,
                    "pagination": "client_side_only_server_ignores_limit_page",
                })
                return result
            except Exception as exc:
                last = exc
                if "10013" in str(exc):
                    break
                if attempt < 2:
                    time.sleep(2 ** (attempt + 1))
        return self._evidence("MassBank", query, 0, [], f"MassBank 请求失败: {_error_message(last or RuntimeError('unknown error'))}")


@register_tool
class MonaSpectrumSearchTool(_RemoteTool):
    tool_name = "mona_spectrum_search"
    version = "1.0.0"
    description = "MoNA Spring Filter 谱图检索"
    parameters = [
        ToolParameter("query", "string", "Spring Filter query", required=True),
        ToolParameter("size", "integer", "Page size", default=10),
        ToolParameter("page", "integer", "Zero-based page", default=0),
    ]

    def validate_input(self, **kwargs) -> list[str]:
        return ["query 不能为空（MoNA 空查询返回 400）"] if not str(kwargs.get("query", "")).strip() else []

    def run(self, **kwargs) -> ToolResult:
        query = str(kwargs.get("query", "")).strip()
        if not query:
            return self._evidence("MoNA", "", 0, [], "MoNA 查询必须至少包含一个 Spring Filter 条件")
        size = min(max(int(kwargs.get("size", 10)), 1), 100)
        page = max(int(kwargs.get("page", 0)), 0)
        try:
            response = requests.get(MONA_URL, params={"query": query, "size": size, "page": page}, headers={"Accept": "application/json"}, timeout=30)
            response.raise_for_status()
            data = response.json()
            rows = data if isinstance(data, list) else data.get("content", data.get("results", data.get("spectra", [])))
            records = [x for x in rows if isinstance(x, dict)] if isinstance(rows, list) else []
            result = self._evidence("MoNA", query, len(records), records)
            result.metadata.update({"page": page, "size": size, "raw_shape": type(data).__name__})
            return result
        except requests.HTTPError as exc:
            hint = ""
            response = getattr(exc, "response", None)
            if response is not None:
                try:
                    hint = str(response.json().get("hint", ""))
                except Exception:
                    hint = ""
            warning = f"MoNA 请求被拒绝 ({response.status_code if response is not None else 'HTTP'}): {_error_message(exc)}"
            if hint:
                warning += f"；服务端语法提示: {hint}"
            return self._evidence("MoNA", query, 0, [], warning)
        except Exception as exc:
            return self._evidence("MoNA", query, 0, [], f"MoNA 请求失败: {_error_message(exc)}")


_PUBCHEM_FAILURES = 0
_PUBCHEM_OPEN_UNTIL = 0.0
_PUBCHEM_LAST_PROBE = 0.0
_PUBCHEM_LOCK = threading.Lock()
def _local_properties(identifier: str) -> dict[str, Any]:
    """Return only facts present in the shared compound profile knowledge base.

    Profiles are MS/MS diagnostic rules, not a complete chemical registry. Do
    not manufacture CAS, SMILES, or other properties that the profile lacks.
    """
    from phyto_reason.knowledge.compound_profiles import find_diagnostic_rules

    rule = find_diagnostic_rules(identifier)
    if not rule:
        return {}
    record: dict[str, Any] = {
        "source": "compound_profiles",
        "actual_source": "compound_profiles",
        "query": identifier,
        "profile_class": rule.get("class"),
        "MolecularFormula": rule.get("formula"),
        "MolecularWeight": rule.get("mass"),
        "pathway": rule.get("pathway", ""),
        "known_metabolites": rule.get("known_metabolites", []),
        "diagnostic_fragments": rule.get("fragments", []),
    }
    return {key: value for key, value in record.items()
            if value not in (None, "", [], {})}


@register_tool
class PubChemCompoundPropertiesTool(_RemoteTool):
    tool_name = "pubchem_compound_properties"
    version = "1.0.0"
    description = "PubChem 分子式/分子量/CAS/SMILES 查询，含 30 分钟熔断降级"
    supported_data_types = ["compound", "metabolite"]
    parameters = [
        ToolParameter("name", "string", "Compound name", required=False),
        ToolParameter("cid", "string", "PubChem CID", required=False),
        ToolParameter("properties", "array", "Property names; CAS is resolved through PubChem synonyms", default=["MolecularFormula", "MolecularWeight", "IsomericSMILES", "InChIKey"]),
        ToolParameter("names", "array", "Batch of compound names", required=False),
    ]

    def validate_input(self, **kwargs) -> list[str]:
        if kwargs.get("names"):
            return []
        return [] if kwargs.get("name") or kwargs.get("cid") else ["需要 name 或 cid"]

    def run(self, **kwargs) -> ToolResult:
        names = [str(n).strip() for n in (kwargs.get("names") or []) if str(n).strip()]
        if names:
            properties = kwargs.get("properties") or ["MolecularFormula", "MolecularWeight", "IsomericSMILES", "InChIKey", "CAS"]
            rows: list[dict[str, Any]] = []
            warnings: list[str] = []
            for name in names:
                item = self.run(name=name, properties=properties)
                record = (item.metadata.get("records") or [{}])[0] if item.metadata.get("records") else {}
                rows.append({"name": name, "properties": record,
                             "source": record.get("source", item.metadata.get("source", "unresolved")),
                             "actual_source": record.get("actual_source", "unresolved"),
                             "warnings": list(item.warnings)})
                warnings.extend([f"{name}: {w}" for w in item.warnings])
            missing = sum(1 for row in rows if not row["properties"])
            return ToolResult(
                status="partial" if warnings or missing else "success",
                evidence_list=[], warnings=warnings,
                metadata={"source": "PubChem batch", "count": len(rows),
                          "resolved": len(rows) - missing, "missing": missing,
                          "records": rows, "properties": properties,
                          "actual_source": "PubChem/compound_profiles"},
            )
        global _PUBCHEM_FAILURES, _PUBCHEM_OPEN_UNTIL, _PUBCHEM_LAST_PROBE
        identifier = str(kwargs.get("cid") or kwargs.get("name") or "").strip()
        namespace = "cid" if kwargs.get("cid") else "name"
        props = kwargs.get("properties") or ["MolecularFormula", "MolecularWeight", "IsomericSMILES", "InChIKey"]
        key = "pchem:" + hashlib.md5(json.dumps([namespace, identifier, props], sort_keys=True).encode()).hexdigest()
        cached = knowledge_cache.get(key)
        if cached is not None:
            cached.setdefault("source", "PubChem")
            cached.setdefault("actual_source", "PubChem")
            return self._evidence("PubChem", identifier, 1, [cached])
        now = time.time()
        with _PUBCHEM_LOCK:
            open_until, last_probe = _PUBCHEM_OPEN_UNTIL, _PUBCHEM_LAST_PROBE
            circuit_open = now < open_until
            probe_due = now - last_probe >= 600
            if circuit_open and not probe_due:
                local = _local_properties(identifier)
                warning = "PubChem 当前不可达（熔断中），已使用 compound_profiles 本地规则/缓存；每 10 分钟探测一次"
                if not local:
                    warning += "；compound_profiles 未收录该化合物"
                return self._evidence("PubChem（本地降级）", identifier, 1 if local else 0, [local] if local else [], warning)
            if circuit_open and probe_due:
                _PUBCHEM_LAST_PROBE = now
        url = f"{PUBCHEM_URL}/{namespace}/{requests.utils.quote(identifier, safe='')}/property/{','.join(props)}/JSON"
        try:
            response = requests.get(url, timeout=15)
            response.raise_for_status()
            data = response.json()
            rows = data.get("PropertyTable", {}).get("Properties", []) if isinstance(data, dict) else []
            record = dict(rows[0]) if rows else {}
            if record:
                record["source"] = "PubChem"
                record["actual_source"] = "PubChem"
            if record and "CAS" in props:
                try:
                    syn_url = f"{PUBCHEM_URL}/{namespace}/{requests.utils.quote(identifier, safe='')}/synonyms/JSON"
                    syn_response = requests.get(syn_url, timeout=15)
                    syn_response.raise_for_status()
                    syn_data = syn_response.json()
                    synonyms = (((syn_data.get("InformationList") or {}).get("Information") or [{}])[0].get("Synonym") or [])
                    import re
                    cas = next((str(s) for s in synonyms if re.fullmatch(r"\d{2,7}-\d{2}-\d", str(s))), "")
                    if cas:
                        record["CAS"] = cas
                        record["CAS (candidate, from synonyms)"] = cas
                except Exception as syn_exc:
                    record["CAS_warning"] = f"CAS synonym lookup failed: {_error_message(syn_exc)}"
            knowledge_cache.set(key, record, ttl=30 * 86400)
            with _PUBCHEM_LOCK:
                _PUBCHEM_FAILURES = 0
                _PUBCHEM_OPEN_UNTIL = 0.0
            return self._evidence("PubChem", identifier, 1 if rows else 0, [record] if rows else [])
        except Exception as exc:
            with _PUBCHEM_LOCK:
                _PUBCHEM_FAILURES += 1
                if _PUBCHEM_FAILURES >= 2:
                    _PUBCHEM_OPEN_UNTIL = time.time() + 1800
                    _PUBCHEM_LAST_PROBE = time.time()
            local = _local_properties(identifier)
            warning = f"PubChem 请求失败，已降级到 compound_profiles 本地规则: {_error_message(exc)}"
            if not local:
                warning += "；compound_profiles 未收录该化合物"
            return self._evidence("PubChem（本地降级）", identifier, 1 if local else 0, [local] if local else [], warning)


def remote_annotation_fallback(query: str, *, precursor_mz: float | None = None, peaks: Iterable[Any] | None = None) -> ToolResult:
    """Execute the declared MassBank -> MoNA -> PubChem -> local chain."""
    warnings: list[str] = []
    mb = MassBankSpectrumSearchTool().run(compound_name=query, precursor_mz=precursor_mz, peaks=list(peaks or []), limit=10)
    if mb.metadata.get("count", 0):
        return mb
    warnings.extend(mb.warnings)
    mn = MonaSpectrumSearchTool().run(query=f"compound.metaData.value~'*{query}*'", size=10, page=0)
    if mn.metadata.get("count", 0):
        mn.warnings = warnings + mn.warnings
        return mn
    warnings.extend(mn.warnings)
    pc = PubChemCompoundPropertiesTool().run(name=query)
    warnings.extend(pc.warnings)
    # Keep the remote property result, then perform the actual local MS2 rule
    # fallback when a spectrum was supplied.
    local_summary: dict[str, Any] = {}
    if precursor_mz is not None and peaks:
        try:
            from phyto_reason.metabolomics.annotator import annotate_spectrum
            from phyto_reason.metabolomics.spectrum_io import Spectrum
            local = annotate_spectrum(Spectrum(precursor_mz=float(precursor_mz), fragments=[
                (float(p[0]), float(p[1])) if isinstance(p, (list, tuple)) else (float(p), 1.0)
                for p in peaks
            ]), target_metabolite=query)
            local_summary = {
                "formula_candidates": local.formula_candidates[:8],
                "candidates": [c.__dict__ for c in local.candidates[:8]],
                "note": local.note,
            }
        except Exception as exc:
            warnings.append(f"本地 compound_profiles 注释失败: {exc}")
    pc.status = "partial"
    pc.metadata["local_compound_profiles"] = local_summary
    pc.warnings = warnings + ["已降级到本地 compound_profiles 规则；远程来源未命中"]
    return pc
