"""
upload.py — POST /upload 文件上传端点 (v4.0).

支持多文件上传到同一 session，数据存入 AgentOrchestrator 供 chat 使用。
"""

from __future__ import annotations

import tempfile
import uuid
import os
from pathlib import Path

from fastapi import APIRouter, UploadFile, File, Form, HTTPException

from phyto_reason.agent.orchestrator import AgentOrchestrator

router = APIRouter(tags=["upload"])

MAX_UPLOAD_BYTES = 100 * 1024 * 1024

_orchestrator = AgentOrchestrator()


def _parse_fasta_text(text: str) -> dict:
    """Parse FASTA text into {id: sequence} dict."""
    result = {}
    current_id, current_seq = None, []
    for line in text.split("\n"):
        line = line.strip()
        if line.startswith(">"):
            if current_id:
                result[current_id] = "".join(current_seq)
            current_id = line[1:].split()[0]
            current_seq = []
        elif current_id:
            current_seq.append(line)
    if current_id:
        result[current_id] = "".join(current_seq)
    return result


def _header_suggests_long_table(content: bytes, encoding: str) -> bool:
    """从首行表头粗判是否为 m/z-RT-intensity 长表（权威判定在解析器内）。"""
    header_line = content.split(b"\n", 1)[0].decode(encoding, errors="replace").lower()
    has_mz = "m/z" in header_line or "mz" in header_line
    has_rt = "rt" in header_line or "retention" in header_line
    has_intensity = "intensity" in header_line or "area" in header_line or "height" in header_line
    has_sample = "sample" in header_line
    return has_mz and has_rt and has_intensity and has_sample


@router.post("/upload")
async def upload_file(
    file: UploadFile = File(...),
    metadata_file: UploadFile | None = File(default=None),
    session_id: str = Form(default=""),
    species: str = Form(default=""),
    target_metabolite: str = Form(default=""),
    data_type: str = Form(default=""),
):
    """上传分析文件并存入 session。

    支持: CSV/TSV/XLSX 表达矩阵、代谢物矩阵（宽表或 m/z-RT-intensity 标准长表）、
    FASTA 启动子序列。多次上传同一 session_id 可累积多份数据。
    data_type 为显式类型（expression/metabolite），优先于文件名推断。
    """
    content = await file.read(MAX_UPLOAD_BYTES + 1)
    if len(content) > MAX_UPLOAD_BYTES:
        raise HTTPException(
            status_code=413,
            detail=f"Uploaded file exceeds the {MAX_UPLOAD_BYTES // (1024 * 1024)} MB limit.",
        )
    filename = file.filename or "unknown"
    name_lower = filename.lower()

    # ── Encoding detection ─────────────────────────────────
    # Try UTF-8 first, then common Chinese encodings, then Latin-1.
    _ENCODINGS_TO_TRY = ["utf-8", "gbk", "gb18030", "gb2312", "latin-1"]
    detected_encoding = "utf-8"
    for enc in _ENCODINGS_TO_TRY:
        try:
            content.decode(enc)
            detected_encoding = enc
            break
        except (UnicodeDecodeError, LookupError):
            continue

    import logging
    logging.getLogger("api.upload").info(
        "File %s: %d bytes, detected encoding=%s", filename, len(content), detected_encoding
    )

    # Create or reuse session
    sid = session_id or f"api-{uuid.uuid4().hex[:12]}"

    # Save raw bytes to temp file (preserve original encoding)
    suffix = os.path.splitext(filename)[1] or ".csv"
    with tempfile.NamedTemporaryFile(delete=False, suffix=suffix) as tmp:
        tmp.write(content)
        tmp_path = tmp.name

    summary = {"filename": filename, "size_bytes": len(content), "file_type": "unknown"}

    # Prepare data dicts for orchestrator
    expr_dict = None
    meta_dict = None
    prom_dict = None
    sample_metadata = None
    metadata_tmp_path = None
    pending_tf_csv = ""
    ms2_csv = ""

    # ── File type detection (supports Chinese + English filenames) ──
    # Chinese keywords:
    #   代谢组/代谢物/代谢 = metabolome/metabolite
    #   转录组/表达/基因 = transcriptome/expression
    #   启动子 = promoter
    _CN_METABOLITE = ["代谢组", "代谢物", "代谢", "质谱", "液质", "气质", "色谱"]
    _CN_EXPRESSION = ["转录组", "表达", "基因", "测序", "rna", "seq"]
    _CN_PROMOTER = ["启动子", "promoter"]

    _is_metabolite = any(
        x in name_lower for x in ["metabolite", "metab", "meta", "metabolom"]
    ) or any(x in filename for x in _CN_METABOLITE)

    _is_promoter = any(x in name_lower for x in [".fasta", ".fa", ".fna", "promoter"])

    # 显式 data_type 优先（Phase 6.4 Step 2：槽位意图不被文件名覆盖）
    explicit_type = (data_type or "").strip().lower()
    if explicit_type in {"expression", "metabolite"}:
        _is_metabolite = explicit_type == "metabolite"
        _is_promoter = False
    elif not _is_metabolite and not _is_promoter and _header_suggests_long_table(content, detected_encoding):
        # 文件名看不出类型时，按表头嗅探标准长表（sample/m/z/RT/intensity）
        _is_metabolite = True
        summary["detection"] = "long_table_header_sniff"

    try:
        if explicit_type == "tf_fasta":
            # 转录本/蛋白 FASTA → TF 预测（Pfam HMM）→ 会话注释注入
            from phyto_reason.tf_prediction.predictor import predict_from_fasta

            summary["file_type"] = "tf_prediction"
            result = predict_from_fasta(tmp_path)
            csv_source = Path(result["annotation_csv"])
            csv_dest = Path(tempfile.gettempdir()) / f"tf_prediction_{sid}.csv"
            csv_source.replace(csv_dest)
            pending_tf_csv = str(csv_dest)
            summary.update({
                "n_sequences": result["n_sequences"],
                "n_predicted_tf": result["n_predicted_tf"],
                "family_counts": result["family_counts"],
                "annotation_csv": str(csv_dest),
                "mode": result["mode"],
            })
        elif explicit_type == "ms2":
            # MS/MS 谱图（MGF/CSV/MSP）→ 注释（分子式/结构类/MSI 分级 + 物种锚定）
            import csv as _csv

            from phyto_reason.metabolomics.annotator import annotate_spectrum
            from phyto_reason.metabolomics.spectrum_io import load_spectrum

            summary["file_type"] = "ms2_spectra"
            session_for_species = _orchestrator._store.get(sid)
            species_name = (
                species
                or (getattr(session_for_species, "species", "") if session_for_species else "")
            )
            spectra = load_spectrum(tmp_path, polarity="")
            if not spectra:
                summary["parse_error"] = "未从文件解析到任何谱图（支持 MGF/CSV/MSP）"
            else:
                limit = min(len(spectra), 50)
                rows = []
                msi_counts: dict[str, int] = {}
                for i, spec in enumerate(spectra[:limit]):
                    res = annotate_spectrum(spec, species=species_name)
                    best = res.candidates[0] if res.candidates else None
                    msi = best.msi_level if best else 4
                    msi_counts[f"Level {msi}"] = msi_counts.get(f"Level {msi}", 0) + 1
                    rows.append({
                        "scan": spec.scan_id or str(i + 1),
                        "precursor_mz": round(spec.precursor_mz, 4),
                        "adduct": res.adduct or spec.adduct,
                        "neutral_mass": round(res.neutral_mass, 4) if res.neutral_mass else "",
                        "msi_level": msi,
                        "best_class": best.class_name if best else "",
                        "best_name": (best.name or "") if best else "",
                        "formula": (best.formula or (res.formula_candidates[0].get("formula", "") if res.formula_candidates else "")),
                        "score": round(best.score, 1) if best else "",
                        "profile_hit": (res.profile_hits[0].get("name", "") if res.profile_hits else ""),
                    })
                summary["n_spectra"] = len(spectra)
                summary["n_annotated"] = len(rows)
                summary["msi_counts"] = msi_counts
                summary["species_anchored"] = species_name
                if len(spectra) > limit:
                    summary["skipped_spectra"] = len(spectra) - limit
                from phyto_reason.platform_paths import output_dir
                out_csv = output_dir() / f"ms2_annotation_{sid}.csv"
                with out_csv.open("w", encoding="utf-8-sig", newline="") as handle:
                    writer = _csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
                    writer.writeheader()
                    writer.writerows(rows)
                summary["annotation_csv"] = str(out_csv)

        elif _is_promoter or any(x in filename for x in _CN_PROMOTER):
            summary["file_type"] = "promoter_fasta"
            prom_dict = _parse_fasta_text(content.decode(detected_encoding))
            summary["n_sequences"] = len(prom_dict)

        elif _is_metabolite:
            from phyto_reason.ingestion.parsers.metabolite_parser import parse_metabolite
            parsed = parse_metabolite(tmp_path, encoding=detected_encoding)
            meta_dict = parsed.to_workflow_dict()
            summary["file_type"] = "metabolite_matrix"
            summary["n_metabolites"] = parsed.n_features
            summary["n_samples"] = parsed.n_samples
            summary["provenance"] = parsed.provenance
            # Include classification details so user can verify column detection
            if parsed.classification:
                summary["sample_columns"] = parsed.classification.get("sample_columns", [])
                summary["excluded_columns"] = parsed.classification.get("excluded_metadata", [])
                summary["content_resolved"] = parsed.classification.get("content_resolved", [])
            if parsed.warnings:
                summary["warnings"] = parsed.warnings

        else:
            # Expression matrix (default) OR auto-detect by content
            # Check if content looks like metabolite data (fewer rows, compound names)
            from phyto_reason.ingestion.parsers.expression_parser import parse_expression
            parsed = parse_expression(tmp_path, encoding=detected_encoding)
            expr_dict = parsed.to_workflow_dict()
            summary["file_type"] = "expression_matrix"
            summary["n_genes"] = parsed.n_features
            summary["n_samples"] = parsed.n_samples
            summary["provenance"] = parsed.provenance

        # Optional sample design/metadata table. It is parsed independently so
        # expression and metabolite matrices can be uploaded in either order.
        if metadata_file is not None:
            metadata_content = await metadata_file.read(MAX_UPLOAD_BYTES + 1)
            if len(metadata_content) > MAX_UPLOAD_BYTES:
                raise HTTPException(status_code=413, detail="Metadata file exceeds the upload limit.")
            metadata_name = metadata_file.filename or "metadata.csv"
            metadata_suffix = os.path.splitext(metadata_name)[1] or ".csv"
            with tempfile.NamedTemporaryFile(delete=False, suffix=metadata_suffix) as metadata_tmp:
                metadata_tmp.write(metadata_content)
                metadata_tmp_path = metadata_tmp.name
            from phyto_reason.ingestion.parsers.metadata_parser import parse_metadata
            parsed_md = parse_metadata(metadata_tmp_path)
            sample_metadata = {
                "source_file": metadata_name,
                "sample_ids": parsed_md.sample_ids,
                "columns_present": parsed_md.columns_present,
                "columns_missing": parsed_md.columns_missing,
                "samples": [sample.model_dump() for sample in parsed_md.samples],
            }
            summary["metadata_file"] = metadata_name
            summary["metadata_warnings"] = parsed_md.warnings

        # ── Anti-overwrite: merge with existing session data ──
        existing = _orchestrator._store.get(sid)
        if existing and existing.has_data:
            if expr_dict is not None and existing.expression_matrix:
                # If session already has expression data, treat new upload as metabolite
                # unless the existing expression is empty
                if existing.metabolite_matrix and not meta_dict:
                    # Both slots are already filled — stack into a combined key
                    import logging
                    logging.getLogger("api.upload").warning(
                        "Session %s: both expr and meta slots filled. "
                        "New expression data replaces old. Old: %d genes → New: %d genes",
                        sid, len(existing.expression_matrix), len(expr_dict),
                    )
            if meta_dict is not None and existing.metabolite_matrix and existing.expression_matrix is None:
                # Metabolite was uploaded before expression — swap
                import logging
                logging.getLogger("api.upload").info(
                    "Session %s: metabolite uploaded before expression. Order corrected.", sid
                )

        # Store into orchestrator session (bridge upload <-> chat)
        _orchestrator.upload_data(
            session_id=sid,
            expression=expr_dict,
            metabolite=meta_dict,
            promoter=prom_dict,
            sample_metadata=sample_metadata,
            species=species,
            target_metabolite=target_metabolite,
        )

        # 会话在 upload_data 内创建——TF 预测注释在会话存在后注入
        if explicit_type == "tf_fasta":
            session_obj = _orchestrator._store.get(sid)
            if session_obj is not None:
                session_obj.tf_annotation_csv = pending_tf_csv
        if explicit_type == "ms2":
            session_obj = _orchestrator._store.get(sid)
            if session_obj is not None:
                session_obj.ms2_annotation_csv = summary.get("annotation_csv", "")

    except Exception as e:
        summary["parse_error"] = str(e)
        import logging
        logging.getLogger("api.upload").error("Upload parse failed: %s — %s", filename, e)
    finally:
        try:
            os.unlink(tmp_path)
        except Exception:
            pass
        if metadata_tmp_path:
            try:
                os.unlink(metadata_tmp_path)
            except Exception:
                pass

    # ── Log what was stored ────────────────────────────────
    import logging
    logger_upload = logging.getLogger("api.upload")
    logger_upload.info(
        "Upload to session %s: file=%s type=%s expr=%s meta=%s prom=%s species=%s target=%s",
        sid, filename, summary.get("file_type"),
        "yes" if expr_dict else "no",
        "yes" if meta_dict else "no",
        "yes" if prom_dict else "no",
        species, target_metabolite,
    )

    # Get session info for response — includes live data stats
    session_info = _orchestrator.get_session_info(sid)

    # ── Diagnostic: verify data is actually stored ──────────
    verify_session = _orchestrator._store.get(sid)
    stored_expr = bool(verify_session.expression_matrix) if verify_session else False
    stored_meta = bool(verify_session.metabolite_matrix) if verify_session else False
    expr_size = len(verify_session.expression_matrix) if (verify_session and verify_session.expression_matrix) else 0
    meta_size = len(verify_session.metabolite_matrix) if (verify_session and verify_session.metabolite_matrix) else 0
    logger_upload.info(
        "Upload VERIFY session %s: has_expr=%s (n=%d) has_meta=%s (n=%d) species=%s",
        sid, stored_expr, expr_size, stored_meta, meta_size, species,
    )

    return {
        "session_id": sid,
        "filename": filename,
        "species": species,
        "target_metabolite": target_metabolite,
        "summary": summary,
        "session": session_info,
        # Include diagnostic fields so the frontend can confirm data receipt
        "diagnostics": {
            "data_stored": stored_expr or stored_meta,
            "expression_genes": expr_size,
            "metabolite_compounds": meta_size,
            "species_stored": species,
            "parse_error": summary.get("parse_error", ""),
            "sample_alignment_report": (
                verify_session.sample_alignment_report if verify_session else {}
            ),
        },
    }
