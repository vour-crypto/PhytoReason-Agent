"""Build the PlantTFDB local TF index (incremental or full refresh)."""
from __future__ import annotations

import argparse
import gzip
import json
import re
import sqlite3
import sys
import time
from datetime import date
from pathlib import Path

import requests

DB_PATH = Path(__file__).resolve().parents[1] / "data" / "planttf.db"
REPORT_PATH = DB_PATH.with_name("build_report.txt")
SPECIES_CACHE = DB_PATH.with_name("species_cache.json")
_UA = {"User-Agent": "PhytoReason/5.1"}
_BATCH = 1000
_SLEEP = 0.5
_FAMILY_HINTS = ("NAC", "MYB", "WRKY", "bZIP", "bHLH", "ERF", "HB", "TCP", "GRAS", "C2H2", "AP2", "LOB", "HD-ZIP", "SBP")
_FALLBACK = ["Os", "Ath", "Zm", "Sl", "Gm", "Ta"]


def discover_all_species() -> dict[str, str]:
    if SPECIES_CACHE.exists():
        cached = json.loads(SPECIES_CACHE.read_text(encoding="utf-8"))
        if cached.get("codes"):
            return cached["codes"]
    codes: dict[str, str] = {}
    for base in ("https://planttfdb.gao-lab.org", "https://plantregmap.gao-lab.org"):
        try:
            response = requests.get(f"{base}/download.php", timeout=15, headers=_UA)
            if not response.ok:
                continue
            for href in re.findall(r'href="([^"]+)"', response.text):
                match = re.search(r"(?:[?&]sp=|/)([A-Z][a-z]{1,2})(?:\.|/|$|_)", href)
                if match and len(match.group(1)) in (2, 3):
                    codes.setdefault(match.group(1), match.group(1))
            if codes:
                SPECIES_CACHE.parent.mkdir(parents=True, exist_ok=True)
                SPECIES_CACHE.write_text(json.dumps({"codes": codes, "source": base, "date": date.today().isoformat()}, indent=2), encoding="utf-8")
                return codes
        except requests.RequestException:
            continue
    return {code: code for code in _FALLBACK}


def discover_url(code: str) -> str | None:
    aliases = {"Os": "Osj", "Zm": "Zma", "Sl": "Sly", "Gm": "Gma", "Ta": "Tae"}
    actual = aliases.get(code, code)
    for url in (f"https://planttfdb.gao-lab.org/download/TF_list/{actual}_TF_list.txt.gz", f"https://plantregmap.gao-lab.org/download/TF_list/{actual}_TF_list.txt.gz"):
        try:
            response = requests.get(url, timeout=15, headers=_UA)
            if response.ok and len(response.content) > 100:
                return url
        except requests.RequestException:
            continue
    return None


def parse_lines(text: str, code: str) -> list[tuple[str, str, str, str]]:
    rows = []
    header: dict[str, int] = {}
    for raw in text.splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t") if "\t" in line else line.split(",")
        if len(parts) < 2:
            continue
        lowered = [part.strip().lower() for part in parts]
        if any("gene" in part or "family" in part for part in lowered):
            header = {part: index for index, part in enumerate(lowered)}
            continue
        if header:
            gene_index = header.get("gene_id", header.get("gene", 1))
            family_index = header.get("family", 2)
            gene_id = parts[gene_index].strip() if gene_index < len(parts) else parts[0].strip()
            family = parts[family_index].strip() if family_index < len(parts) else ""
            rows.append((code, gene_id, family, ""))
            continue
        first, second = parts[0].strip(), parts[1].strip()
        first_family = any(hint in first.upper() for hint in _FAMILY_HINTS)
        second_family = any(hint in second.upper() for hint in _FAMILY_HINTS)
        gene_id, family = (second, first) if first_family and not second_family else (first, second)
        rows.append((code, gene_id, family, parts[2].strip() if len(parts) > 2 else ""))
    return rows


def valid(rows: list[tuple[str, str, str, str]]) -> bool:
    families = {row[2].upper() for row in rows}
    return len(rows) >= 10 and any(any(hint in family for family in families) for hint in _FAMILY_HINTS)


def upsert_species(conn: sqlite3.Connection, code: str, rows: list[tuple[str, str, str, str]]) -> None:
    conn.execute("DELETE FROM tf WHERE species=?", (code,))
    for offset in range(0, len(rows), _BATCH):
        conn.executemany("INSERT INTO tf VALUES(?,?,?,?)", rows[offset:offset + _BATCH])
    conn.commit()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--all", action="store_true")
    parser.add_argument("--species", default="")
    args = parser.parse_args()
    if args.all:
        codes = discover_all_species()
    elif args.species:
        codes = {code.strip(): code.strip() for code in args.species.split(",") if code.strip()}
    else:
        print("用法：--all 或 --species Os,Ath")
        return 1
    DB_PATH.parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(DB_PATH)
    conn.execute("CREATE TABLE IF NOT EXISTS tf (species TEXT, gene_id TEXT, family TEXT, description TEXT)")
    conn.execute("CREATE INDEX IF NOT EXISTS idx_tf_gene ON tf(gene_id)")
    conn.execute("CREATE INDEX IF NOT EXISTS idx_tf_family ON tf(family)")
    conn.execute("CREATE INDEX IF NOT EXISTS idx_tf_description ON tf(description)")
    conn.execute("CREATE INDEX IF NOT EXISTS idx_tf_gene_family_desc ON tf(gene_id,family,description)")
    report = ["# PlantTFDB 本地索引构建报告", f"# 构建日期: {date.today().isoformat()}", f"# 模式: {'全库' if args.all else '增量'} · 物种数 {len(codes)}"]
    ok = fail = 0
    for code, label in codes.items():
        url = discover_url(code)
        if not url:
            report.append(f"[FAIL] {label}({code}): 未找到 TF 列表 URL"); fail += 1; continue
        try:
            response = requests.get(url, timeout=30, headers=_UA); response.raise_for_status()
            payload = gzip.decompress(response.content) if url.endswith(".gz") else response.content
            rows = parse_lines(payload.decode("utf-8", errors="replace"), code)
        except requests.RequestException as exc:
            report.append(f"[FAIL] {label}({code}): 下载失败 {exc} · {url}"); fail += 1; continue
        if not valid(rows):
            report.append(f"[FAIL] {label}({code}): 解析异常 · {url}"); fail += 1; continue
        upsert_species(conn, code, rows); report.append(f"[OK] {label}({code}): {len(rows)} TF · {len({row[2] for row in rows})} 家族 · {url}"); ok += 1; time.sleep(_SLEEP)
    conn.close(); report.append(f"# 合计: 成功 {ok} · 失败 {fail}"); REPORT_PATH.write_text("\n".join(report), encoding="utf-8"); print("\n".join(report)); return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
