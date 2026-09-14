"""sequences.py — FASTA 读取与 ORF 翻译（TF 预测的序列前处理）。

纯 Python 标准密码子表实现，无外部依赖：
- 读取 FASTA 为 {id: sequence}
- 自动判别核酸 / 蛋白
- 核酸输入按 6 阅读框找最长 ORF 翻译（优先 ATG 起始）
"""

from __future__ import annotations

from pathlib import Path

_STOP = {"TAA", "TAG", "TGA"}
_DNA_ALPHABET = set("ACGTUNRYKMSWBDHVacgtunrykmswbdhv")

_CODON = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L", "TCT": "S", "TCC": "S",
    "TCA": "S", "TCG": "S", "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W", "CTT": "L", "CTC": "L",
    "CTA": "L", "CTG": "L", "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q", "CGT": "R", "CGC": "R",
    "CGA": "R", "CGG": "R", "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T", "AAT": "N", "AAC": "N",
    "AAA": "K", "AAG": "K", "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V", "GCT": "A", "GCC": "A",
    "GCA": "A", "GCG": "A", "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

_COMPLEMENT = str.maketrans("ACGTUNacgtun", "TGCAANtgcaan")


def read_fasta(path: str | Path) -> dict[str, str]:
    """读取 FASTA 为 {id: sequence}；重复 id 自动加后缀保留全部记录。"""
    sequences: dict[str, str] = {}
    current: str | None = None
    chunks: list[str] = []
    with Path(path).open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                if current is not None:
                    _store(sequences, current, "".join(chunks))
                header = line[1:].split()[0] if len(line) > 1 else ""
                current = header or f"seq_{len(sequences) + 1}"
                chunks = []
            elif current is not None:
                chunks.append(line.upper())
    if current is not None:
        _store(sequences, current, "".join(chunks))
    return sequences


def _store(sequences: dict[str, str], key: str, value: str) -> None:
    if key in sequences:
        suffix = 2
        while f"{key}_{suffix}" in sequences:
            suffix += 1
        key = f"{key}_{suffix}"
    sequences[key] = value


def is_nucleotide(sequence: str) -> bool:
    if not sequence:
        return False
    dna = sum(1 for ch in sequence[:500] if ch in _DNA_ALPHABET)
    return dna / min(len(sequence), 500) > 0.9


def reverse_complement(sequence: str) -> str:
    return sequence.translate(_COMPLEMENT)[::-1]


def _translate_frame(dna: str, offset: int) -> list[tuple[int, str]]:
    """单帧翻译：返回 (起始核苷酸位置, 肽段) 列表，肽段按终止密码子切分。"""
    peptides: list[tuple[int, str]] = []
    amino: list[str] = []
    start = offset
    for i in range(offset, len(dna) - 2, 3):
        codon = dna[i:i + 3]
        residue = _CODON.get(codon, "X")
        if residue == "*":
            if len(amino) >= 1:
                peptides.append((start, "".join(amino)))
            amino = []
            start = i + 3
        else:
            if not amino:
                start = i
            amino.append(residue)
    if len(amino) >= 1:
        peptides.append((start, "".join(amino)))
    return peptides


def longest_orf_protein(dna: str, min_aa: int = 30) -> str | None:
    """6 阅读框最长 ORF 翻译。

    优先取最长的 ATG 起始完整 ORF（≥min_aa）；无符合者退化为
    最长的终止密码子间片段（de novo 转录本 5' 常不完整）。
    """
    dna = dna.upper().replace("U", "T")
    frames = [dna, dna[1:], dna[2:], reverse_complement(dna),
              reverse_complement(dna)[1:], reverse_complement(dna)[2:]]
    best_complete: tuple[int, str] | None = None  # (aa_len, protein)
    best_open: tuple[int, str] | None = None
    for frame in frames:
        for _nt_pos, peptide in _translate_frame(frame, 0):
            if len(peptide) < min_aa:
                continue
            if best_open is None or len(peptide) > len(best_open[1]):
                best_open = (len(peptide), peptide)
            atg_index = peptide.find("M")
            if atg_index >= 0 and len(peptide) - atg_index >= min_aa:
                complete = peptide[atg_index:]
                if best_complete is None or len(complete) > len(best_complete[1]):
                    best_complete = (len(complete), complete)
    if best_complete is not None:
        return best_complete[1]
    if best_open is not None:
        return best_open[1]
    return None


def protein_sequences_from_fasta(path: str | Path, min_aa: int = 30) -> tuple[dict[str, str], str]:
    """FASTA → 蛋白序列集合。核酸输入自动六框翻译，返回 (序列, 模式 aa|nt)。"""
    raw = read_fasta(path)
    if not raw:
        raise ValueError(f"FASTA 为空或无法解析: {path}")
    sample = next(iter(raw.values()))
    if is_nucleotide(sample):
        proteins: dict[str, str] = {}
        missing = 0
        for gene_id, dna in raw.items():
            protein = longest_orf_protein(dna, min_aa=min_aa)
            if protein:
                proteins[gene_id] = protein
            else:
                missing += 1
        if not proteins:
            raise ValueError("核酸 FASTA 中未找到任何 ≥min_aa 的 ORF；请确认输入是转录本/CDS 序列")
        return proteins, "nt"
    return dict(raw), "aa"
