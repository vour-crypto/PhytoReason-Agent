# annotation_mapper.py
from pathlib import Path
import pandas as pd
class AnnotationIndex:
    def __init__(self):
        self.term_to_genes = {}
        self.gene_to_terms = {}
        self.n_genes = 0
        self.warnings = []
    def search(self, keyword):
        kw = keyword.lower()
        results = set()
        for term, genes in self.term_to_genes.items():
            if kw in term.lower():
                results.update(genes)
        return results

    def search_substring(self, keyword):
        kw = keyword.lower()
        results = set()
        for gid, terms in self.gene_to_terms.items():
            full = chr(32).join(terms).lower()
            if kw in full:
                results.add(gid)
        return results
    def get_terms(self, gene_id):
        return self.gene_to_terms.get(gene_id, [])


def load_annotation(path, text_columns=None):
    path = Path(path)
    idx = AnnotationIndex()
    ext = path.suffix.lower()
    if ext == chr(46) + chr(120) + chr(108) + chr(115) + chr(120):
        df = pd.read_excel(path)
    else:
        sep = chr(9) if ext == chr(46) + chr(116) + chr(115) + chr(118) else chr(44)
        df = pd.read_csv(path, sep=sep, low_memory=False)
    id_col = None
    for col in df.columns:
        cl = col.lower().strip()
        if cl in (chr(103)+chr(101)+chr(110)+chr(101)+chr(105)+chr(100), chr(103)+chr(101)+chr(110)+chr(101)+chr(95)+chr(105)+chr(100)):
            id_col = col
            break
    if id_col is None:
        idx.warnings.append(chr(78)+chr(111)+chr(32)+chr(103)+chr(101)+chr(110)+chr(101)+chr(32)+chr(73)+chr(68)+chr(32)+chr(99)+chr(111)+chr(108))
        return idx
    if text_columns is None:
        desc = ["kegg","annotation","description","function","pfam","swissprot","trembl","nr","go","kog","eggnog"]
        text_columns = [c for c in df.columns if any(k in c.lower() for k in desc)]
    for idx_r, row in df.iterrows():
        gid = str(row[id_col]).strip()
        if not gid or gid in (chr(110)+chr(97)+chr(110), chr(34)+chr(34)): continue
        terms = []
        for col in text_columns:
            if col in df.columns:
                v = str(row[col]).strip()
                if v and v not in (chr(110)+chr(97)+chr(110), chr(34)+chr(34), chr(45)+chr(45), chr(45)): terms.append(v)
        if not terms: continue
        idx.n_genes += 1
        if gid not in idx.gene_to_terms:
            idx.gene_to_terms[gid] = []
        idx.gene_to_terms[gid].extend(terms)
        for term in terms:
            for token in term.replace(chr(124),chr(32)).replace(chr(59),chr(32)).split():
                t = token.strip()
                if len(t) >= 2:
                    idx.term_to_genes.setdefault(t, []).append(gid)
    idx.warnings.append(chr(76)+chr(111)+chr(97)+chr(100)+chr(101)+chr(100)+chr(32)+str(idx.n_genes)+chr(32)+chr(97)+chr(110)+chr(110)+chr(111)+chr(116)+chr(97)+chr(116)+chr(105)+chr(111)+chr(110)+chr(115))
    return idx