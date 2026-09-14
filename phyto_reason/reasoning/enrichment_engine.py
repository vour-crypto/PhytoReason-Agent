# enrichment_engine.py
from __future__ import annotations
import logging, math
from scipy.stats import fisher_exact
from phyto_reason.utils.stats_utils import compute_fdr
logger = logging.getLogger(__name__)

class EnrichmentEngine:
    @staticmethod
    def _tokenize_annotations(annotation_index):
        if annotation_index is None: return {}
        GO_TERMS = [
            "biological_process","molecular_function","cellular_component",
            "metabolic","biosynthetic","regulation","response",
            "oxidation","reduction","transferase","hydrolase",
            "transcription","phosphorylation","transport",
        ]
        return {g: [t for t in GO_TERMS if any(t in ann.lower() for ann in terms)]
                for g, terms in annotation_index.gene_to_terms.items()}

    @staticmethod
    def enrich(gene_list, background_list, annotation_index):
        gene_set = set(gene_list)
        bg_set = set(background_list)
        go_map = EnrichmentEngine._tokenize_annotations(annotation_index)
        results = []
        def _gk(g):
            if g in go_map: return g
            if g+".1" in go_map: return g+".1"
            c = g.rpartition(".")[0]
            return c if c in go_map else ""
        all_terms_seen = set()
        for g in (gene_set | bg_set):
            gk = _gk(g)
            if not gk: continue
            for t in go_map.get(gk, []):
                all_terms_seen.add(t)
        for term in sorted(all_terms_seen):
            a = sum(1 for g in gene_set if term in go_map.get(_gk(g), []))
            b = len(gene_set) - a
            c = sum(1 for g in bg_set if term in go_map.get(_gk(g), []))
            d = len(bg_set) - c
            if a < 2: continue
            try:
                _, p = fisher_exact([[a,b],[c,d]], alternative="greater")
            except: continue
            if not (math.isnan(p) or p < 1e-300):
                results.append((term, a, b, c, d, p))
        if not results: return []
        pvals = [r[5] for r in results]
        qvals = compute_fdr(pvals)
        enriched = []
        for (term, a, b, c, d, p), q in zip(results, qvals):
            if q < 0.05:
                enriched.append({
                    "term":term,"pvalue":p,"fdr":q,
                    "in_gene_set":a,"not_in_gene_set":b,
                    "in_background":c,"not_in_background":d,
                })
        return enriched