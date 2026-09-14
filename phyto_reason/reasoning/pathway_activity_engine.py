# pathway_activity_engine.py
from __future__ import annotations
import logging, math, itertools
import numpy as np
from phyto_reason.reasoning.pathway_reasoner import PATHWAY_KNOWLEDGE, _detect_pathway_class
logger = logging.getLogger(__name__)

class PathwayActivityEngine:
    @staticmethod
    def compute_pathway_activity(metabolite, expression_matrix=None, annotation_index=None):
        result = {}
        result["pathway_name"] = "unknown"
        result["enzymes_found"] = []
        result["is_active"] = False
        if not expression_matrix:
            return result
        static_enzymes = []
        pw_key = _detect_pathway_class(metabolite)
        if pw_key and pw_key in PATHWAY_KNOWLEDGE:
            result["pathway_name"] = PATHWAY_KNOWLEDGE[pw_key]["name"]
            static_enzymes = PATHWAY_KNOWLEDGE[pw_key]["enzymes"]
        kegg = PathwayActivityEngine._fetch_kegg_enzymes(metabolite)
        if kegg:
            enzymes = kegg
        else:
            enzymes = static_enzymes
        if not enzymes:
            return result
        result["enzymes_total"] = len(enzymes)
        found = PathwayActivityEngine._search_enzymes(enzymes, expression_matrix, annotation_index)
        result["enzymes_found"] = found
        if len(found) < 2:
            return result
        PathwayActivityEngine._compute_scores(result, found, expression_matrix, annotation_index)
        return result
    @staticmethod
    def _fetch_kegg_enzymes(metabolite):
        import requests, time
        BASE = "https://rest.kegg.jp"
        for attempt in range(2):
            try:
                # Phase 6.4 Step 2: 移除关闭证书校验的写法，恢复 TLS 验证；
                # 网络失败降级为静态通路知识（PATHWAY_KNOWLEDGE），不静默关闭校验。
                r = requests.get(BASE+"/find/compound/"+metabolite, timeout=10)
                if r.status_code != 200: continue
                parts = r.text.strip().split(chr(10))
                if not parts: continue
                cpd = parts[0].split(chr(9))[0]
                if not cpd: continue
                time.sleep(0.3)
                d = requests.get(BASE+"/get/"+cpd, timeout=10).text
                pids = []
                for line in d.split(chr(10)):
                    if line.strip().startswith("PATHWAY"):
                        for p in line.split():
                            if p.startswith("map") or p.startswith("ko"): pids.append(p); break
                all_kos = set()
                for pid in pids[:2]:
                    time.sleep(0.3)
                    kdata = requests.get(BASE+"/link/ko/"+pid, timeout=10).text
                    for line in kdata.strip().split(chr(10)):
                        segs = line.split(chr(9))
                        if len(segs)>=2:
                            ko = segs[1].strip()
                            if ko.startswith("ko:"): all_kos.add(ko)
                if all_kos: return sorted(all_kos)
            except Exception as e:
                if attempt==0:
                    time.sleep(1)
                    continue
                logger.warning(
                    "KEGG lookup degraded to static pathway knowledge for %r after TLS/network failure: %s",
                    metabolite, e,
                )
        return []
    def _search_enzymes(enzyme_names, expression_matrix, annotation_index):
        found = []
        for e in enzyme_names:
            if not e or len(e)<2: continue
            if annotation_index is not None:
                ah = annotation_index.search_substring(e)
                if ah:
                    for g in ah:
                        if g in expression_matrix or g.rpartition(".")[0] in expression_matrix:
                            found.append(e); break
                    continue
            for g in expression_matrix:
                if e.upper() in g.upper():
                    found.append(e); break
        return found
    @staticmethod
    def _compute_scores(result, found_enzymes, expression_matrix, annotation_index):
        import math
        all_means = []
        for g, vals in expression_matrix.items():
            nums = [v for v in vals.values() if isinstance(v,(int,float)) and not math.isnan(v)]
            if nums: all_means.append(float(np.mean(nums)))
        gm = float(np.mean(all_means)) if all_means else 0.0
        gs = float(np.std(all_means)) if all_means else 1.0
        zscores = []
        for e in found_enzymes:
            ek = PathwayActivityEngine._gene_for_enzyme(e,expression_matrix,annotation_index)
            if ek:
                nums = [v for v in expression_matrix[ek].values() if isinstance(v,(int,float)) and not math.isnan(v)]
                if nums: zscores.append((float(np.mean(nums))-gm)/(gs+1e-10))
        if zscores: result["mean_z_score"] = round(float(np.mean(zscores)),3)
        pairs = []
        for e1,e2 in itertools.combinations(found_enzymes,2):
            g1 = PathwayActivityEngine._gene_for_enzyme(e1,expression_matrix,annotation_index)
            g2 = PathwayActivityEngine._gene_for_enzyme(e2,expression_matrix,annotation_index)
            if g1 and g2:
                v1 = [v for v in expression_matrix[g1].values() if isinstance(v,(int,float)) and not math.isnan(v)]
                v2 = [v for v in expression_matrix[g2].values() if isinstance(v,(int,float)) and not math.isnan(v)]
                if len(v1)>=3 and len(v2)>=3:
                    a1 = (np.array(v1)-np.mean(v1))/(np.std(v1)+1e-10)
                    a2 = (np.array(v2)-np.mean(v2))/(np.std(v2)+1e-10)
                    r = float(np.dot(a1,a2))/(len(a1)-1)
                    pairs.append(max(-1.0,min(1.0,r)))
        if pairs: result["coherence_score"] = round(max(0.0,float(np.mean(pairs))),3)
        coh = result.get("coherence_score",0.0)
        mz = result.get("mean_z_score",0.0)
        ne = len(found_enzymes)
        if coh >= 0.20:
            result["is_active"] = True
        result["activity_interpretation"] = " enz="+str(ne)+" coh="+str(round(coh,2))+" z="+str(round(mz,2))
    @staticmethod
    def _gene_for_enzyme(e, expression_matrix, annotation_index):
        for g in expression_matrix:
            if e.upper() in g.upper(): return g
        if annotation_index is not None:
            for hit in annotation_index.search_substring(e):
                if hit in expression_matrix: return hit
                c = hit.rpartition(".")[0]
                if c in expression_matrix: return c
        return None
