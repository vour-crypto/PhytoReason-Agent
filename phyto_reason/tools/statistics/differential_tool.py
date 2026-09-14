# differential_tool.py
from __future__ import annotations
import logging, math
import numpy as np
from scipy.stats import ttest_ind
from sklearn.decomposition import PCA
from phyto_reason.models.tool_result import ToolResult
from phyto_reason.models.candidate_gene import CandidateGene
from phyto_reason.models.evidence import Evidence,EvidenceType
from phyto_reason.tools.base_tool import BaseTool,ToolParameter
from phyto_reason.tools.tool_registry import register_tool
from phyto_reason.utils.stats_utils import compute_fdr
logger=logging.getLogger(__name__)

@register_tool
class DifferentialTool(BaseTool):
    tool_name='differential_analysis'
    def validate_input(self,**kwargs):return []
    description='Fold change t-test FDR PCA'
    parameters=[ToolParameter(name='expression_matrix',type='object',required=True),ToolParameter(name='group_a',type='array',required=True),ToolParameter(name='group_b',type='array',required=True)]
    supported_data_types=['expression','metabolite']

    def run(self,**kwargs):
        mat=kwargs.get('expression_matrix',{})
        ga=kwargs.get('group_a',[]);gb=kwargs.get('group_b',[])
        if len(ga)<2 or len(gb)<2:
            return ToolResult(candidates=[],evidence_list=[],warnings=['need >=2 per group'])
        results=[];all_p=[]
        candidates=[]
        for fid,svals in mat.items():
            a=np.array([svals.get(s,np.nan) for s in ga],dtype=float)
            b=np.array([svals.get(s,np.nan) for s in gb],dtype=float)
            a=a[~np.isnan(a)];b=b[~np.isnan(b)]
            if len(a)<2 or len(b)<2: continue
            fc=np.mean(b)/max(np.mean(a),1e-10)
            lfc=math.log2(fc) if fc>0 else 0
            _,p=ttest_ind(a,b,equal_var=False)
            if math.isnan(p) or p<1e-300: p=1e-300
            results.append((fid,lfc,p));all_p.append(p)
        if not results:
            return ToolResult(candidates=[],evidence_list=[],warnings=['no valid features'])
        qvals=compute_fdr(all_p)
        volcano=[];n_sig=0
        for (fid,lfc,p),q in zip(results,qvals):
            np10=-math.log10(max(p,1e-100))
            volcano.append({'gene':fid,'lfc':round(lfc,3),'p':round(p,6),'q':round(q,6),'np10':round(np10,3)})
            if q<0.05 and abs(lfc)>=1.0:
                n_sig+=1
                cg=CandidateGene(gene_id=fid,pathway_score=round(min(abs(lfc)/5,1),3))
                e=Evidence(evidence_type=EvidenceType.DEG,source='differential',score=round(min(abs(lfc)/5,1),3),confidence=round(1-min(q*10,0.95),3),description='log2FC='+str(round(lfc,2))+' q='+str(round(q,4)))
                candidates.append(cg)
        pca_data=None
        try:
            ss=list(set(ga+gb))
            X=[]
            for s in ss: X.append([mat[fid].get(s,0) for fid in list(mat.keys())[:500]])
            if len(X)>=3:
                pca=PCA(n_components=3);scores=pca.fit_transform(X)
                pca_data={'samples':ss,'pc1':[float(x[0]) for x in scores],'pc2':[float(x[1]) for x in scores],'var1':round(float(pca.explained_variance_ratio_[0])*100,1),'var2':round(float(pca.explained_variance_ratio_[1])*100,1)}
        except: pass
        return ToolResult(candidates=candidates,evidence_list=[],warnings=[],metadata={'n_sig':n_sig,'n_tested':len(results),'volcano':volcano,'pca':pca_data})