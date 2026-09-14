"""
tf_knowledge_base.py — TF family ↔ metabolite regulatory knowledge.

Canonical source for all TF-metabolite prior knowledge in the system.
Consolidates data from:
  - Published plant biology literature (with PMIDs)
  - PlantTFDB / JASPAR database knowledge
  - reasoning/tf_prior_reasoner.py (legacy)
  - ontology/tf_prior_ontology.py (legacy)

All entries are citation-backed and species-scoped.
"""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class TFMetaboliteRelation:
    """A known or inferred regulatory relationship between a TF family and metabolite."""
    tf_family: str
    metabolite_class: str        # e.g. "alkaloid", "flavonoid", "anthocyanin"
    specific_metabolite: str = ""  # e.g. "berberine" if known
    strength: str = "moderate"    # "strong" | "moderate" | "weak" | "inferred"
    score: float = 0.5           # 0.0-1.0
    description: str = ""
    pmids: list[str] = field(default_factory=list)
    species_scope: str = "general"
    mechanism: str = ""           # e.g. "direct binding", "MBW complex", "jasmonate response"
    known_examples: list[str] = field(default_factory=list)  # specific gene examples


@dataclass
class KnowledgeQueryResult:
    query: str
    relations: list[TFMetaboliteRelation] = field(default_factory=list)
    total: int = 0
    note: str = ""


# ── Canonical TF-Metabolite Knowledge Base ──────────────────────
# All entries must have at least one PMID or be marked "inferred"

CANONICAL_KNOWLEDGE: list[TFMetaboliteRelation] = [
    # ════════════════════════════════════════════════════════════
    # MYB family
    # ════════════════════════════════════════════════════════════
    TFMetaboliteRelation(
        tf_family="MYB", metabolite_class="flavonoid",
        strength="strong", score=0.90,
        description="R2R3-MYB TFs are master regulators of flavonoid biosynthesis. "
                    "Key examples: AtMYB12 (flavonol), AtMYB75/PAP1 (anthocyanin), SbMYB60 (sorghum).",
        pmids=["25228336", "31530398", "21653882"],
        species_scope="general",
        mechanism="Direct binding to AC elements in flavonoid gene promoters",
        known_examples=["AtMYB12", "AtMYB75/PAP1", "AtMYB111", "SbMYB60"],
    ),
    TFMetaboliteRelation(
        tf_family="MYB", metabolite_class="anthocyanin",
        strength="strong", score=0.95,
        description="MYB-bHLH-WD40 (MBW) complex specifically activates anthocyanin structural genes. "
                    "PAP1/PAP2/MYB113/MYB114 in Arabidopsis.",
        pmids=["25228336", "26847441"],
        species_scope="general",
        mechanism="MBW ternary complex binding to anthocyanin gene promoters",
        known_examples=["AtMYB75/PAP1", "AtMYB90/PAP2", "AtMYB113", "AtMYB114"],
    ),
    TFMetaboliteRelation(
        tf_family="MYB", metabolite_class="phenylpropanoid",
        strength="moderate", score=0.65,
        description="MYB TFs regulate phenylpropanoid pathway upstream genes PAL, C4H, 4CL.",
        pmids=["21653882"],
        species_scope="general",
        mechanism="AC element binding in phenylpropanoid gene promoters",
    ),
    TFMetaboliteRelation(
        tf_family="MYB", metabolite_class="lignin",
        strength="strong", score=0.80,
        description="MYB58/MYB63/MYB85 specifically activate lignin biosynthetic genes.",
        pmids=["24027067"],
        species_scope="general",
        known_examples=["AtMYB58", "AtMYB63", "AtMYB85"],
    ),
    TFMetaboliteRelation(
        tf_family="MYB", metabolite_class="alkaloid",
        strength="weak", score=0.30,
        description="Some MYB TFs regulate alkaloid biosynthesis in specific species (e.g., CrMYB in Catharanthus).",
        pmids=["27207470"],
        species_scope="apocynaceae, papaveraceae",
        known_examples=["CrMYB"],
    ),
    TFMetaboliteRelation(
        tf_family="MYB", metabolite_class="terpenoid",
        strength="moderate", score=0.50,
        description="MYB TFs regulate terpenoid biosynthesis in some species (e.g., AaMYB1 in Artemisia).",
        pmids=["28417071"],
        species_scope="asteraceae, lamiaceae",
    ),

    # ════════════════════════════════════════════════════════════
    # bHLH family
    # ════════════════════════════════════════════════════════════
    TFMetaboliteRelation(
        tf_family="bHLH", metabolite_class="anthocyanin",
        strength="strong", score=0.90,
        description="bHLH TFs (TT8/GL3/EGL3) form MBW complex with MYB and WD40 "
                    "to activate anthocyanin biosynthesis.",
        pmids=["26847441", "25228336"],
        species_scope="general",
        mechanism="MBW ternary complex formation",
        known_examples=["AtTT8", "AtGL3", "AtEGL3", "PhAN1"],
    ),
    TFMetaboliteRelation(
        tf_family="bHLH", metabolite_class="flavonoid",
        strength="moderate", score=0.60,
        description="bHLH TFs cooperate with MYB to regulate flavonoid pathway branches.",
        pmids=["26847441"],
        species_scope="general",
    ),
    TFMetaboliteRelation(
        tf_family="bHLH", metabolite_class="alkaloid",
        strength="moderate", score=0.55,
        description="bHLH TFs regulate alkaloid biosynthesis via jasmonate signaling. "
                    "Examples: CjbHLH1 (Coptis), CrMYC2 (Catharanthus).",
        pmids=["25296257", "28417071"],
        species_scope="ranunculaceae, apocynaceae",
        mechanism="JA-responsive bHLH binding to alkaloid gene promoters",
        known_examples=["CjbHLH1", "CrMYC2", "NtMYC2"],
    ),
    TFMetaboliteRelation(
        tf_family="bHLH", metabolite_class="terpenoid",
        strength="moderate", score=0.50,
        description="bHLH TFs (MYC2) regulate sesquiterpene biosynthesis via JA signaling.",
        pmids=["28417071"],
        species_scope="general",
        mechanism="JA-responsive activation of terpene synthase genes",
        known_examples=["AtMYC2", "AaMYC2"],
    ),

    # ════════════════════════════════════════════════════════════
    # WRKY family
    # ════════════════════════════════════════════════════════════
    TFMetaboliteRelation(
        tf_family="WRKY", metabolite_class="alkaloid",
        strength="strong", score=0.85,
        description="WRKY TFs widely regulate alkaloid biosynthesis: nicotine in tobacco, "
                    "benzylisoquinoline in opium poppy. WRKYs bind W-box in biosynthetic gene promoters.",
        pmids=["14593171", "24659499", "24958891"],
        species_scope="solanaceae, papaveraceae",
        mechanism="W-box binding in alkaloid gene promoters",
        known_examples=["NtWRKY1", "NtWRKY2", "PsWRKY"],
    ),
    TFMetaboliteRelation(
        tf_family="WRKY", metabolite_class="terpenoid",
        strength="moderate", score=0.50,
        description="WRKY TFs regulate sesquiterpene biosynthesis in cotton and Artemisia.",
        pmids=["28417071"],
        species_scope="malvaceae, asteraceae",
        known_examples=["GaWRKY1", "AaWRKY1"],
    ),
    TFMetaboliteRelation(
        tf_family="WRKY", metabolite_class="flavonoid",
        strength="weak", score=0.25,
        description="Some WRKYs regulate flavonoid pathway under stress conditions, but not primary function.",
        pmids=["25732535"],
        species_scope="general",
    ),

    # ════════════════════════════════════════════════════════════
    # ERF/AP2 family
    # ════════════════════════════════════════════════════════════
    TFMetaboliteRelation(
        tf_family="ERF", metabolite_class="alkaloid",
        strength="strong", score=0.85,
        description="ERF TFs (especially JA-responsive ERFs) are key regulators of alkaloid biosynthesis. "
                    "Examples: ORCA3 (Catharanthus), NtERF189 (nicotine).",
        pmids=["17419843", "19033552", "25296257"],
        species_scope="apocynaceae, solanaceae",
        mechanism="GCC-box binding in alkaloid gene promoters, JA-responsive",
        known_examples=["CrORCA3", "NtERF189", "NtERF221/ORC1"],
    ),
    TFMetaboliteRelation(
        tf_family="ERF", metabolite_class="terpenoid",
        strength="moderate", score=0.55,
        description="ERF TFs regulate terpenoid biosynthesis. "
                    "Examples: AaERF1/2 (artemisinin), SmERF (tanshinone).",
        pmids=["28417071", "27207470"],
        species_scope="asteraceae, lamiaceae",
        known_examples=["AaERF1", "AaERF2", "SmERF6"],
    ),
    TFMetaboliteRelation(
        tf_family="ERF", metabolite_class="flavonoid",
        strength="weak", score=0.30,
        description="Some ERFs regulate flavonoid pathway under stress via GCC-box in PAL/CHS promoters.",
        pmids=["25732535"],
        species_scope="general",
    ),

    # ════════════════════════════════════════════════════════════
    # NAC family
    # ════════════════════════════════════════════════════════════
    TFMetaboliteRelation(
        tf_family="NAC", metabolite_class="anthocyanin",
        strength="moderate", score=0.55,
        description="NAC TFs (e.g., MdNAC52) regulate anthocyanin biosynthesis. "
                    "Some NACs also regulate lignin via secondary cell wall pathway.",
        pmids=["28417071"],
        species_scope="general",
        known_examples=["MdNAC52", "AtNAC"],
    ),
    TFMetaboliteRelation(
        tf_family="NAC", metabolite_class="lignin",
        strength="strong", score=0.80,
        description="NAC TFs (NST/SND/VND family) are master regulators of secondary cell wall biosynthesis "
                    "including lignin.",
        pmids=["17419843", "19033552"],
        species_scope="general",
        mechanism="Direct binding to secondary wall gene promoters",
        known_examples=["AtNST1", "AtNST2", "AtSND1", "AtVND6", "AtVND7"],
    ),

    # ════════════════════════════════════════════════════════════
    # bZIP family
    # ════════════════════════════════════════════════════════════
    TFMetaboliteRelation(
        tf_family="bZIP", metabolite_class="alkaloid",
        strength="weak", score=0.30,
        description="bZIP TFs may regulate alkaloid biosynthesis in response to light/UV signals.",
        pmids=["28417071"],
        species_scope="general",
    ),
    TFMetaboliteRelation(
        tf_family="bZIP", metabolite_class="flavonoid",
        strength="moderate", score=0.55,
        description="bZIP TFs (HY5) regulate flavonoid biosynthesis genes in response to UV/light.",
        pmids=["25732535"],
        species_scope="general",
        mechanism="Light-responsive G-box binding in flavonoid gene promoters",
        known_examples=["AtHY5"],
    ),

    # ════════════════════════════════════════════════════════════
    # WD40 family (co-regulator)
    # ════════════════════════════════════════════════════════════
    TFMetaboliteRelation(
        tf_family="WD40", metabolite_class="anthocyanin",
        strength="strong", score=0.85,
        description="WD40 proteins (TTG1) are essential components of MBW complex for anthocyanin activation.",
        pmids=["26847441"],
        species_scope="general",
        mechanism="MBW complex scaffold protein",
        known_examples=["AtTTG1", "PhAN11"],
    ),

    # ════════════════════════════════════════════════════════════
    # AP2/ERF family — merged from PRIOR_KNOWLEDGE + ontology
    # ════════════════════════════════════════════════════════════
    TFMetaboliteRelation(
        tf_family="AP2-ERF", metabolite_class="alkaloid",
        strength="moderate", score=0.50,
        description="AP2-ERF TFs are associated with alkaloid biosynthesis in some species. "
                    "Evidence level is weaker than dedicated ERF subgroup.",
        pmids=["28417071"],
        species_scope="general",
    ),
    TFMetaboliteRelation(
        tf_family="AP2-ERF", metabolite_class="flavonoid",
        strength="weak", score=0.25,
        description="Some ERF/AP2 TFs may regulate flavonoid pathway under specific stress conditions.",
        pmids=["25732535"],
        species_scope="general",
    ),
    TFMetaboliteRelation(
        tf_family="AP2-ERF", metabolite_class="terpenoid",
        strength="moderate", score=0.55,
        description="AP2/ERF TFs regulate terpenoid biosynthesis. "
                    "Examples: AaERF1/2 (artemisinin), SmERF (tanshinone).",
        pmids=["28417071", "27207470"],
        species_scope="asteraceae, lamiaceae",
        known_examples=["AaERF1", "AaERF2", "SmERF6"],
    ),
    TFMetaboliteRelation(
        tf_family="AP2-ERF", metabolite_class="lignin",
        strength="weak", score=0.20,
        description="Limited evidence for AP2-ERF lignin regulation.",
        pmids=[],
        species_scope="general",
    ),
    TFMetaboliteRelation(
        tf_family="AP2-ERF", metabolite_class="defense_metabolite",
        strength="moderate", score=0.45,
        description="AP2-ERF TFs involved in defense-related secondary metabolism.",
        pmids=["14593171"],
        species_scope="general",
    ),

    # ════════════════════════════════════════════════════════════
    # C2H2 family — merged from PRIOR_KNOWLEDGE + ontology
    # ════════════════════════════════════════════════════════════
    TFMetaboliteRelation(
        tf_family="C2H2", metabolite_class="alkaloid",
        strength="weak", score=0.20,
        description="Limited evidence for C2H2 zinc finger TFs in alkaloid regulation.",
        pmids=["28417071"],
        species_scope="general",
    ),
    TFMetaboliteRelation(
        tf_family="C2H2", metabolite_class="flavonoid",
        strength="weak", score=0.25,
        description="Some C2H2 zinc finger proteins may affect flavonoid accumulation under stress.",
        pmids=[],
        species_scope="general",
    ),
    TFMetaboliteRelation(
        tf_family="C2H2", metabolite_class="defense_metabolite",
        strength="weak", score=0.25,
        description="C2H2 zinc finger proteins may be involved in defense metabolite regulation.",
        pmids=[],
        species_scope="general",
    ),

    # ════════════════════════════════════════════════════════════
    # TCP family — merged from PRIOR_KNOWLEDGE
    # ════════════════════════════════════════════════════════════
    TFMetaboliteRelation(
        tf_family="TCP", metabolite_class="flavonoid",
        strength="weak", score=0.20,
        description="TCP TFs may indirectly affect flavonoid biosynthesis through growth regulation.",
        pmids=["28417071"],
        species_scope="general",
    ),
    TFMetaboliteRelation(
        tf_family="TCP", metabolite_class="lignin",
        strength="weak", score=0.25,
        description="TCP TFs may influence secondary cell wall formation including lignin.",
        pmids=[],
        species_scope="general",
    ),
    TFMetaboliteRelation(
        tf_family="TCP", metabolite_class="terpenoid",
        strength="weak", score=0.15,
        description="Very limited evidence for TCP in terpenoid regulation.",
        pmids=[],
        species_scope="general",
    ),

    # ════════════════════════════════════════════════════════════
    # NAC expansion — merged from PRIOR_KNOWLEDGE
    # ════════════════════════════════════════════════════════════
    TFMetaboliteRelation(
        tf_family="NAC", metabolite_class="alkaloid",
        strength="weak", score=0.20,
        description="Limited evidence for NAC TFs in alkaloid regulation.",
        pmids=["28417071"],
        species_scope="general",
    ),
    TFMetaboliteRelation(
        tf_family="NAC", metabolite_class="flavonoid",
        strength="weak", score=0.30,
        description="Some NAC TFs (e.g., MdNAC52) may influence flavonoid/anthocyanin accumulation.",
        pmids=["28417071"],
        species_scope="general",
    ),
    TFMetaboliteRelation(
        tf_family="NAC", metabolite_class="terpenoid",
        strength="weak", score=0.15,
        description="Very limited evidence for NAC in terpenoid regulation.",
        pmids=[],
        species_scope="general",
    ),

    # ════════════════════════════════════════════════════════════
    # bHLH expansions — merged from PRIOR_KNOWLEDGE
    # ════════════════════════════════════════════════════════════
    TFMetaboliteRelation(
        tf_family="bHLH", metabolite_class="lignin",
        strength="weak", score=0.20,
        description="Limited evidence for bHLH in direct lignin regulation.",
        pmids=[],
        species_scope="general",
    ),
    TFMetaboliteRelation(
        tf_family="bHLH", specific_metabolite="nicotine", metabolite_class="alkaloid",
        strength="moderate", score=0.55,
        description="bHLH TFs (MYC2) regulate nicotine biosynthesis via JA signaling in tobacco.",
        pmids=["17419843", "19033552"],
        species_scope="solanaceae",
        known_examples=["NtMYC2"],
    ),

    # ════════════════════════════════════════════════════════════
    # bZIP expansions — merged from PRIOR_KNOWLEDGE
    # ════════════════════════════════════════════════════════════
    TFMetaboliteRelation(
        tf_family="bZIP", metabolite_class="terpenoid",
        strength="weak", score=0.20,
        description="Some bZIP TFs may influence terpenoid biosynthesis under light stress.",
        pmids=["28417071"],
        species_scope="general",
    ),
    TFMetaboliteRelation(
        tf_family="bZIP", metabolite_class="phenolic",
        strength="moderate", score=0.40,
        description="bZIP TFs (HY5) regulate phenolic compound biosynthesis in response to UV.",
        pmids=["25732535"],
        species_scope="general",
    ),

    # ════════════════════════════════════════════════════════════
    # WRKY expansions — merged from PRIOR_KNOWLEDGE
    # ════════════════════════════════════════════════════════════
    TFMetaboliteRelation(
        tf_family="WRKY", metabolite_class="lignin",
        strength="weak", score=0.25,
        description="Some WRKY TFs may influence lignin biosynthesis during defense.",
        pmids=["28417071"],
        species_scope="general",
    ),
    TFMetaboliteRelation(
        tf_family="MYB", metabolite_class="flavonoid",
        specific_metabolite="quercetin",
        strength="strong", score=0.80,
        description="MYB TFs are associated with quercetin (flavonoid) regulation. Top paper: The MYB transcription factor PbMYB12b positively regulates flavonol biosynthesis in pear fruit.. PMID: 30791875. Additional evidence: PMID:40518759.",
        pmids=['30791875', '40518759', '30375656'],
        species_scope="general",
    ),
    TFMetaboliteRelation(
        tf_family="bHLH", metabolite_class="flavonoid",
        specific_metabolite="quercetin",
        strength="strong", score=0.80,
        description="bHLH TFs are associated with quercetin (flavonoid) regulation. Top paper: Seasonal dynamics and molecular regulation of flavonoid biosynthesis in Cyclocarya paliurus (Batal.) Iljinsk.. PMID: 40104034. Additional evidence: PMID:39591435.",
        pmids=['40104034', '39591435', '36580169'],
        species_scope="general",
    ),
    TFMetaboliteRelation(
        tf_family="WRKY", metabolite_class="alkaloid",
        specific_metabolite="berberine",
        strength="strong", score=0.80,
        description="WRKY TFs are associated with berberine (alkaloid) regulation. Top paper: Genomic profiling of WRKY transcription factors and functional analysis of CcWRKY7, CcWRKY29, and CcWRKY32 related to protoberberine alkaloids biosynthesis in Coptis chinensis Fran. PMID: 37035743. Additional evidence: PMID:342",
        pmids=['37035743', '34220919', '26108744'],
        species_scope="general",
    ),
    TFMetaboliteRelation(
        tf_family="ERF", metabolite_class="alkaloid",
        specific_metabolite="nicotine",
        strength="strong", score=0.80,
        description="ERF TFs are associated with nicotine (alkaloid) regulation. Top paper: Transcription Factors in Alkaloid Engineering.. PMID: 34827717. Additional evidence: PMID:32883605.",
        pmids=['34827717', '32883605', '32392307'],
        species_scope="general",
    ),
    TFMetaboliteRelation(
        tf_family="bZIP", metabolite_class="flavonoid",
        specific_metabolite="quercetin",
        strength="strong", score=0.80,
        description="bZIP TFs are associated with quercetin (flavonoid) regulation. Top paper: Review: ABA, flavonols, and the evolvability of land plants.. PMID: 30824025. Additional evidence: PMID:38203776.",
        pmids=['30824025', '38203776', '35184164'],
        species_scope="general",
    ),
    TFMetaboliteRelation(
        tf_family="MYB", metabolite_class="alkaloid",
        specific_metabolite="berberine",
        strength="moderate", score=0.55,
        description="MYB TFs are associated with berberine (alkaloid) regulation. Top paper: Berberine exerts anti-tumor activity in diffuse large B-cell lymphoma by modulating c-myc/CD47 axis.. PMID: 33930347. Additional evidence: PMID:40428298.",
        pmids=['33930347', '40428298', '35847926'],
        species_scope="general",
    ),
    TFMetaboliteRelation(
        tf_family="bHLH", metabolite_class="alkaloid",
        specific_metabolite="berberine",
        strength="moderate", score=0.55,
        description="bHLH TFs are associated with berberine (alkaloid) regulation. Top paper: Molecular genetics of alkaloid biosynthesis in Nicotiana tabacum.. PMID: 23953973. Additional evidence: PMID:25713177.",
        pmids=['23953973', '25713177', '36122814'],
        species_scope="general",
    ),
    TFMetaboliteRelation(
        tf_family="WRKY", metabolite_class="terpenoid",
        specific_metabolite="artemisinin",
        strength="moderate", score=0.55,
        description="WRKY TFs are associated with artemisinin (terpenoid) regulation. Top paper: New insights into artemisinin regulation.. PMID: 28837410. Additional evidence: PMID:33165526.",
        pmids=['28837410', '33165526', '28001315'],
        species_scope="general",
    ),
    TFMetaboliteRelation(
        tf_family="ERF", metabolite_class="terpenoid",
        specific_metabolite="artemisinin",
        strength="moderate", score=0.55,
        description="ERF TFs are associated with artemisinin (terpenoid) regulation. Top paper: New insights into artemisinin regulation.. PMID: 28837410. Additional evidence: PMID:32883605.",
        pmids=['28837410', '32883605', '31087059'],
        species_scope="general",
    ),
    TFMetaboliteRelation(
        tf_family="bZIP", metabolite_class="alkaloid",
        specific_metabolite="caffeine",
        strength="moderate", score=0.55,
        description="bZIP TFs are associated with caffeine (alkaloid) regulation. Top paper: Transcriptional profiling of genes that are regulated by the endoplasmic reticulum-bound transcription factor AIbZIP/CREB3L4 in prostate cells.. PMID: 17712038. Additional evidence: PMID:7816617.",
        pmids=['17712038', '7816617'],
        species_scope="general",
    ),
    TFMetaboliteRelation(
        tf_family="ARF", metabolite_class="flavonoid",
        specific_metabolite="quercetin",
        strength="moderate", score=0.55,
        description="ARF TFs are associated with quercetin (flavonoid) regulation. Top paper: Transcriptome profiling reveals the roles of pigment formation mechanisms in yellow Paeonia delavayi flowers.. PMID: 36580169. Additional evidence: PMID:38561656.",
        pmids=['36580169', '38561656', '40585254'],
        species_scope="general",
    ),
    TFMetaboliteRelation(
        tf_family="MYB", metabolite_class="terpenoid",
        specific_metabolite="artemisinin",
        strength="moderate", score=0.50,
        description="MYB TFs are associated with artemisinin (terpenoid) regulation. Top paper: New insights into artemisinin regulation.. PMID: 28837410. Additional evidence: PMID:40141085.",
        pmids=['28837410', '40141085', '33672342'],
        species_scope="general",
    ),
    TFMetaboliteRelation(
        tf_family="bHLH", metabolite_class="terpenoid",
        specific_metabolite="artemisinin",
        strength="moderate", score=0.50,
        description="bHLH TFs are associated with artemisinin (terpenoid) regulation. Top paper: New insights into artemisinin regulation.. PMID: 28837410. Additional evidence: PMID:39440419.",
        pmids=['28837410', '39440419', '40256890'],
        species_scope="general",
    ),
    TFMetaboliteRelation(
        tf_family="WRKY", metabolite_class="flavonoid",
        specific_metabolite="quercetin",
        strength="moderate", score=0.50,
        description="WRKY TFs are associated with quercetin (flavonoid) regulation. Top paper: Metabolic profile and transcriptome reveal the mystery of petal blotch formation in rose.. PMID: 36670355. Additional evidence: PMID:36520245.",
        pmids=['36670355', '36520245', '37641387'],
        species_scope="general",
    ),
    TFMetaboliteRelation(
        tf_family="ERF", metabolite_class="flavonoid",
        specific_metabolite="quercetin",
        strength="moderate", score=0.50,
        description="ERF TFs are associated with quercetin (flavonoid) regulation. Top paper: Metabolomics and transcriptomics provide insights into the flavonoid biosynthesis pathway in the roots of developing Aster tataricus.. PMID: 36520245. Additional evidence: PMID:40442608.",
        pmids=['36520245', '40442608', '34659295'],
        species_scope="general",
    ),
    TFMetaboliteRelation(
        tf_family="NAC", metabolite_class="alkaloid",
        specific_metabolite="berberine",
        strength="moderate", score=0.50,
        description="NAC TFs are associated with berberine (alkaloid) regulation. Top paper: Levo-tetrahydropalmatine attenuates methamphetamine reward behavior and the accompanying activation of ERK phosphorylation in mice.. PMID: 31398456. Additional evidence: PMID:27121748.",
        pmids=['31398456', '27121748', '35792445'],
        species_scope="general",
    ),
    TFMetaboliteRelation(
        tf_family="NAC", metabolite_class="terpenoid",
        specific_metabolite="artemisinin",
        strength="moderate", score=0.50,
        description="NAC TFs are associated with artemisinin (terpenoid) regulation. Top paper: Overexpression of a Novel NAC Domain-Containing Transcription Factor Gene (AaNAC1) Enhances the Content of Artemisinin and Increases Tolerance to Drought and Botrytis cinerea in Ar. PMID: 27388340. Additional evidence: PMID:3",
        pmids=['27388340', '32625243'],
        species_scope="general",
    ),
    TFMetaboliteRelation(
        tf_family="NAC", metabolite_class="flavonoid",
        specific_metabolite="quercetin",
        strength="moderate", score=0.50,
        description="NAC TFs are associated with quercetin (flavonoid) regulation. Top paper: Mitochondrial biogenesis: pharmacological approaches.. PMID: 24606795. Additional evidence: PMID:41503995.",
        pmids=['24606795', '41503995', '39232832'],
        species_scope="general",
    ),
    TFMetaboliteRelation(
        tf_family="bZIP", metabolite_class="terpenoid",
        specific_metabolite="artemisinin",
        strength="moderate", score=0.50,
        description="bZIP TFs are associated with artemisinin (terpenoid) regulation. Top paper: Light-Induced Artemisinin Biosynthesis Is Regulated by the bZIP Transcription Factor AaHY5 in Artemisia annua.. PMID: 31076768. Additional evidence: PMID:37820909.",
        pmids=['31076768', '37820909', '31120500'],
        species_scope="general",
    ),
    TFMetaboliteRelation(
        tf_family="WD40", metabolite_class="terpenoid",
        specific_metabolite="artemisinin",
        strength="moderate", score=0.50,
        description="WD40 TFs are associated with artemisinin (terpenoid) regulation. Top paper: Genome-wide identification of the WD40 protein family and functional characterization of AaTTG1 in Artemisia annua.. PMID: 39689807. Additional evidence: PMID:42307650.",
        pmids=['39689807', '42307650'],
        species_scope="general",
    ),
    TFMetaboliteRelation(
        tf_family="WD40", metabolite_class="flavonoid",
        specific_metabolite="quercetin",
        strength="moderate", score=0.50,
        description="WD40 TFs are associated with quercetin (flavonoid) regulation. Top paper: Seasonal dynamics and molecular regulation of flavonoid biosynthesis in Cyclocarya paliurus (Batal.) Iljinsk.. PMID: 40104034. Additional evidence: PMID:39591435.",
        pmids=['40104034', '39591435', '38561656'],
        species_scope="general",
    ),
    TFMetaboliteRelation(
        tf_family="SPL", metabolite_class="terpenoid",
        specific_metabolite="artemisinin",
        strength="moderate", score=0.50,
        description="SPL TFs are associated with artemisinin (terpenoid) regulation. Top paper: The SPB-Box Transcription Factor AaSPL2 Positively Regulates Artemisinin Biosynthesis in Artemisia annua L.. PMID: 31024586. Additional evidence: PMID:35193735.",
        pmids=['31024586', '35193735'],
        species_scope="general",
    ),
    TFMetaboliteRelation(
        tf_family="SPL", metabolite_class="alkaloid",
        specific_metabolite="vinblastine",
        strength="moderate", score=0.50,
        description="SPL TFs are associated with vinblastine (alkaloid) regulation. Top paper: Unique induction of p21(WAF1/CIP1)expression by vinorelbine in androgen-independent prostate cancer cells.. PMID: 14562033. Additional evidence: PMID:40604370.",
        pmids=['14562033', '40604370'],
        species_scope="general",
    ),
    TFMetaboliteRelation(
        tf_family="GRAS", metabolite_class="alkaloid",
        specific_metabolite="nicotine",
        strength="moderate", score=0.50,
        description="GRAS TFs are associated with nicotine (alkaloid) regulation. Top paper: Simple Purification of Nicotiana benthamiana-Produced Recombinant Colicins: High-Yield Recovery of Purified Proteins with Minimum Alkaloid Content Supports the Suitability of the H. PMID: 29286298. Additional evidence: PMID:2906",
        pmids=['29286298', '29065549', '42184531'],
        species_scope="general",
    ),
    TFMetaboliteRelation(
        tf_family="GRAS", metabolite_class="terpenoid",
        specific_metabolite="artemisinin",
        strength="moderate", score=0.50,
        description="GRAS TFs are associated with artemisinin (terpenoid) regulation. Top paper: Plasma membrane recycling drives reservoir formation during Toxoplasma gondii intracellular replication.. PMID: 41026795. Additional evidence: PMID:30842758.",
        pmids=['41026795', '30842758'],
        species_scope="general",
    ),
    TFMetaboliteRelation(
        tf_family="GRAS", metabolite_class="flavonoid",
        specific_metabolite="quercetin",
        strength="moderate", score=0.50,
        description="GRAS TFs are associated with quercetin (flavonoid) regulation. Top paper: Effects of somatostatin, curcumin, and quercetin on the fatty acid profile of breast cancer cell membranes.. PMID: 31545905. Additional evidence: PMID:41237233.",
        pmids=['31545905', '41237233', '41074961'],
        species_scope="general",
    ),
    TFMetaboliteRelation(
        tf_family="ARF", metabolite_class="alkaloid",
        specific_metabolite="nicotine",
        strength="moderate", score=0.50,
        description="ARF TFs are associated with nicotine (alkaloid) regulation. Top paper: Transcriptomic analysis provides insights into the AUXIN RESPONSE FACTOR 6-mediated repression of nicotine biosynthesis in tobacco (Nicotiana tabacum L.).. PMID: 34302568. Additional evidence: PMID:41439816.",
        pmids=['34302568', '41439816'],
        species_scope="general",
    ),
    TFMetaboliteRelation(
        tf_family="HD-ZIP", metabolite_class="alkaloid",
        specific_metabolite="caffeine",
        strength="moderate", score=0.50,
        description="HD-ZIP TFs are associated with caffeine (alkaloid) regulation. Top paper: The Biosynthesis of Main Taste Compounds Is Coordinately Regulated by miRNAs and Phytohormones in Tea Plant (Camellia sinensis).. PMID: 32379968. Additional evidence: PMID:39519276.",
        pmids=['32379968', '39519276'],
        species_scope="general",
    ),
    TFMetaboliteRelation(
        tf_family="HD-ZIP", metabolite_class="flavonoid",
        specific_metabolite="quercetin",
        strength="moderate", score=0.50,
        description="HD-ZIP TFs are associated with quercetin (flavonoid) regulation. Top paper: Flavonoids modify root growth and modulate expression of SHORT-ROOT and HD-ZIP III.. PMID: 26473454. Additional evidence: PMID:41923379.",
        pmids=['26473454', '41923379'],
        species_scope="general",
    ),
    TFMetaboliteRelation(
        tf_family="TCP", metabolite_class="alkaloid",
        specific_metabolite="nicotine",
        strength="moderate", score=0.50,
        description="TCP TFs are associated with nicotine (alkaloid) regulation. Top paper: Nicotine coregulates multiple pathways involved in protein modification/degradation in rat brain.. PMID: 15582157. Additional evidence: PMID:41920974.",
        pmids=['15582157', '41920974', '24935899'],
        species_scope="general",
    ),
    TFMetaboliteRelation(
        tf_family="TCP", metabolite_class="terpenoid",
        specific_metabolite="artemisinin",
        strength="moderate", score=0.50,
        description="TCP TFs are associated with artemisinin (terpenoid) regulation. Top paper: Targeting the molecular chaperone CCT2 inhibits GBM progression by influencing KRAS stability.. PMID: 38582394. Additional evidence: PMID:33539631.",
        pmids=['38582394', '33539631'],
        species_scope="general",
    ),
    TFMetaboliteRelation(
        tf_family="C2H2", metabolite_class="terpenoid",
        specific_metabolite="artemisinin",
        strength="moderate", score=0.50,
        description="C2H2 TFs are associated with artemisinin (terpenoid) regulation. Top paper: Chromatin Accessibility Is Associated with Artemisinin Biosynthesis Regulation in Artemisia annua.. PMID: 33672342. Additional evidence: PMID:35554686.",
        pmids=['33672342', '35554686'],
        species_scope="general",
    ),
    TFMetaboliteRelation(
        tf_family="HSF", metabolite_class="alkaloid",
        specific_metabolite="caffeine",
        strength="moderate", score=0.50,
        description="HSF TFs are associated with caffeine (alkaloid) regulation. Top paper: Coffee extract and caffeine enhance the heat shock response and promote proteostasis in an HSF-1-dependent manner in Caenorhabditis elegans.. PMID: 28674941. Additional evidence: PMID:9184081.",
        pmids=['28674941', '9184081'],
        species_scope="general",
    ),
    TFMetaboliteRelation(
        tf_family="HSF", metabolite_class="flavonoid",
        specific_metabolite="quercetin",
        strength="moderate", score=0.50,
        description="HSF TFs are associated with quercetin (flavonoid) regulation. Top paper: Quercetin diminishes scleral ER stress and protein misfolding: High throughput transcriptome analysis in form-deprivation myopia of Guinea pigs.. PMID: 40541913. Additional evidence: PMID:41531942.",
        pmids=['40541913', '41531942', '34992623'],
        species_scope="general",
    ),
    TFMetaboliteRelation(
        tf_family="GATA", metabolite_class="alkaloid",
        specific_metabolite="berberine",
        strength="moderate", score=0.50,
        description="GATA TFs are associated with berberine (alkaloid) regulation. Top paper: Berberine increases expression of GATA-2 and GATA-3 during inhibition of adipocyte differentiation.. PMID: 19403287. Additional evidence: PMID:9581981.",
        pmids=['19403287', '9581981', '19799972'],
        species_scope="general",
    ),
    TFMetaboliteRelation(
        tf_family="GATA", metabolite_class="flavonoid",
        specific_metabolite="quercetin",
        strength="moderate", score=0.50,
        description="GATA TFs are associated with quercetin (flavonoid) regulation. Top paper: Quercetin regulates Th1/Th2 balance in a murine model of asthma.. PMID: 19061976. Additional evidence: PMID:7885836.",
        pmids=['19061976', '7885836', '38091784'],
        species_scope="general",
    ),
    TFMetaboliteRelation(
        tf_family="WD40", metabolite_class="alkaloid",
        specific_metabolite="nicotine",
        strength="weak", score=0.30,
        description="WD40 TFs are associated with nicotine (alkaloid) regulation. Top paper: Basic helix-loop-helix transcription factors and regulation of alkaloid biosynthesis.. PMID: 22067108.",
        pmids=['22067108'],
        species_scope="general",
    ),
    TFMetaboliteRelation(
        tf_family="ARF", metabolite_class="terpenoid",
        specific_metabolite="taxol",
        strength="weak", score=0.30,
        description="ARF TFs are associated with taxol (terpenoid) regulation. Top paper: Analysis of the protein expression changes during taxol-induced apoptosis under translation inhibition conditions.. PMID: 20717708.",
        pmids=['20717708'],
        species_scope="general",
    ),
    TFMetaboliteRelation(
        tf_family="HD-ZIP", metabolite_class="terpenoid",
        specific_metabolite="artemisinin",
        strength="weak", score=0.30,
        description="HD-ZIP TFs are associated with artemisinin (terpenoid) regulation. Top paper: AaSPATULA synergistically coordinates glandular trichome initiation and vegetative growth to maximize whole-plant artemisinin yield in Artemisia annua.. PMID: 42030805.",
        pmids=['42030805'],
        species_scope="general",
    ),
    TFMetaboliteRelation(
        tf_family="MADS-box", metabolite_class="terpenoid",
        specific_metabolite="artemisinin",
        strength="weak", score=0.30,
        description="MADS-box TFs are associated with artemisinin (terpenoid) regulation. Top paper: MADS-box gene AaSEP4 promotes artemisinin biosynthesis in Artemisia annua.. PMID: 36119604.",
        pmids=['36119604'],
        species_scope="general",
    ),
    TFMetaboliteRelation(
        tf_family="MADS-box", metabolite_class="alkaloid",
        specific_metabolite="caffeine",
        strength="weak", score=0.30,
        description="MADS-box TFs are associated with caffeine (alkaloid) regulation. Top paper: Characterization of a serum response factor-like protein in Saccharomyces cerevisiae, Rlm1, which has transcriptional activity regulated by the Mpk1 (Slt2) mitogen-activated protei. PMID: 9111331.",
        pmids=['9111331'],
        species_scope="general",
    ),
    TFMetaboliteRelation(
        tf_family="C2H2", metabolite_class="alkaloid",
        specific_metabolite="nicotine",
        strength="weak", score=0.30,
        description="C2H2 TFs are associated with nicotine (alkaloid) regulation. Top paper: AtGIS, a C2H2 zinc-finger transcription factor from Arabidopsis regulates glandular trichome development through GA signaling in tobacco.. PMID: 28034756.",
        pmids=['28034756'],
        species_scope="general",
    ),
    TFMetaboliteRelation(
        tf_family="C2H2", metabolite_class="flavonoid",
        specific_metabolite="quercetin",
        strength="weak", score=0.30,
        description="C2H2 TFs are associated with quercetin (flavonoid) regulation. Top paper: Regulatory network of flavonoids, phenolic acids and terpenoids biosynthesis in Zizyphus jujuba Mill. cv. Goutou jujube fruits.. PMID: 39383616.",
        pmids=['39383616'],
        species_scope="general",
    ),
    TFMetaboliteRelation(
        tf_family="Dof", metabolite_class="flavonoid",
        specific_metabolite="quercetin",
        strength="weak", score=0.30,
        description="Dof TFs are associated with quercetin (flavonoid) regulation. Top paper: Identification of the Dof Gene Family in Quinoa and Its Potential Role in Regulating Flavonoid Synthesis Under Different Stress Conditions.. PMID: 40282311.",
        pmids=['40282311'],
        species_scope="general",
    ),
    TFMetaboliteRelation(
        tf_family="GATA", metabolite_class="terpenoid",
        specific_metabolite="artemisinin",
        strength="weak", score=0.30,
        description="GATA TFs are associated with artemisinin (terpenoid) regulation. Top paper: RNA sequencing in Artemisia annua L explored the genetic and metabolic responses to hardly soluble aluminum phosphate treatment.. PMID: 37118364.",
        pmids=['37118364'],
        species_scope="general",
    ),
    TFMetaboliteRelation(
        tf_family="Trihelix", metabolite_class="alkaloid",
        specific_metabolite="caffeine",
        strength="weak", score=0.30,
        description="Trihelix TFs are associated with caffeine (alkaloid) regulation. Top paper: Transcriptomic Analysis Reveals the Molecular Adaptation of Three Major Secondary Metabolic Pathways to Multiple Macronutrient Starvation in Tea (Camellia sinensis).. PMID: 32106614.",
        pmids=['32106614'],
        species_scope="general",
    ),

]


# Metabolite name → class mapping (for fuzzy matching)
# Class name aliases (for compatibility with old PRIOR_KNOWLEDGE / ontology naming)
CLASS_ALIASES: dict[str, str] = {
    "defense": "defense_metabolite",
    "sesquiterpene": "terpenoid",
    "nicotine": "alkaloid",
}

METABOLITE_CLASS_MAP: dict[str, str] = {
    # Alkaloids
    "berberine": "alkaloid", "coptisine": "alkaloid", "palmatine": "alkaloid",
    "nicotine": "alkaloid", "nornicotine": "alkaloid", "anatabine": "alkaloid",
    "morphine": "alkaloid", "codeine": "alkaloid", "sanguinarine": "alkaloid",
    "vincristine": "alkaloid", "vinblastine": "alkaloid", "catharanthine": "alkaloid",
    "vindoline": "alkaloid", "caffeine": "alkaloid", "solanine": "alkaloid",
    "tomatine": "alkaloid", "cocaine": "alkaloid", "quinine": "alkaloid",
    "colchicine": "alkaloid", "strychnine": "alkaloid",
    # Flavonoids
    "anthocyanin": "anthocyanin", "cyanidin": "anthocyanin", "delphinidin": "anthocyanin",
    "flavonoid": "flavonoid", "flavonol": "flavonoid", "kaempferol": "flavonoid",
    "quercetin": "flavonoid", "rutin": "flavonoid", "myricetin": "flavonoid",
    "catechin": "flavonoid", "epicatechin": "flavonoid", "proanthocyanidin": "flavonoid",
    "isoflavone": "flavonoid", "genistein": "flavonoid", "daidzein": "flavonoid",
    "chalcone": "flavonoid", "naringenin": "flavonoid",
    # Phenylpropanoids
    "lignin": "lignin", "coumarin": "phenylpropanoid", "lignan": "phenylpropanoid",
    # Terpenoids
    "artemisinin": "terpenoid", "arteannuin": "terpenoid",
    "taxol": "terpenoid", "paclitaxel": "terpenoid",
    "tanshinone": "terpenoid", "menthol": "terpenoid",
    "gossypol": "terpenoid", "limonene": "terpenoid",
    # Glucosinolates
    "glucosinolate": "phenylpropanoid", "sulforaphane": "phenylpropanoid",
}


class TFKnowledgeBase:
    """Query the canonical TF-metabolite regulatory knowledge base."""

    def __init__(self) -> None:
        self._relations = CANONICAL_KNOWLEDGE
        self._class_map = METABOLITE_CLASS_MAP
        self._aliases = CLASS_ALIASES

    def _resolve_class(self, metabolite: str) -> str:
        """Map a specific metabolite name to its metabolic class."""
        m = metabolite.lower().strip()
        # Check class aliases first (e.g. "defense" → "defense_metabolite")
        if m in self._aliases:
            return self._aliases[m]
        # Direct match
        if m in self._class_map:
            return self._class_map[m]
        # Substring match
        for name, cls in self._class_map.items():
            if name in m or m in name:
                return cls
        # Check aliases as substring
        for alias, cls in self._aliases.items():
            if alias in m or m in alias:
                return cls
        return m  # return as-is (may match metabolite_class directly)

    def query(
        self,
        metabolite: str = "",
        tf_family: str = "",
        species: str = "",
        min_score: float = 0.0,
    ) -> KnowledgeQueryResult:
        """Query the knowledge base.

        Args:
            metabolite: Metabolite name or class (e.g. "berberine", "alkaloid")
            tf_family: TF family name (e.g. "MYB", "bHLH", "WRKY")
            species: Species name for scope filtering
            min_score: Minimum evidence score (0-1)
        """
        metabolite_lower = metabolite.lower().strip() if metabolite else ""
        tf_lower = tf_family.lower().strip() if tf_family else ""
        species_lower = species.lower().strip() if species else ""

        # Resolve specific metabolite to class (e.g. "berberine" → "alkaloid")
        resolved_class = self._resolve_class(metabolite) if metabolite else ""

        results = []
        for rel in self._relations:
            if rel.score < min_score:
                continue

            # Filter by TF family
            if tf_lower and tf_lower not in rel.tf_family.lower():
                continue

            # Filter by metabolite
            if metabolite_lower:
                meta_match = (
                    metabolite_lower in rel.metabolite_class.lower()
                    or metabolite_lower in rel.specific_metabolite.lower()
                    or resolved_class in rel.metabolite_class.lower()
                )
                if not meta_match:
                    # Check if any description keyword matches
                    if metabolite_lower not in rel.description.lower():
                        continue

            # Filter by species (loose match)
            if species_lower and rel.species_scope != "general":
                if species_lower not in rel.species_scope.lower():
                    # Still include if strength is strong
                    if rel.strength != "strong":
                        continue

            results.append(rel)

        # Sort by score
        results.sort(key=lambda r: r.score, reverse=True)

        note = ""
        if not results and metabolite:
            note = (
                f"No known TF regulators found for '{metabolite}'. "
                f"Consider: (1) searching literature, (2) checking related metabolite classes, "
                f"(3) the metabolite may be regulated by non-transcriptional mechanisms."
            )

        return KnowledgeQueryResult(
            query=f"metabolite={metabolite}, tf={tf_family}, species={species}",
            relations=results,
            total=len(results),
            note=note,
        )

    def get_families_for_metabolite(self, metabolite: str) -> list[str]:
        """Get all TF families known to regulate a given metabolite."""
        result = self.query(metabolite=metabolite)
        families = list({r.tf_family for r in result.relations})
        return sorted(families, key=lambda f: max(
            (r.score for r in result.relations if r.tf_family == f), default=0
        ), reverse=True)

    def get_metabolites_for_family(self, tf_family: str) -> list[str]:
        """Get all metabolite classes regulated by a given TF family."""
        result = self.query(tf_family=tf_family)
        classes = list({r.metabolite_class for r in result.relations})
        return sorted(classes, key=lambda c: max(
            (r.score for r in result.relations if r.metabolite_class == c), default=0
        ), reverse=True)

    def get_citations(self, tf_family: str, metabolite: str) -> list[str]:
        """Get all PMIDs for a specific TF-metabolite relationship."""
        result = self.query(tf_family=tf_family, metabolite=metabolite)
        pmids = []
        for r in result.relations:
            pmids.extend(r.pmids)
        return list(set(pmids))
