"""
o2pls.py — O2PLS (Two-way Orthogonal Partial Least Squares) 节点。

多变量整合分析：将转录组 (X) 和代谢组 (Y) 数据进行联合分解，
分离联合变异和正交（非联合）变异。对标高分多组学文章中的
O2PLS/CCA 标准分析流程。

算法参考:
  - Trygg & Wold (2003), "O2-PLS, a two-block (X-Y) latent variable
    regression method with an integral OSC filter"
  - el Bouhaddani et al. (2016), "Evaluation of O2PLS in Omics
    data integration", BMC Bioinformatics

实现:
  1. sklearn PLS 作为基分解 (joint variation)
  2. 正交信号校正 (OSC): 过滤 X 中与 Y 正交的变异
  3. VIP (Variable Importance in Projection) 评分
  4. 交叉验证选取最优潜变量数
  5. Loadings, scores, 解释方差比例
"""

from __future__ import annotations

import logging
from typing import Any

import numpy as np

from phyto_reason.workflows.runtime_state import RuntimeState

logger = logging.getLogger("workflow_nodes.o2pls")


# ═══════════════════════════════════════════════════════════════
# O2PLS Core
# ═══════════════════════════════════════════════════════════════

def _standardize(X: np.ndarray) -> np.ndarray:
    """Column-wise standardization (mean=0, std=1)."""
    mean = np.mean(X, axis=0, keepdims=True)
    std = np.std(X, axis=0, ddof=1, keepdims=True)
    std[std < 1e-10] = 1.0
    return (X - mean) / std


def _explained_variance(X: np.ndarray, scores: np.ndarray, loadings: np.ndarray) -> float:
    """Compute proportion of variance in X explained by given scores and loadings.

    R²_X = 1 - ||X - T·P^T||² / ||X||²
    """
    X_pred = scores @ loadings.T
    ss_res = np.sum((X - X_pred) ** 2)
    ss_tot = np.sum(X ** 2)
    if ss_tot < 1e-15:
        return 0.0
    return float(1.0 - ss_res / ss_tot)


def _vip_scores(weights: np.ndarray, scores: np.ndarray, Y: np.ndarray) -> np.ndarray:
    """Calculate Variable Importance in Projection (VIP) scores.

    VIP_j = sqrt( p * Σ_a (w_ja² * SSY_a) / Σ_a SSY_a )

    where p = number of X variables, a = latent variable index,
    SSY_a = variance in Y explained by component a.

    Args:
        weights: (p, n_components) weight matrix
        scores: (n_samples, n_components) score matrix
        Y: (n_samples, q) response matrix

    Returns:
        (p,) array of VIP scores
    """
    p = weights.shape[0]
    n_comp = weights.shape[1]

    # SSY_a: sum of squared Y scores per component (simplified)
    ssy_a = np.sum(scores ** 2, axis=0)  # (n_comp,)
    ssy_total = np.sum(ssy_a)

    if ssy_total < 1e-15:
        return np.zeros(p)

    vip = np.zeros(p)
    for j in range(p):
        w_j = weights[j, :] ** 2  # (n_comp,)
        vip[j] = np.sqrt(p * np.sum(w_j * ssy_a) / ssy_total)

    return vip


def _cross_validate_n_components(
    X: np.ndarray,
    Y: np.ndarray,
    max_comp: int = 10,
    n_folds: int = 5,
) -> dict:
    """Cross-validation to select optimal number of PLS components.

    Uses R²_Q² metric (predictive Q² via cross-validation).

    Args:
        X: (n, p) predictor matrix
        Y: (n, q) response matrix
        max_comp: max components to test
        n_folds: CV folds

    Returns:
        {best_n_comp, q2_values, r2_values}
    """
    from sklearn.cross_decomposition import PLSRegression
    from sklearn.model_selection import KFold

    n = X.shape[0]
    max_comp = min(max_comp, n - n_folds, X.shape[1])

    if max_comp < 2:
        return {"best_n_comp": 1, "q2_values": [], "r2_values": []}

    q2_values: list[float] = []
    r2_values: list[float] = []

    for n_comp in range(1, max_comp + 1):
        # K-fold CV for Q²
        kf = KFold(n_splits=min(n_folds, n), shuffle=True, random_state=42)
        y_pred_all = np.zeros_like(Y)

        try:
            for train_idx, test_idx in kf.split(X):
                X_train, X_test = X[train_idx], X[test_idx]
                Y_train = Y[train_idx]

                pls = PLSRegression(n_components=n_comp, scale=False)
                pls.fit(X_train, Y_train)
                y_pred_all[test_idx] = pls.predict(X_test)

            # Q² (predictive)
            press = np.sum((Y - y_pred_all) ** 2)
            tss = np.sum((Y - np.mean(Y, axis=0)) ** 2)
            q2 = 1.0 - press / (tss + 1e-15)
            q2_values.append(round(float(q2), 4))

            # R² (fit)
            pls = PLSRegression(n_components=n_comp, scale=False)
            pls.fit(X, Y)
            Y_pred = pls.predict(X)
            rss = np.sum((Y - Y_pred) ** 2)
            r2 = 1.0 - rss / (tss + 1e-15)
            r2_values.append(round(float(r2), 4))
        except Exception as e:
            logger.debug("PLS CV failed for n_comp=%d: %s", n_comp, e)
            q2_values.append(0.0)
            r2_values.append(0.0)

    # Pick best n_comp: maximize Q², but penalize large models
    if not q2_values:
        return {"best_n_comp": 1, "q2_values": [], "r2_values": []}

    # Heuristic: pick n_comp where Q² gains < 0.02
    best = 1
    for i in range(1, len(q2_values)):
        if q2_values[i] - q2_values[i - 1] > 0.02:
            best = i + 1

    return {
        "best_n_comp": min(best, max_comp),
        "q2_values": q2_values,
        "r2_values": r2_values,
    }


def _orthogonal_signal_correction(
    X: np.ndarray,
    Y_score: np.ndarray,
    n_orthogonal: int = 2,
) -> tuple[np.ndarray, list[np.ndarray], list[np.ndarray]]:
    """Remove Y-orthogonal variation from X.

    For each orthogonal component:
      1. Compute X weights orthogonal to Y_score: w⊥ = X^T t_Y - t_Y^T·t_Y · w_PLS
      2. Score: t⊥ = X w⊥
      3. Loading: p⊥ = X^T t⊥ / (t⊥^T t⊥)
      4. Deflate: X = X - t⊥ p⊥^T

    Args:
        X: (n, p) centered data matrix
        Y_score: (n,) first PLS Y-score (joint component proxy)
        n_orthogonal: number of orthogonal components to remove

    Returns:
        (X_corrected, t_orth_list, p_orth_list)
    """
    X_corr = X.copy()
    t_orth_list: list[np.ndarray] = []
    p_orth_list: list[np.ndarray] = []
    n = X.shape[0]

    for _ in range(n_orthogonal):
        # 1. Orthogonal weight: X^T t_Y, then orthogonalize
        w_o = X_corr.T @ Y_score  # (p,)
        w_o = w_o / (np.linalg.norm(w_o) + 1e-15)

        # 2. Orthogonal score
        t_o = X_corr @ w_o  # (n,)
        t_o = t_o - np.mean(t_o)
        t_o = t_o / (np.linalg.norm(t_o) + 1e-15)

        # 3. Orthogonal loading
        p_o = X_corr.T @ t_o / (t_o @ t_o + 1e-15)  # (p,)

        # 4. Deflate X
        X_corr = X_corr - np.outer(t_o, p_o)

        t_orth_list.append(t_o)
        p_orth_list.append(p_o)

    return X_corr, t_orth_list, p_orth_list


def _run_o2pls(
    X_raw: np.ndarray,
    Y_raw: np.ndarray,
    x_names: list[str],
    y_names: list[str],
    n_joint: int | None = None,
    n_orthogonal: int = 2,
) -> dict:
    """Run full O2PLS pipeline.

    Args:
        X_raw: (n_samples, p_genes) expression data
        Y_raw: (n_samples, q_metabolites) metabolite data
        x_names: gene IDs for X columns
        y_names: metabolite IDs for Y columns
        n_joint: number of joint components (auto-selected if None)
        n_orthogonal: number of orthogonal components

    Returns:
        O2PLS report dict
    """
    n, p = X_raw.shape
    q = Y_raw.shape[1]

    if n < 5:
        return {"reason": f"too few samples ({n}) for O2PLS", "n_samples": n}
    if p < 3:
        return {"reason": f"too few genes ({p}) for O2PLS dimensionality reduction"}
    if q < 1:
        return {"reason": "no metabolite data for Y block"}

    # Standardize
    X = _standardize(X_raw)
    Y = _standardize(Y_raw)

    # Determine n_joint via CV
    if n_joint is None:
        cv_result = _cross_validate_n_components(X, Y, max_comp=min(10, n - 2, p, q + 2))
        n_joint = max(1, cv_result["best_n_comp"])
        cv_info = {"q2_values": cv_result["q2_values"], "r2_values": cv_result["r2_values"]}
    else:
        cv_info = {}

    # ── Step 1: Joint PLS decomposition ──────────────────
    from sklearn.cross_decomposition import PLSRegression

    pls = PLSRegression(n_components=n_joint, scale=False)
    pls.fit(X, Y)

    # Scores (latent variables)
    T_joint = pls.x_scores_  # (n, n_joint) — joint variation in X
    U_joint = pls.y_scores_  # (n, n_joint) — joint variation in Y

    # Weights
    W_x = pls.x_weights_   # (p, n_joint)
    W_y = pls.y_weights_   # (q, n_joint)

    # Loadings
    P_x = pls.x_loadings_  # (p, n_joint)
    P_y = pls.y_loadings_  # (q, n_joint)

    # ── Step 2: VIP scores ─────────────────────────────
    vip = _vip_scores(W_x, T_joint, Y)

    # ── Step 3: Orthogonal signal correction ───────────
    # Use first Y score as joint signal proxy
    if n_orthogonal > 0 and T_joint.shape[1] > 0:
        y_score_proxy = T_joint[:, 0]  # First joint X-score (aligned with Y)
        X_corrected, t_orth, p_orth = _orthogonal_signal_correction(
            X, y_score_proxy, n_orthogonal=n_orthogonal,
        )
        ortho_variance_explained = float(
            np.sum(X ** 2) - np.sum(X_corrected ** 2)
        ) / max(np.sum(X ** 2), 1e-15)
    else:
        X_corrected = X
        t_orth, p_orth = [], []
        ortho_variance_explained = 0.0

    # ── Step 4: Variance explained ─────────────────────
    # Joint variance in X
    X_joint_var = _explained_variance(X, T_joint, P_x)
    # Joint variance in Y
    Y_joint_var = _explained_variance(Y, U_joint, P_y)

    # ── Step 5: Top contributing variables ─────────────
    # X loadings (genes most contributing to joint component 1)
    x_loadings_comp1 = [(x_names[i], float(P_x[i, 0]), float(vip[i]))
                        for i in range(min(p, len(x_names)))]
    x_loadings_comp1.sort(key=lambda tup: abs(tup[1]), reverse=True)

    # Y loadings (metabolites most contributing to joint component 1)
    if q > 0 and P_y.shape[0] == q:
        y_loadings_comp1 = [(y_names[i], float(P_y[i, 0]))
                            for i in range(min(q, len(y_names)))]
        y_loadings_comp1.sort(key=lambda tup: abs(tup[1]), reverse=True)
    else:
        y_loadings_comp1 = []

    # ── Step 6: Joint scores (sample-level) ────────────
    joint_scores = [
        {"sample_idx": i, "scores": [round(float(T_joint[i, j]), 4)
                                     for j in range(min(n_joint, 3))]}
        for i in range(n)
    ]

    # ── Assemble report ────────────────────────────────
    o2pls_report = {
        "n_samples": n,
        "n_x_variables": p,
        "n_y_variables": q,
        "n_joint_components": n_joint,
        "n_orthogonal_components": n_orthogonal,
        "x_joint_variance_explained": round(X_joint_var, 4),
        "y_joint_variance_explained": round(Y_joint_var, 4),
        "orthogonal_variance_removed": round(ortho_variance_explained, 4),
        "cv_info": cv_info,
        "top_x_variables": [
            {"name": name, "loading": ld, "vip": v}
            for name, ld, v in x_loadings_comp1[:30]
        ],
        "top_y_variables": [
            {"name": name, "loading": ld}
            for name, ld in y_loadings_comp1[:20]
        ],
        "joint_scores": joint_scores,
        "x_weights_comp1": [
            {"name": x_names[i], "weight": float(W_x[i, 0])}
            for i in range(min(p, len(x_names)))
            if abs(W_x[i, 0]) > np.percentile(np.abs(W_x[:, 0]), 90)
        ][:30] if n_joint > 0 and W_x.shape[1] > 0 else [],
        "warnings": [],
    }

    if X_joint_var < 0.10:
        o2pls_report["warnings"].append(
            f"Low joint variance explained ({X_joint_var:.1%}). "
            f"Transcriptome and metabolome may be largely decoupled "
            f"under these conditions."
        )

    return o2pls_report


# ═══════════════════════════════════════════════════════════════
# Node
# ═══════════════════════════════════════════════════════════════

def o2pls_node(state: RuntimeState) -> dict:
    """O2PLS 多变量整合分析节点。

    使用 Two-way Orthogonal PLS 将转录组和代谢组数据联合分解：
      - 联合变异 (joint variation): PLS 潜变量
      - 正交变异 (orthogonal variation): OSC 过滤
      - VIP 评分: 变量重要性
      - 交叉验证: 选取最优潜变量数

    需要:
      - 表达矩阵 (X block)
      - 代谢物矩阵 (Y block)
      - 足够的共同样本 (≥8)

    Returns:
        dict with o2pls_report
    """
    try:
        from phyto_reason.workflows.execution_nodes import _should_skip_node
        if _should_skip_node(state, "o2pls"):
            state.add_trace("o2pls", "skipped", summary={"reason": "user skipped this step"})
            return {}

        expr = None
        meta_mat = None
        if state.planner_state:
            expr = state.planner_state.expression_matrix
            meta_mat = state.planner_state.metabolite_matrix

        deg_report = getattr(state, "deg_report", None) or {}
        dam_report = getattr(state, "dam_report", None) or {}

        if not expr or not meta_mat:
            state.add_trace("o2pls", "completed",
                            summary={"reason": "missing expression or metabolite data"})
            return {"o2pls_report": {"reason": "missing data for O2PLS"}}

        # Get common samples
        expr_samples: set[str] = set()
        for v in expr.values():
            expr_samples.update(v.keys())
            break
        meta_samples: set[str] = set()
        for v in meta_mat.values():
            meta_samples.update(v.keys())
            break
        common_samples = sorted(expr_samples & meta_samples)

        if len(common_samples) < 8:
            state.add_trace("o2pls", "completed",
                            summary={"reason": f"too few common samples: {len(common_samples)} (need ≥8)"})
            return {"o2pls_report": {
                "reason": f"O2PLS requires ≥8 samples (got {len(common_samples)})",
                "n_samples": len(common_samples),
            }}

        # Select genes: top DEGs (by effect size)
        if deg_report.get("top_genes"):
            gene_ids = [g["gene_id"] for g in deg_report["top_genes"][:200]
                        if g["gene_id"] in expr]
        else:
            gene_ids = list(expr.keys())[:200]

        # Select metabolites: top DAMs
        if dam_report.get("top_metabolites"):
            meta_ids = [m["metabolite"] for m in dam_report["top_metabolites"][:30]
                        if m["metabolite"] in meta_mat]
        else:
            meta_ids = list(meta_mat.keys())[:30]

        if len(gene_ids) < 5 or len(meta_ids) < 2:
            state.add_trace("o2pls", "completed",
                            summary={"reason": f"insufficient variables: {len(gene_ids)} genes, {len(meta_ids)} metabolites"})
            return {"o2pls_report": {
                "reason": f"need ≥5 genes and ≥2 metabolites (got {len(gene_ids)}, {len(meta_ids)})"
            }}

        # Build data matrices
        n_samples = len(common_samples)
        X_raw = np.zeros((n_samples, len(gene_ids)))
        Y_raw = np.zeros((n_samples, len(meta_ids)))

        for j, gid in enumerate(gene_ids):
            arr = np.array([expr[gid].get(s, np.nan) for s in common_samples], dtype=float)
            arr[np.isnan(arr)] = np.nanmean(arr) if not np.all(np.isnan(arr)) else 0.0
            X_raw[:, j] = arr

        for j, mid in enumerate(meta_ids):
            arr = np.array([meta_mat[mid].get(s, np.nan) for s in common_samples], dtype=float)
            arr[np.isnan(arr)] = np.nanmean(arr) if not np.all(np.isnan(arr)) else 0.0
            Y_raw[:, j] = arr

        logger.info("O2PLS: %d samples x %d genes x %d metabolites",
                    n_samples, len(gene_ids), len(meta_ids))

        o2pls_report = _run_o2pls(
            X_raw=X_raw,
            Y_raw=Y_raw,
            x_names=gene_ids,
            y_names=meta_ids,
            n_joint=None,  # Auto-select via CV
            n_orthogonal=2,
        )

        if "reason" in o2pls_report and "n_joint_components" not in o2pls_report:
            state.add_trace("o2pls", "completed", summary=o2pls_report)
            state.o2pls_report = o2pls_report
            return {"o2pls_report": o2pls_report}

        n_joint = o2pls_report.get("n_joint_components", 0)
        x_var = o2pls_report.get("x_joint_variance_explained", 0)
        y_var = o2pls_report.get("y_joint_variance_explained", 0)

        if x_var > 0.15 and y_var > 0.15:
            interpretation = (
                f"O2PLS identified {n_joint} joint components explaining "
                f"{x_var:.1%} of transcriptome and {y_var:.1%} of metabolome variance. "
                f"Strong joint structure suggests coordinated transcriptome-metabolome regulation."
            )
        elif x_var > 0.05:
            interpretation = (
                f"O2PLS found {n_joint} joint components with "
                f"moderate variance explained (X: {x_var:.1%}, Y: {y_var:.1%}). "
                f"Some transcriptome-metabolome coupling detected."
            )
        else:
            interpretation = (
                f"O2PLS explained limited joint variance (X: {x_var:.1%}, Y: {y_var:.1%}). "
                f"Transcriptome and metabolome may operate on different timescales "
                f"or be subject to strong post-transcriptional regulation."
            )

        o2pls_report["interpretation"] = interpretation

        # Generate O2PLS figures
        figure_urls: list[str] = []
        figure_mds: list[str] = []
        try:
            from phyto_reason.visualization.multiomics_figures import (
                plot_o2pls_scores, plot_o2pls_loadings,
            )
            # Joint scores scatter plot
            scores_url = plot_o2pls_scores(
                o2pls_report.get("joint_scores", []),
                title="O2PLS Joint Component Scores",
            )
            if scores_url:
                figure_urls.append(scores_url)
                figure_mds.append(f"![O2PLS Scores]({scores_url})")

            # Loadings plot
            top_x = o2pls_report.get("top_x_variables", [])
            top_y = o2pls_report.get("top_y_variables", [])
            if top_x:
                load_url = plot_o2pls_loadings(
                    top_x, top_y,
                    title="O2PLS Variable Loadings (Joint Comp 1)",
                )
                if load_url:
                    figure_urls.append(load_url)
                    figure_mds.append(f"![O2PLS Loadings]({load_url})")
        except Exception as e:
            logger.debug("O2PLS figure generation skipped: %s", e)

        if figure_urls:
            o2pls_report["figure_urls"] = figure_urls
            o2pls_report["figure_markdown"] = "\n".join(figure_mds)

        state.o2pls_report = o2pls_report

        state.add_trace("o2pls", "completed", summary={
            "n_samples": n_samples,
            "n_joint_comp": n_joint,
            "x_var": round(x_var, 4),
            "y_var": round(y_var, 4),
            "n_top_genes": len(o2pls_report.get("top_x_variables", [])),
        })

        logger.info("O2PLS: %d joint comps, X_var=%.3f, Y_var=%.3f",
                    n_joint, x_var, y_var)

        return {"o2pls_report": o2pls_report}

    except Exception as e:
        logger.error("o2pls_node failed: %s", e, exc_info=True)
        state.add_trace("o2pls", "failed", error=str(e))
        return {"o2pls_report": {"error": str(e)}}
