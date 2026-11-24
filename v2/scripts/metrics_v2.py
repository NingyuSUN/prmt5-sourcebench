"""Deterministic, tie-aware ranking, calibration, concordance, and paired scaffold-bootstrap metrics for PRMT5 v2."""
from __future__ import annotations

import math
from typing import Any

import numpy as np
from scipy.stats import spearmanr
from sklearn.metrics import average_precision_score, roc_auc_score


def choose_validation_threshold(y_val: np.ndarray, score_val: np.ndarray) -> float:
    """Select threshold maximizing validation balanced accuracy using exact integer tie-free comparison."""
    y = np.asarray(y_val, dtype=int)
    s = np.asarray(score_val, dtype=float)
    if len(np.unique(y)) != 2:
        raise ValueError("Both classes required to choose validation threshold")
    thresholds = np.r_[np.unique(s), np.nextafter(s.max(), np.inf)]
    pos = int(y.sum())
    neg = int((1 - y).sum())
    best = None
    for thr in thresholds:
        pred = s >= thr
        tp = int(np.sum(pred & (y == 1)))
        tn = int(np.sum((~pred) & (y == 0)))
        cand = (tp * neg + tn * pos, float(thr))
        if best is None or cand > best:
            best = cand
    return float(best[1])


def binary_threshold_metrics(y: np.ndarray, score: np.ndarray, threshold: float) -> dict[str, Any]:
    y = np.asarray(y, dtype=int)
    s = np.asarray(score, dtype=float)
    pred = s >= float(threshold)
    tp = int(np.sum(pred & (y == 1)))
    fp = int(np.sum(pred & (y == 0)))
    fn = int(np.sum((~pred) & (y == 1)))
    tn = int(np.sum((~pred) & (y == 0)))
    sens = tp / (tp + fn) if (tp + fn) > 0 else None
    spec = tn / (tn + fp) if (tn + fp) > 0 else None
    bal_acc = 0.5 * (sens + spec) if (sens is not None and spec is not None) else None
    return {
        "threshold": float(threshold),
        "n": int(len(y)),
        "tp": tp,
        "fp": fp,
        "fn": fn,
        "tn": tn,
        "sensitivity": sens,
        "specificity": spec,
        "balanced_accuracy": bal_acc,
    }


def top_k_screening_metrics(y: np.ndarray, score: np.ndarray, ids: list[str], k: int) -> dict[str, Any]:
    """Tie-aware top-K hit enrichment (expected hits over uniform tie resolution + deterministic count)."""
    y = np.asarray(y, dtype=int)
    s = np.asarray(score, dtype=float)
    n = len(y)
    k = min(n, max(1, int(k)))
    pos = int(y.sum())
    boundary = float(np.sort(s)[-k])
    above = s > boundary
    tied = s == boundary
    above_n = int(above.sum())
    above_hits = int(y[above].sum())
    tie_n = int(tied.sum())
    tie_hits = int(y[tied].sum())
    slots = k - above_n
    exp_hits = above_hits + (slots * tie_hits / tie_n if tie_n > 0 else 0.0)
    rand_hits = k * pos / n if n > 0 else 0.0
    ef = (exp_hits / rand_hits) if rand_hits > 0 else None
    order = sorted(range(n), key=lambda i: (-s[i], str(ids[i])))[:k]
    det_hits = int(y[order].sum())
    return {
        "k": k,
        "n": n,
        "positives_total": pos,
        "expected_hits": float(exp_hits),
        "deterministic_hits": det_hits,
        "random_expected_hits": float(rand_hits),
        "enrichment_factor": float(ef) if ef is not None else None,
        "boundary_tie_size": tie_n,
    }


def ranking_summary(y: np.ndarray, score: np.ndarray, ids: list[str], top_k: int = 20) -> dict[str, Any]:
    y = np.asarray(y, dtype=int)
    s = np.asarray(score, dtype=float)
    pos = int(y.sum())
    n = len(y)
    ap = float(average_precision_score(y, s)) if 0 < pos <= n else 0.0
    auc = float(roc_auc_score(y, s)) if 0 < pos < n else None
    return {
        "n": n,
        "positives": pos,
        "negatives": n - pos,
        "prevalence": pos / n if n else 0.0,
        "average_precision": ap,
        "roc_auc": auc,
        "screening": top_k_screening_metrics(y, s, ids, top_k),
    }


def censor_aware_concordance_index(
    p_target: np.ndarray,
    relations: list[str] | np.ndarray,
    score: np.ndarray,
    min_diff: float = 0.05,
) -> tuple[float | None, int]:
    """Concordance index (C-index) respecting censored activity bounds.

    `p_target` is higher-is-more-potent (e.g. pIC50/pKi or delta_pActivity).
    For relation:
      '='  : exact p_target
      '<'  : IC50 < val => true p_target >= p_target (right-censored upward in potency)
      '>'  : IC50 > val => true p_target <= p_target (left-censored downward in potency)
    Pair (i, j) with p_i > p_j + min_diff is comparable iff:
      i has rel in ('=', '<', '<=') AND j has rel in ('=', '>', '>=').
    """
    p = np.asarray(p_target, dtype=float)
    s = np.asarray(score, dtype=float)
    rels = [str(r) for r in relations]
    n = len(p)
    concordant = 0.0
    total_pairs = 0
    can_be_higher = np.array([r in ("=", "<", "<=", "exact", "interval") for r in rels])
    can_be_lower = np.array([r in ("=", ">", ">=", "exact", "interval") for r in rels])

    for i in range(n):
        if not can_be_higher[i]:
            continue
        for j in range(n):
            if i == j or not can_be_lower[j]:
                continue
            if p[i] - p[j] >= min_diff:
                total_pairs += 1
                if s[i] > s[j]:
                    concordant += 1.0
                elif s[i] == s[j]:
                    concordant += 0.5
    if total_pairs == 0:
        return None, 0
    return float(concordant / total_pairs), int(total_pairs)


def within_group_weighted_c_index(
    groups: list[str] | np.ndarray,
    p_target: np.ndarray,
    relations: list[str] | np.ndarray,
    score: np.ndarray,
) -> dict[str, Any]:
    """Compute pair-weighted and compound-weighted C-index strictly within homogeneous assay/study groups."""
    g_arr = np.asarray([str(g) for g in groups])
    p_arr = np.asarray(p_target, dtype=float)
    r_arr = np.asarray([str(r) for r in relations])
    s_arr = np.asarray(score, dtype=float)

    per_group = {}
    weighted_num = 0.0
    weighted_den = 0
    for g in sorted(set(g_arr)):
        idx = np.flatnonzero(g_arr == g)
        ci, npairs = censor_aware_concordance_index(p_arr[idx], r_arr[idx], s_arr[idx])
        # Exact-only Spearman rho if >= 5 exact points with variance
        exact_mask = (r_arr[idx] == "=")
        rho = None
        if exact_mask.sum() >= 5 and len(np.unique(p_arr[idx][exact_mask])) >= 3 and len(np.unique(s_arr[idx][exact_mask])) >= 2:
            rho_val = spearmanr(p_arr[idx][exact_mask], s_arr[idx][exact_mask]).statistic
            rho = float(rho_val) if np.isfinite(rho_val) else None
        per_group[g] = {
            "n": int(len(idx)),
            "comparable_pairs": int(npairs),
            "c_index": ci,
            "spearman_rho_exact": rho,
        }
        if ci is not None and npairs > 0:
            weighted_num += ci * npairs
            weighted_den += npairs

    return {
        "weighted_c_index": float(weighted_num / weighted_den) if weighted_den > 0 else None,
        "total_comparable_pairs": int(weighted_den),
        "per_group": per_group,
    }


def paired_scaffold_bootstrap(
    scaffolds: list[str] | np.ndarray,
    scores_by_model: dict[str, np.ndarray],
    metric_fn,
    comparisons: list[tuple[str, str]],
    n_draws: int = 2000,
    seed: int = 20261120,
) -> dict[str, Any]:
    """Paired cluster bootstrap resampling Murcko scaffolds with replacement."""
    sc_arr = np.asarray([str(x) for x in scaffolds])
    uniq_sc = np.unique(sc_arr)
    members = [np.flatnonzero(sc_arr == sc) for sc in uniq_sc]
    rng = np.random.default_rng(seed)

    point_est = {m: metric_fn(np.arange(len(sc_arr)), s) for m, s in scores_by_model.items()}
    samples = {m: [] for m in scores_by_model}
    valid_draws = 0
    invalid_draws = 0

    for _ in range(n_draws):
        chosen = rng.integers(0, len(uniq_sc), size=len(uniq_sc))
        idx = np.concatenate([members[c] for c in chosen])
        vals = {}
        ok = True
        for m, s in scores_by_model.items():
            v = metric_fn(idx, s[idx])
            if v is None or not math.isfinite(v):
                ok = False
                break
            vals[m] = float(v)
        if not ok:
            invalid_draws += 1
            continue
        valid_draws += 1
        for m, v in vals.items():
            samples[m].append(v)

    model_ci = {}
    for m, arr_list in samples.items():
        arr = np.asarray(arr_list, dtype=float)
        model_ci[m] = {
            "point": float(point_est[m]) if point_est[m] is not None else None,
            "ci95": [float(np.quantile(arr, 0.025)), float(np.quantile(arr, 0.975))] if len(arr) else [None, None],
        }

    diff_ci = {}
    for m_a, m_b in comparisons:
        arr_a = np.asarray(samples[m_a], dtype=float)
        arr_b = np.asarray(samples[m_b], dtype=float)
        diff = arr_a - arr_b
        pt_diff = float(point_est[m_a] - point_est[m_b])
        p_two = float(2.0 * min(np.mean(diff <= 0.0), np.mean(diff >= 0.0))) if len(diff) else None
        diff_ci[f"{m_a}_minus_{m_b}"] = {
            "point_diff": pt_diff,
            "ci95": [float(np.quantile(diff, 0.025)), float(np.quantile(diff, 0.975))] if len(diff) else [None, None],
            "bootstrap_two_sided_p": p_two,
        }

    return {
        "scaffold_clusters": int(len(uniq_sc)),
        "requested_draws": int(n_draws),
        "valid_draws": int(valid_draws),
        "invalid_draws": int(invalid_draws),
        "seed": int(seed),
        "models": model_ci,
        "paired_differences": diff_ci,
    }
