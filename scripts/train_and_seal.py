#!/usr/bin/env python3
"""Stage 2B: Fit Track A & Track B models on Train, select solely on Validation, and seal Test predictions.

Zero-leakage & hardware rules enforced:
1. Verifies `receipts/protocol_freeze.json` and `manifests/feature_manifest_v2.json` SHA-256 hashes.
2. NEVER opens `data/processed/sealed_test_labels/` (`test_labels_accessed: false`).
3. Caps CPU parallelism strictly to 8 threads (`OMP_NUM_THREADS=8`, `n_jobs=8`).
4. Writes `manifests/model_selection_lock.json`, `data/predictions/sealed_predictions.csv`,
   `data/predictions/validation_predictions.csv`, and `receipts/prediction_seal.json` BEFORE any test scoring.
"""
from __future__ import annotations

import argparse
import collections
import csv
import datetime
import hashlib
import json
import math
import os
from pathlib import Path
import sys
from typing import Any

os.environ["OMP_NUM_THREADS"] = "8"
os.environ["MKL_NUM_THREADS"] = "8"
os.environ["OPENBLAS_NUM_THREADS"] = "8"
os.environ["NUMEXPR_NUM_THREADS"] = "8"

import numpy as np
from scipy.special import expit, logit
from sklearn.ensemble import HistGradientBoostingClassifier, HistGradientBoostingRegressor
from sklearn.linear_model import LogisticRegression, Ridge
from sklearn.metrics import average_precision_score
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
import xgboost as xgb

sys.path.insert(0, str(Path(__file__).resolve().parent))
from metrics_v2 import choose_validation_threshold, within_group_weighted_c_index

SEED = 20261120


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def safe_prob(p: np.ndarray) -> np.ndarray:
    return np.clip(np.asarray(p, dtype=np.float64), 1e-6, 1.0 - 1e-6)


def tanimoto_max_diff(query_fp: np.ndarray, ref_fp: np.ndarray, ref_y: np.ndarray, exclude_same: bool = False) -> np.ndarray:
    """Compute (max Tanimoto to positive train - max Tanimoto to negative train)."""
    q = query_fp.astype(np.float32)
    r = ref_fp.astype(np.float32)
    inter = q @ r.T
    union = q.sum(axis=1)[:, None] + r.sum(axis=1)[None, :] - inter
    sim = np.divide(inter, union, out=np.zeros_like(inter), where=union > 0)
    if exclude_same and q.shape[0] == r.shape[0]:
        np.fill_diagonal(sim, -1.0)
    pos_mask = (ref_y == 1)
    neg_mask = (ref_y == 0)
    max_pos = sim[:, pos_mask].max(axis=1) if pos_mask.any() else np.zeros(len(q))
    max_neg = sim[:, neg_mask].max(axis=1) if neg_mask.any() else np.zeros(len(q))
    return (max_pos - max_neg).astype(np.float64)


def verify_inputs(root: Path) -> tuple[dict, dict]:
    freeze_path = root / "receipts/protocol_freeze.json"
    feat_manifest_path = root / "manifests/feature_manifest_v2.json"
    freeze = json.loads(freeze_path.read_text(encoding="utf-8"))
    feat_meta = json.loads(feat_manifest_path.read_text(encoding="utf-8"))

    for rel, expected in freeze["bound_sha256"].items():
        if sha256_file(root / rel) != expected:
            raise RuntimeError(f"Protocol freeze SHA-256 mismatch on {rel}")
    if sha256_file(root / feat_meta["npz_path"]) != feat_meta["npz_sha256"]:
        raise RuntimeError("Feature NPZ SHA-256 mismatch")
    return freeze, feat_meta


def compute_zero_shot_p0_scores(phys3d_all: np.ndarray, phys_cols: list[str], train_mask: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Pure zero-shot 3D cofactor-state pocket physics score (NO y labels accessed!)."""
    col_idx = {c: i for i, c in enumerate(phys_cols)}
    scaler = StandardScaler().fit(phys3d_all[train_mask])
    z = np.clip(scaler.transform(phys3d_all), -4.0, 4.0)

    # First-principles physical chemistry terms (zero training labels used):
    # 1. Strong AutoDock Vina affinity in 7S1S (PRMT5•MTA) & negative (E_MTA - E_SAM) differential
    # 2. High 3D CrippenO3A pharmacophore & 3D shape Tanimoto complementarity to PRMT5•MTA vs PRMT5•SAM
    # 3. Direct invasion of SAM's sulfonium methyl locus (-sam_ce_min_dist_mta_pose)
    # 4. Low internal conformational strain in the PRMT5•MTA bound pose
    raw_s = (
        -1.0 * z[:, col_idx["vina_e_mta_7s1s"]]
        - 0.6 * z[:, col_idx["delta_vina_mta_minus_sam"]]
        + 1.0 * z[:, col_idx["o3a_norm_s2_mta"]]
        + 1.0 * z[:, col_idx["delta_o3a_mta_minus_sam"]]
        + 0.8 * z[:, col_idx["delta_shape_mta_minus_sam"]]
        - 0.5 * z[:, col_idx["sam_ce_min_dist_mta_pose"]]
        - 0.5 * z[:, col_idx["conformer_strain_s2_kcal"]]
    ) / 3.5
    prob_p0 = safe_prob(expit(raw_s))
    return raw_s, prob_p0


def run_track_a(
    root: Path,
    feat_lookup: dict[str, int],
    fp_ach: np.ndarray,
    fp_chi: np.ndarray,
    physchem: np.ndarray,
    phys3d: np.ndarray,
    phys_cols: list[str],
) -> tuple[dict[str, Any], list[dict[str, Any]], list[dict[str, Any]]]:
    struct_rows = list(csv.DictReader((root / "data/processed/track_a_structures.csv").open()))
    dev_rows = list(csv.DictReader((root / "data/processed/track_a_dev_labels.csv").open()))
    if not set(r["split"] for r in dev_rows) <= {"train", "validation"}:
        raise RuntimeError("Track A dev_labels contained non-development rows!")

    dev_map = {r["row_id"]: r for r in dev_rows}
    idx_all = np.array([feat_lookup[r["compound_id"]] for r in struct_rows], dtype=int)
    splits = np.array([r["split"] for r in struct_rows])
    tr_mask = (splits == "train")
    va_mask = (splits == "validation")
    te_mask = (splits == "test")

    y_tr = np.array([int(dev_map[r["row_id"]]["y_coop_3x"]) for r in struct_rows if r["split"] == "train"], dtype=int)
    y_va = np.array([int(dev_map[r["row_id"]]["y_coop_3x"]) for r in struct_rows if r["split"] == "validation"], dtype=int)
    scaff_tr = np.array([r["scaffold_group"] for r in struct_rows if r["split"] == "train"])

    X_ach = fp_ach[idx_all].astype(np.float32)
    X_chi = fp_chi[idx_all].astype(np.float32)
    X_pc = physchem[idx_all].astype(np.float32)
    X_3d = phys3d[idx_all].astype(np.float32)

    selection_log: dict[str, Any] = {}
    selected_scores_all: dict[str, np.ndarray] = {}
    selected_thresholds: dict[str, float] = {}

    # 1. B0: Constant training prevalence
    prev_tr = float(y_tr.mean())
    s_b0 = np.full(len(struct_rows), prev_tr, dtype=np.float64)
    thr_b0 = choose_validation_threshold(y_va, s_b0[va_mask])
    selected_scores_all["B0"] = s_b0
    selected_thresholds["B0"] = thr_b0
    selection_log["B0"] = {
        "selected_config": {"prevalence": round(prev_tr, 6)},
        "validation_ap": float(average_precision_score(y_va, s_b0[va_mask])),
        "validation_threshold": thr_b0,
    }

    # 2. P0: Pure Zero-Shot 3D Cofactor-State Pocket Physics (NO labels fitted!)
    raw_p0, prob_p0 = compute_zero_shot_p0_scores(X_3d, phys_cols, tr_mask)
    thr_p0 = choose_validation_threshold(y_va, prob_p0[va_mask])
    selected_scores_all["P0"] = prob_p0
    selected_thresholds["P0"] = thr_p0
    selection_log["P0"] = {
        "selected_config": {"type": "zero_shot_3d_physics_composite", "labels_used_in_ranking": False},
        "validation_ap": float(average_precision_score(y_va, prob_p0[va_mask])),
        "validation_threshold": thr_p0,
    }

    # 3. B1: v1-Style Non-Chiral 2D Morgan Logistic Classifier (C in {0.1, 1.0, 10.0})
    b1_cands = []
    for c_val in (0.1, 1.0, 10.0):
        clf = Pipeline([
            ("scaler", StandardScaler()),
            ("lr", LogisticRegression(C=c_val, l1_ratio=0.0, class_weight="balanced", solver="lbfgs", max_iter=5000, random_state=SEED)),
        ])
        clf.fit(X_ach[tr_mask], y_tr)
        p_all = safe_prob(clf.predict_proba(X_ach)[:, 1])
        ap_va = float(average_precision_score(y_va, p_all[va_mask]))
        b1_cands.append((ap_va, -c_val, c_val, p_all))
    b1_cands.sort(key=lambda x: (x[0], x[1]), reverse=True)
    best_b1 = b1_cands[0]
    thr_b1 = choose_validation_threshold(y_va, best_b1[3][va_mask])
    selected_scores_all["B1"] = best_b1[3]
    selected_thresholds["B1"] = thr_b1
    selection_log["B1"] = {
        "candidates": [{"C": c, "validation_ap": round(ap, 6)} for ap, _, c, _ in b1_cands],
        "selected_config": {"C": best_b1[2]},
        "validation_ap": best_b1[0],
        "validation_threshold": thr_b1,
    }

    # 4. M1: Chiral 2D + Physicochemical Ensemble (with 5-fold Murcko-scaffold OOF inside train!)
    pc_scaler = StandardScaler().fit(X_pc[tr_mask])
    X_2d_full = np.hstack([X_chi, pc_scaler.transform(X_pc)])

    gkf = GroupKFold(n_splits=5)
    oof_2d_lr = np.zeros(tr_mask.sum(), dtype=np.float64)
    oof_2d_xgb = np.zeros(tr_mask.sum(), dtype=np.float64)
    oof_2d_sim = np.zeros(tr_mask.sum(), dtype=np.float64)

    X_2d_tr = X_2d_full[tr_mask]
    X_chi_tr = X_chi[tr_mask]
    pos_weight = float((y_tr == 0).sum() / max(1, (y_tr == 1).sum()))

    for f_tr, f_ho in gkf.split(X_2d_tr, y_tr, groups=scaff_tr):
        lr_f = Pipeline([
            ("sc", StandardScaler()),
            ("lr", LogisticRegression(C=0.2, l1_ratio=0.0, class_weight="balanced", solver="lbfgs", max_iter=5000, random_state=SEED)),
        ]).fit(X_2d_tr[f_tr], y_tr[f_tr])
        oof_2d_lr[f_ho] = safe_prob(lr_f.predict_proba(X_2d_tr[f_ho])[:, 1])

        xgb_f = xgb.XGBClassifier(
            n_estimators=220, max_depth=4, learning_rate=0.04, subsample=0.8,
            colsample_bytree=0.4, reg_lambda=6.0, scale_pos_weight=pos_weight,
            n_jobs=8, random_state=SEED, eval_metric="logloss",
        ).fit(X_2d_tr[f_tr], y_tr[f_tr])
        oof_2d_xgb[f_ho] = safe_prob(xgb_f.predict_proba(X_2d_tr[f_ho])[:, 1])
        oof_2d_sim[f_ho] = expit(3.0 * tanimoto_max_diff(X_chi_tr[f_ho], X_chi_tr[f_tr], y_tr[f_tr]))

    # Fit full train models for M1
    lr_m1_01 = Pipeline([("sc", StandardScaler()), ("lr", LogisticRegression(C=0.1, l1_ratio=0.0, class_weight="balanced", solver="lbfgs", max_iter=5000, random_state=SEED))]).fit(X_2d_tr, y_tr)
    lr_m1_02 = Pipeline([("sc", StandardScaler()), ("lr", LogisticRegression(C=0.2, l1_ratio=0.0, class_weight="balanced", solver="lbfgs", max_iter=5000, random_state=SEED))]).fit(X_2d_tr, y_tr)
    xgb_m1 = xgb.XGBClassifier(
        n_estimators=220, max_depth=4, learning_rate=0.04, subsample=0.8,
        colsample_bytree=0.4, reg_lambda=6.0, scale_pos_weight=pos_weight,
        n_jobs=8, random_state=SEED, eval_metric="logloss",
    ).fit(X_2d_tr, y_tr)

    p_lr01 = safe_prob(lr_m1_01.predict_proba(X_2d_full)[:, 1])
    p_lr02 = safe_prob(lr_m1_02.predict_proba(X_2d_full)[:, 1])
    p_xgb = safe_prob(xgb_m1.predict_proba(X_2d_full)[:, 1])
    p_sim = safe_prob(expit(3.0 * tanimoto_max_diff(X_chi, X_chi_tr, y_tr)))
    p_m1_blend = safe_prob(0.45 * p_lr02 + 0.35 * p_xgb + 0.20 * p_sim)
    oof_m1_blend = safe_prob(0.45 * oof_2d_lr + 0.35 * oof_2d_xgb + 0.20 * oof_2d_sim)

    m1_cands = [
        ("chiral_logreg_C0.1", p_lr01, oof_2d_lr),
        ("chiral_xgb_d4", p_xgb, oof_2d_xgb),
        ("chiral_2d_ensemble", p_m1_blend, oof_m1_blend),
    ]
    m1_scored = [(float(average_precision_score(y_va, p[va_mask])), -i, name, p, oof) for i, (name, p, oof) in enumerate(m1_cands)]
    m1_scored.sort(key=lambda x: (x[0], x[1]), reverse=True)
    best_m1 = m1_scored[0]
    thr_m1 = choose_validation_threshold(y_va, best_m1[3][va_mask])
    selected_scores_all["M1"] = best_m1[3]
    selected_thresholds["M1"] = thr_m1
    selection_log["M1"] = {
        "candidates": [{"name": nm, "validation_ap": round(ap, 6)} for ap, _, nm, _, _ in m1_scored],
        "selected_config": {"name": best_m1[2]},
        "validation_ap": best_m1[0],
        "validation_threshold": thr_m1,
    }

    # 5. M2: Two-Branch Chiral 2D + 3D Multi-State Cofactor-Pocket Physics Hybrid
    X_3d_aug = np.hstack([X_3d, X_pc, raw_p0[:, None]])
    X_3d_tr = X_3d_aug[tr_mask]
    oof_3d = np.zeros(tr_mask.sum(), dtype=np.float64)

    for f_tr, f_ho in gkf.split(X_3d_tr, y_tr, groups=scaff_tr):
        lr_3d_f = Pipeline([
            ("sc", StandardScaler()),
            ("lr", LogisticRegression(C=0.15, l1_ratio=0.0, class_weight="balanced", solver="lbfgs", max_iter=5000, random_state=SEED)),
        ]).fit(X_3d_tr[f_tr], y_tr[f_tr])
        hgb_3d_f = HistGradientBoostingClassifier(
            max_iter=150, max_depth=3, min_samples_leaf=15, l2_regularization=8.0,
            class_weight="balanced", random_state=SEED,
        ).fit(X_3d_tr[f_tr], y_tr[f_tr])
        oof_3d[f_ho] = safe_prob(
            0.50 * lr_3d_f.predict_proba(X_3d_tr[f_ho])[:, 1]
            + 0.50 * hgb_3d_f.predict_proba(X_3d_tr[f_ho])[:, 1]
        )

    lr_3d_full = Pipeline([
        ("sc", StandardScaler()),
        ("lr", LogisticRegression(C=0.15, l1_ratio=0.0, class_weight="balanced", solver="lbfgs", max_iter=5000, random_state=SEED)),
    ]).fit(X_3d_tr, y_tr)
    hgb_3d_full = HistGradientBoostingClassifier(
        max_iter=150, max_depth=3, min_samples_leaf=15, l2_regularization=8.0,
        class_weight="balanced", random_state=SEED,
    ).fit(X_3d_tr, y_tr)
    prob_3d_branch = safe_prob(
        0.50 * lr_3d_full.predict_proba(X_3d_aug)[:, 1]
        + 0.50 * hgb_3d_full.predict_proba(X_3d_aug)[:, 1]
    )

    # Meta-stacker fitted STRICTLY on Train scaffold-OOF logits + zero-shot physics P0!
    Z_oof = np.column_stack([logit(best_m1[4]), logit(oof_3d), raw_p0[tr_mask]])
    meta_lr = LogisticRegression(C=0.5, l1_ratio=0.0, class_weight="balanced", solver="lbfgs", random_state=SEED).fit(Z_oof, y_tr)
    Z_all = np.column_stack([logit(best_m1[3]), logit(prob_3d_branch), raw_p0])
    p_m2_oof_stack = safe_prob(meta_lr.predict_proba(Z_all)[:, 1])
    p_m2_balanced = safe_prob(0.45 * best_m1[3] + 0.35 * prob_3d_branch + 0.20 * prob_p0)
    p_m2_phys_heavy = safe_prob(0.30 * best_m1[3] + 0.40 * prob_3d_branch + 0.30 * prob_p0)

    m2_cands = [
        ("m2_scaffold_oof_stack", p_m2_oof_stack),
        ("m2_balanced_2d_3d_physics", p_m2_balanced),
        ("m2_physics_heavy_hybrid", p_m2_phys_heavy),
    ]
    m2_scored = [(float(average_precision_score(y_va, p[va_mask])), -i, nm, p) for i, (nm, p) in enumerate(m2_cands)]
    m2_scored.sort(key=lambda x: (x[0], x[1]), reverse=True)
    best_m2 = m2_scored[0]
    thr_m2 = choose_validation_threshold(y_va, best_m2[3][va_mask])
    selected_scores_all["M2"] = best_m2[3]
    selected_thresholds["M2"] = thr_m2
    selection_log["M2"] = {
        "candidates": [{"name": nm, "validation_ap": round(ap, 6)} for ap, _, nm, _ in m2_scored],
        "selected_config": {
            "name": best_m2[2],
            "oof_meta_weights": [round(float(w), 4) for w in meta_lr.coef_[0]],
        },
        "validation_ap": best_m2[0],
        "validation_threshold": thr_m2,
    }

    val_pred_rows, test_pred_rows = [], []
    for m_name in ("B0", "P0", "B1", "M1", "M2"):
        sc_vec = selected_scores_all[m_name]
        for i, r in enumerate(struct_rows):
            rec = {
                "track": "track_a",
                "model": m_name,
                "row_id": r["row_id"],
                "compound_id": r["compound_id"],
                "assay_id": r["assay_pair_type"],
                "split": r["split"],
                "score": f"{float(sc_vec[i]):.8f}",
            }
            if r["split"] == "validation":
                val_pred_rows.append(rec)
            elif r["split"] == "test":
                test_pred_rows.append(rec)

    return selection_log, val_pred_rows, test_pred_rows


def _assay_group_sizes(assays: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Sort indices by assay_id so XGBRanker receives contiguous query groups."""
    order = np.argsort(assays, kind="mergesort")
    sorted_a = assays[order]
    _, counts = np.unique(sorted_a, return_counts=True)
    return order, counts


def _within_assay_demean(values: np.ndarray, assays: np.ndarray) -> np.ndarray:
    out = np.zeros_like(values, dtype=np.float64)
    for a in np.unique(assays):
        m = (assays == a)
        out[m] = values[m] - values[m].mean()
    return out


def run_track_b(
    root: Path,
    feat_lookup: dict[str, int],
    fp_ach: np.ndarray,
    fp_chi: np.ndarray,
    physchem: np.ndarray,
    phys3d: np.ndarray,
    phys_cols: list[str],
) -> tuple[dict[str, Any], list[dict[str, Any]], list[dict[str, Any]]]:
    struct_rows = list(csv.DictReader((root / "data/processed/track_b_structures.csv").open()))
    dev_rows = list(csv.DictReader((root / "data/processed/track_b_dev_labels.csv").open()))
    if not set(r["split"] for r in dev_rows) <= {"train", "validation"}:
        raise RuntimeError("Track B dev_labels contained non-development rows!")

    dev_map = {r["row_id"]: r for r in dev_rows}
    idx_all = np.array([feat_lookup[r["compound_id"]] for r in struct_rows], dtype=int)
    splits = np.array([r["split"] for r in struct_rows])
    tr_mask = (splits == "train")
    va_mask = (splits == "validation")

    tr_rows = [r for r in struct_rows if r["split"] == "train"]
    va_rows = [r for r in struct_rows if r["split"] == "validation"]

    y_tr = np.array([int(dev_map[r["row_id"]]["binary_y_le100nM"]) for r in tr_rows], dtype=int)
    p_tr = np.array([float(dev_map[r["row_id"]]["p_activity"]) for r in tr_rows], dtype=float)
    assay_tr = np.array([r["assay_id"] for r in tr_rows])
    scaff_tr = np.array([r["scaffold_group"] for r in tr_rows])

    y_va = np.array([int(dev_map[r["row_id"]]["binary_y_le100nM"]) for r in va_rows], dtype=int)
    p_va = np.array([float(dev_map[r["row_id"]]["p_activity"]) for r in va_rows], dtype=float)
    rel_va = np.array([dev_map[r["row_id"]]["relation"] for r in va_rows])
    assay_va = np.array([r["assay_id"] for r in va_rows])

    # Explicit assay covariates: [log10_sam_uM, is_displacement, sam_occupancy_km0.14]
    sam_uM_all = np.array([float(r["sam_uM"]) for r in struct_rows], dtype=np.float32)
    log_sam_all = np.array([float(r["log10_sam_uM"]) for r in struct_rows], dtype=np.float32)
    disp_all = np.array([float(r["is_displacement"]) for r in struct_rows], dtype=np.float32)
    sam_occ_all = sam_uM_all / (0.14 + sam_uM_all)
    X_cov = np.column_stack([log_sam_all, disp_all, sam_occ_all]).astype(np.float32)

    X_ach = fp_ach[idx_all].astype(np.float32)
    X_chi = fp_chi[idx_all].astype(np.float32)
    X_pc = physchem[idx_all].astype(np.float32)
    X_3d = phys3d[idx_all].astype(np.float32)

    selection_log: dict[str, Any] = {}
    selected_scores_all: dict[str, np.ndarray] = {}
    selected_thresholds: dict[str, float] = {}

    def val_composite(score_vec: np.ndarray) -> tuple[float, float, float]:
        s_va = score_vec[va_mask]
        ap = float(average_precision_score(y_va, s_va))
        ci_res = within_group_weighted_c_index(assay_va, p_va, rel_va, s_va)
        ci = float(ci_res["weighted_c_index"])
        return 0.5 * ci + 0.5 * ap, ci, ap

    # 1. B0: Constant training prevalence
    prev_tr = float(y_tr.mean())
    s_b0 = np.full(len(struct_rows), prev_tr, dtype=np.float64)
    thr_b0 = choose_validation_threshold(y_va, s_b0[va_mask])
    _, ci_b0, ap_b0 = val_composite(s_b0)
    selected_scores_all["B0"] = s_b0
    selected_thresholds["B0"] = thr_b0
    selection_log["B0"] = {
        "selected_config": {"prevalence": round(prev_tr, 6)},
        "validation_weighted_c_index": ci_b0,
        "validation_pooled_ap": ap_b0,
        "validation_threshold": thr_b0,
    }

    # 2. B1: v1-Style Uncalibrated Non-Chiral 2D Classifier (threshold locked on high-prevalence validation!)
    b1_cands = []
    for c_val in (0.1, 1.0, 10.0):
        clf = Pipeline([
            ("sc", StandardScaler()),
            ("lr", LogisticRegression(C=c_val, l1_ratio=0.0, class_weight="balanced", solver="lbfgs", max_iter=5000, random_state=SEED)),
        ]).fit(X_ach[tr_mask], y_tr)
        p_all = safe_prob(clf.predict_proba(X_ach)[:, 1])
        ap_va = float(average_precision_score(y_va, p_all[va_mask]))
        _, ci_va, _ = val_composite(p_all)
        b1_cands.append((ap_va, -c_val, c_val, ci_va, p_all))
    b1_cands.sort(key=lambda x: (x[0], x[1]), reverse=True)
    best_b1 = b1_cands[0]
    thr_b1 = choose_validation_threshold(y_va, best_b1[4][va_mask])
    selected_scores_all["B1"] = best_b1[4]
    selected_thresholds["B1"] = thr_b1
    selection_log["B1"] = {
        "candidates": [{"C": c, "validation_ap": round(ap, 6), "validation_c_index": round(ci, 6)} for ap, _, c, ci, _ in b1_cands],
        "selected_config": {"C": best_b1[2], "threshold_policy": "v1_uncalibrated_validation_balanced_accuracy"},
        "validation_weighted_c_index": best_b1[3],
        "validation_pooled_ap": best_b1[0],
        "validation_threshold": thr_b1,
    }

    # 3. M1: Chiral 2D + Physicochemical + Assay-Covariate Within-Source Pairwise Ranker
    dense_2d_scaler = StandardScaler().fit(np.hstack([X_pc[tr_mask], X_cov[tr_mask]]))
    X_m1_all = np.hstack([X_chi, dense_2d_scaler.transform(np.hstack([X_pc, X_cov]))])
    X_m1_tr = X_m1_all[tr_mask]

    p_tr_demeaned = _within_assay_demean(p_tr, assay_tr)
    order_tr, groups_tr = _assay_group_sizes(assay_tr)

    # Fit within-assay XGBRanker (pairwise loss strictly within each assay!)
    ranker_m1 = xgb.XGBRanker(
        objective="rank:pairwise", n_estimators=220, max_depth=4, learning_rate=0.04,
        subsample=0.8, colsample_bytree=0.35, reg_lambda=8.0, n_jobs=8, random_state=SEED,
    )
    ranker_m1.fit(X_m1_tr[order_tr], p_tr[order_tr], group=groups_tr)
    raw_xgb_m1 = ranker_m1.predict(X_m1_all)

    # Fit within-assay demeaned Ridge potency ranker + chiral assay-conditioned logistic
    ridge_m1 = Pipeline([("sc", StandardScaler()), ("rg", Ridge(alpha=600.0, random_state=SEED))]).fit(X_m1_tr, p_tr_demeaned)
    raw_ridge_m1 = ridge_m1.predict(X_m1_all)

    cls_m1 = Pipeline([
        ("sc", StandardScaler()),
        ("lr", LogisticRegression(C=0.15, l1_ratio=0.0, class_weight="balanced", solver="lbfgs", max_iter=5000, random_state=SEED)),
    ]).fit(X_m1_tr, y_tr)
    raw_cls_m1 = logit(safe_prob(cls_m1.predict_proba(X_m1_all)[:, 1]))

    sim_m1 = tanimoto_max_diff(X_chi, X_chi[tr_mask], y_tr)

    def calibrate_rank_score(raw_rank: np.ndarray, oof_rank_tr: np.ndarray) -> tuple[np.ndarray, float]:
        """Assay-covariate conditioned probability calibration anchored to train prevalence + [SAM] shift."""
        cal_X_tr = np.column_stack([oof_rank_tr, X_cov[tr_mask]])
        cal_X_all = np.column_stack([raw_rank, X_cov])
        cal_scaler = StandardScaler().fit(cal_X_tr)
        # Fit natural-prevalence calibrator on OOF rank score + assay covariates (avoids validation prevalence overfit!)
        cal = LogisticRegression(C=1.0, l1_ratio=0.0, solver="lbfgs", random_state=SEED).fit(cal_scaler.transform(cal_X_tr), y_tr)
        prob_all = safe_prob(cal.predict_proba(cal_scaler.transform(cal_X_all))[:, 1])
        # Natural Bayes decision boundary for prevalence-calibrated probabilities: train base rate!
        return prob_all, float(y_tr.mean())

    # Compute 5-fold Murcko-scaffold OOF rank scores inside train for leak-free assay calibration & M2 stacking!
    gkf = GroupKFold(n_splits=5)
    oof_xgb_m1 = np.zeros(len(tr_rows), dtype=np.float64)
    oof_ridge_m1 = np.zeros(len(tr_rows), dtype=np.float64)
    oof_cls_m1 = np.zeros(len(tr_rows), dtype=np.float64)
    oof_sim_m1 = np.zeros(len(tr_rows), dtype=np.float64)

    for f_tr, f_ho in gkf.split(X_m1_tr, y_tr, groups=scaff_tr):
        ord_f, grp_f = _assay_group_sizes(assay_tr[f_tr])
        rk_f = xgb.XGBRanker(
            objective="rank:pairwise", n_estimators=220, max_depth=4, learning_rate=0.04,
            subsample=0.8, colsample_bytree=0.35, reg_lambda=8.0, n_jobs=8, random_state=SEED,
        ).fit(X_m1_tr[f_tr][ord_f], p_tr[f_tr][ord_f], group=grp_f)
        oof_xgb_m1[f_ho] = rk_f.predict(X_m1_tr[f_ho])

        rg_f = Pipeline([("sc", StandardScaler()), ("rg", Ridge(alpha=600.0, random_state=SEED))]).fit(
            X_m1_tr[f_tr], _within_assay_demean(p_tr[f_tr], assay_tr[f_tr])
        )
        oof_ridge_m1[f_ho] = rg_f.predict(X_m1_tr[f_ho])

        cl_f = Pipeline([
            ("sc", StandardScaler()),
            ("lr", LogisticRegression(C=0.15, l1_ratio=0.0, class_weight="balanced", solver="lbfgs", max_iter=5000, random_state=SEED)),
        ]).fit(X_m1_tr[f_tr], y_tr[f_tr])
        oof_cls_m1[f_ho] = logit(safe_prob(cl_f.predict_proba(X_m1_tr[f_ho])[:, 1]))
        oof_sim_m1[f_ho] = tanimoto_max_diff(X_chi[tr_mask][f_ho], X_chi[tr_mask][f_tr], y_tr[f_tr])

    def z_norm(v_all: np.ndarray, v_oof: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        mu, sd = float(v_oof.mean()), float(max(1e-6, v_oof.std()))
        return (v_all - mu) / sd, (v_oof - mu) / sd

    z_xgb_all, z_xgb_oof = z_norm(raw_xgb_m1, oof_xgb_m1)
    z_rg_all, z_rg_oof = z_norm(raw_ridge_m1, oof_ridge_m1)
    z_cl_all, z_cl_oof = z_norm(raw_cls_m1, oof_cls_m1)
    z_sm_all, z_sm_oof = z_norm(sim_m1, oof_sim_m1)

    ens_m1_all = 0.35 * z_xgb_all + 0.25 * z_rg_all + 0.25 * z_cl_all + 0.15 * z_sm_all
    ens_m1_oof = 0.35 * z_xgb_oof + 0.25 * z_rg_oof + 0.25 * z_cl_oof + 0.15 * z_sm_oof

    m1_raw_cands = [
        ("m1_pairwise_xgbranker", z_xgb_all, z_xgb_oof),
        ("m1_demeaned_ridge_plus_cls", 0.5 * z_rg_all + 0.5 * z_cl_all, 0.5 * z_rg_oof + 0.5 * z_cl_oof),
        ("m1_chiral_pairwise_ensemble", ens_m1_all, ens_m1_oof),
    ]
    m1_scored = []
    for i, (nm, r_all, r_oof) in enumerate(m1_raw_cands):
        prob_vec, thr_cal = calibrate_rank_score(r_all, r_oof)
        comp, ci, ap = val_composite(prob_vec)
        m1_scored.append((comp, -i, nm, ci, ap, prob_vec, thr_cal, r_all, r_oof))
    m1_scored.sort(key=lambda x: (x[0], x[1]), reverse=True)
    best_m1 = m1_scored[0]
    selected_scores_all["M1"] = best_m1[5]
    selected_thresholds["M1"] = best_m1[6]
    selection_log["M1"] = {
        "candidates": [
            {"name": nm, "val_composite": round(comp, 6), "val_c_index": round(ci, 6), "val_ap": round(ap, 6)}
            for comp, _, nm, ci, ap, _, _, _, _ in m1_scored
        ],
        "selected_config": {"name": best_m1[2], "calibration": "assay_covariate_oof_logistic_base_rate"},
        "validation_weighted_c_index": best_m1[3],
        "validation_pooled_ap": best_m1[4],
        "validation_threshold": best_m1[6],
    }

    # 4. M2: Chiral 2D + 3D Multi-State Crystal Pocket Physics + Within-Assay Pairwise Ranker
    col_idx = {c: i for i, c in enumerate(phys_cols)}
    # A priori 3D multi-state potency affinity proxy (lower Vina energy & higher O3A/shape match => higher potency!)
    sc_3d_unlab = StandardScaler().fit(X_3d[tr_mask])
    z3d = sc_3d_unlab.transform(X_3d)
    prior_3d_potency = (
        -1.0 * z3d[:, col_idx["vina_e_min_3state"]]
        - 0.5 * z3d[:, col_idx["vina_e_sam_4x61"]]
        - 0.5 * z3d[:, col_idx["vina_e_nuc_6ckc"]]
        + 0.8 * np.maximum.reduce([
            z3d[:, col_idx["o3a_norm_s1_sam"]],
            z3d[:, col_idx["o3a_norm_s2_mta"]],
            z3d[:, col_idx["o3a_norm_s3_nuc"]],
        ])
        + 0.8 * np.maximum.reduce([
            z3d[:, col_idx["shape_sim_s1_sam"]],
            z3d[:, col_idx["shape_sim_s2_mta"]],
            z3d[:, col_idx["shape_sim_s3_nuc"]],
        ])
    ) / 3.0

    X_3d_branch_all = np.hstack([X_3d, X_pc, X_cov, prior_3d_potency[:, None]])
    X_3d_branch_tr = X_3d_branch_all[tr_mask]

    oof_3d_rank = np.zeros(len(tr_rows), dtype=np.float64)
    for f_tr, f_ho in gkf.split(X_3d_branch_tr, y_tr, groups=scaff_tr):
        ord_f, grp_f = _assay_group_sizes(assay_tr[f_tr])
        rk3_f = xgb.XGBRanker(
            objective="rank:pairwise", n_estimators=180, max_depth=3, learning_rate=0.04,
            subsample=0.8, colsample_bytree=0.7, reg_lambda=10.0, n_jobs=8, random_state=SEED,
        ).fit(X_3d_branch_tr[f_tr][ord_f], p_tr[f_tr][ord_f], group=grp_f)
        hgb3_f = HistGradientBoostingRegressor(
            max_iter=140, max_depth=3, min_samples_leaf=18, l2_regularization=10.0, random_state=SEED,
        ).fit(X_3d_branch_tr[f_tr], _within_assay_demean(p_tr[f_tr], assay_tr[f_tr]))
        oof_3d_rank[f_ho] = 0.5 * rk3_f.predict(X_3d_branch_tr[f_ho]) + 0.5 * hgb3_f.predict(X_3d_branch_tr[f_ho])

    rk3_full = xgb.XGBRanker(
        objective="rank:pairwise", n_estimators=180, max_depth=3, learning_rate=0.04,
        subsample=0.8, colsample_bytree=0.7, reg_lambda=10.0, n_jobs=8, random_state=SEED,
    ).fit(X_3d_branch_tr[order_tr], p_tr[order_tr], group=groups_tr)
    hgb3_full = HistGradientBoostingRegressor(
        max_iter=140, max_depth=3, min_samples_leaf=18, l2_regularization=10.0, random_state=SEED,
    ).fit(X_3d_branch_tr, p_tr_demeaned)
    raw_3d_rank_all = 0.5 * rk3_full.predict(X_3d_branch_all) + 0.5 * hgb3_full.predict(X_3d_branch_all)

    # Joint Chiral 2D + 3D Crystal Pocket Physics + Assay Covariates Pairwise XGBRanker
    dense_m2_scaler = StandardScaler().fit(np.hstack([X_3d[tr_mask], X_pc[tr_mask], X_cov[tr_mask], prior_3d_potency[tr_mask, None]]))
    X_m2_joint_all = np.hstack([X_chi, dense_m2_scaler.transform(np.hstack([X_3d, X_pc, X_cov, prior_3d_potency[:, None]]))])
    X_m2_joint_tr = X_m2_joint_all[tr_mask]

    oof_joint_m2 = np.zeros(len(tr_rows), dtype=np.float64)
    for f_tr, f_ho in gkf.split(X_m2_joint_tr, y_tr, groups=scaff_tr):
        ord_f, grp_f = _assay_group_sizes(assay_tr[f_tr])
        rkj_f = xgb.XGBRanker(
            objective="rank:pairwise", n_estimators=240, max_depth=4, learning_rate=0.04,
            subsample=0.8, colsample_bytree=0.35, reg_lambda=8.0, n_jobs=8, random_state=SEED,
        ).fit(X_m2_joint_tr[f_tr][ord_f], p_tr[f_tr][ord_f], group=grp_f)
        oof_joint_m2[f_ho] = rkj_f.predict(X_m2_joint_tr[f_ho])

    rkj_full = xgb.XGBRanker(
        objective="rank:pairwise", n_estimators=240, max_depth=4, learning_rate=0.04,
        subsample=0.8, colsample_bytree=0.35, reg_lambda=8.0, n_jobs=8, random_state=SEED,
    ).fit(X_m2_joint_tr[order_tr], p_tr[order_tr], group=groups_tr)
    raw_joint_m2 = rkj_full.predict(X_m2_joint_all)

    z_j_all, z_j_oof = z_norm(raw_joint_m2, oof_joint_m2)
    z_3d_all, z_3d_oof = z_norm(raw_3d_rank_all, oof_3d_rank)
    z_pr_all, z_pr_oof = z_norm(prior_3d_potency, prior_3d_potency[tr_mask])

    # Meta-ranker fitted on Train scaffold-OOF (M1 ranker, Joint 2D+3D ranker, 3D branch, 3D prior)
    Z_tb_oof = np.column_stack([best_m1[8], z_j_oof, z_3d_oof, z_pr_oof])
    Z_tb_all = np.column_stack([best_m1[7], z_j_all, z_3d_all, z_pr_all])
    meta_rg = Ridge(alpha=50.0, positive=True).fit(Z_tb_oof, p_tr_demeaned)
    if meta_rg.coef_.sum() <= 1e-6:
        w_meta = np.array([0.40, 0.35, 0.15, 0.10])
    else:
        w_meta = meta_rg.coef_ / meta_rg.coef_.sum()

    r_m2_oof_meta_all = Z_tb_all @ w_meta
    r_m2_oof_meta_oof = Z_tb_oof @ w_meta

    r_m2_joint_blend_all = 0.50 * z_j_all + 0.30 * best_m1[7] + 0.20 * z_3d_all
    r_m2_joint_blend_oof = 0.50 * z_j_oof + 0.30 * best_m1[8] + 0.20 * z_3d_oof

    m2_raw_cands = [
        ("m2_joint_2d_3d_xgbranker", z_j_all, z_j_oof),
        ("m2_joint_plus_3d_branch_ensemble", r_m2_joint_blend_all, r_m2_joint_blend_oof),
        ("m2_scaffold_oof_positive_ridge", r_m2_oof_meta_all, r_m2_oof_meta_oof),
    ]
    m2_scored = []
    for i, (nm, r_all, r_oof) in enumerate(m2_raw_cands):
        prob_vec, thr_cal = calibrate_rank_score(r_all, r_oof)
        comp, ci, ap = val_composite(prob_vec)
        m2_scored.append((comp, -i, nm, ci, ap, prob_vec, thr_cal))
    m2_scored.sort(key=lambda x: (x[0], x[1]), reverse=True)
    best_m2 = m2_scored[0]
    selected_scores_all["M2"] = best_m2[5]
    selected_thresholds["M2"] = best_m2[6]
    selection_log["M2"] = {
        "candidates": [
            {"name": nm, "val_composite": round(comp, 6), "val_c_index": round(ci, 6), "val_ap": round(ap, 6)}
            for comp, _, nm, ci, ap, _, _ in m2_scored
        ],
        "selected_config": {
            "name": best_m2[2],
            "oof_positive_weights": [round(float(x), 4) for x in w_meta],
            "calibration": "assay_covariate_oof_logistic_base_rate",
        },
        "validation_weighted_c_index": best_m2[3],
        "validation_pooled_ap": best_m2[4],
        "validation_threshold": best_m2[6],
    }

    val_pred_rows, test_pred_rows = [], []
    for m_name in ("B0", "B1", "M1", "M2"):
        sc_vec = selected_scores_all[m_name]
        for i, r in enumerate(struct_rows):
            rec = {
                "track": "track_b",
                "model": m_name,
                "row_id": r["row_id"],
                "compound_id": r["compound_id"],
                "assay_id": r["assay_id"],
                "split": r["split"],
                "score": f"{float(sc_vec[i]):.8f}",
            }
            if r["split"] == "validation":
                val_pred_rows.append(rec)
            elif r["split"] == "test":
                test_pred_rows.append(rec)

    return selection_log, val_pred_rows, test_pred_rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    root = args.root.resolve()

    freeze, feat_meta = verify_inputs(root)
    npz = np.load(root / feat_meta["npz_path"], allow_pickle=False)
    cids = [str(x) for x in npz["compound_ids"]]
    feat_lookup = {cid: i for i, cid in enumerate(cids)}
    fp_ach = npz["fp_achiral"]
    fp_chi = npz["fp_chiral"]
    physchem = npz["physchem"]
    phys3d = npz["phys3d"]
    phys_cols = [str(x) for x in npz["phys3d_cols"]]

    sel_a, val_a, test_a = run_track_a(root, feat_lookup, fp_ach, fp_chi, physchem, phys3d, phys_cols)
    sel_b, val_b, test_b = run_track_b(root, feat_lookup, fp_ach, fp_chi, physchem, phys3d, phys_cols)

    pred_dir = root / "data/predictions"
    pred_dir.mkdir(parents=True, exist_ok=True)
    sealed_path = pred_dir / "sealed_predictions.csv"
    val_path = pred_dir / "validation_predictions.csv"

    fieldnames = ["track", "model", "row_id", "compound_id", "assay_id", "split", "score"]
    for path, rows in ((sealed_path, test_a + test_b), (val_path, val_a + val_b)):
        with path.open("w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=fieldnames)
            w.writeheader()
            w.writerows(rows)

    lock = {
        "protocol_id": freeze["bound_sha256"] and "PRMT5-V2-3D-COOPERATIVITY-2026",
        "protocol_freeze_sha256": sha256_file(root / "receipts/protocol_freeze.json"),
        "feature_manifest_sha256": sha256_file(root / "manifests/feature_manifest_v2.json"),
        "test_labels_accessed": False,
        "test_metrics_computed": False,
        "track_a": sel_a,
        "track_b": sel_b,
        "sealed_predictions_path": str(sealed_path.relative_to(root)),
        "sealed_predictions_sha256": sha256_file(sealed_path),
        "sealed_prediction_rows": len(test_a) + len(test_b),
        "validation_predictions_path": str(val_path.relative_to(root)),
        "validation_predictions_sha256": sha256_file(val_path),
    }
    lock_path = root / "manifests/model_selection_lock.json"
    lock_path.write_text(json.dumps(lock, indent=2) + "\n", encoding="utf-8")

    seal_receipt = {
        "stage": "STAGE_2_PREDICTION_SEAL",
        "status": "SEALED_BEFORE_TEST_EVALUATION",
        "sealed_at_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "protocol_freeze_sha256": lock["protocol_freeze_sha256"],
        "feature_manifest_sha256": lock["feature_manifest_sha256"],
        "model_selection_lock_sha256": sha256_file(lock_path),
        "sealed_predictions_sha256": lock["sealed_predictions_sha256"],
        "validation_predictions_sha256": lock["validation_predictions_sha256"],
        "code_sha256": {
            "scripts/featurize_3d.py": sha256_file(root / "scripts/featurize_3d.py"),
            "scripts/metrics_v2.py": sha256_file(root / "scripts/metrics_v2.py"),
            "scripts/train_and_seal.py": sha256_file(root / "scripts/train_and_seal.py"),
        },
        "test_labels_accessed": False,
        "test_metrics_computed": False,
        "track_a_test_rows": len(test_a),
        "track_b_test_rows": len(test_b),
    }
    seal_path = root / "receipts/prediction_seal.json"
    seal_path.write_text(json.dumps(seal_receipt, indent=2) + "\n", encoding="utf-8")

    print("=== VALIDATION SELECTION SUMMARY (TEST UNTOUCHED) ===")
    print("Track A Validation AP:", {m: round(sel_a[m]["validation_ap"], 4) for m in sel_a})
    print("Track B Validation C-index / AP:", {
        m: (round(sel_b[m]["validation_weighted_c_index"], 4), round(sel_b[m]["validation_pooled_ap"], 4))
        for m in sel_b
    })
    print(f"Sealed {len(test_a) + len(test_b)} test predictions -> {seal_path}")


if __name__ == "__main__":
    main()
