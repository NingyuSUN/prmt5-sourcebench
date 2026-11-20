#!/usr/bin/env python3
"""Stage 3: Score sealed PRMT5 v2 predictions, run 2,000 paired scaffold bootstraps, and audit biophysical drivers.

Strictly enforces zero-leakage replay contracts:
1. Verifies `receipts/protocol_freeze.json`, `receipts/prediction_seal.json`, and `manifests/model_selection_lock.json` SHA-256 hashes.
2. Refuses to overwrite existing output files (use `--output-dir` for deterministic replay).
3. Never fits, tunes, or rescales any model or threshold.
"""
from __future__ import annotations

import argparse
import collections
import csv
import hashlib
import json
import math
from pathlib import Path
import sys
from typing import Any

import numpy as np
from scipy.stats import spearmanr
from sklearn.metrics import average_precision_score, roc_auc_score

sys.path.insert(0, str(Path(__file__).resolve().parent))
from metrics_v2 import (
    binary_threshold_metrics,
    censor_aware_concordance_index,
    paired_scaffold_bootstrap,
    ranking_summary,
    within_group_weighted_c_index,
)


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def require(cond: bool, msg: str) -> None:
    if not cond:
        raise ValueError(msg)


def verify_full_seal_chain(root: Path) -> tuple[dict, dict, dict]:
    freeze_path = root / "receipts/protocol_freeze.json"
    seal_path = root / "receipts/prediction_seal.json"
    lock_path = root / "manifests/model_selection_lock.json"
    feat_path = root / "manifests/feature_manifest_v2.json"

    for p in (freeze_path, seal_path, lock_path, feat_path):
        require(p.exists(), f"Missing required manifest/receipt: {p}")

    freeze = json.loads(freeze_path.read_text(encoding="utf-8"))
    seal = json.loads(seal_path.read_text(encoding="utf-8"))
    lock = json.loads(lock_path.read_text(encoding="utf-8"))
    feat = json.loads(feat_path.read_text(encoding="utf-8"))

    require(freeze.get("models_fitted") is False, "protocol_freeze violated models_fitted=False")
    require(seal.get("test_labels_accessed") is False, "prediction_seal violated test_labels_accessed=False")
    require(seal.get("test_metrics_computed") is False, "prediction_seal violated test_metrics_computed=False")
    require(lock.get("test_labels_accessed") is False, "model_selection_lock violated test_labels_accessed=False")

    require(seal["protocol_freeze_sha256"] == sha256_file(freeze_path), "protocol_freeze hash mismatch in prediction_seal")
    require(lock["protocol_freeze_sha256"] == sha256_file(freeze_path), "protocol_freeze hash mismatch in model_selection_lock")
    require(seal["feature_manifest_sha256"] == sha256_file(feat_path), "feature_manifest hash mismatch in prediction_seal")
    require(seal["model_selection_lock_sha256"] == sha256_file(lock_path), "model_selection_lock hash mismatch in prediction_seal")

    for rel, expected in freeze["bound_sha256"].items():
        actual = sha256_file(root / rel)
        require(actual == expected, f"Tampered frozen artifact {rel}: {actual} != {expected}")

    require(sha256_file(root / feat["npz_path"]) == feat["npz_sha256"], "Tampered feature NPZ bytes")
    require(sha256_file(root / feat["csv_path"]) == feat["csv_sha256"], "Tampered 3D physics CSV bytes")

    pred_file = root / lock["sealed_predictions_path"]
    require(sha256_file(pred_file) == lock["sealed_predictions_sha256"], "Tampered sealed_predictions.csv vs lock")
    require(sha256_file(pred_file) == seal["sealed_predictions_sha256"], "Tampered sealed_predictions.csv vs seal")

    for script_rel, expected_sha in seal["code_sha256"].items():
        require(sha256_file(root / script_rel) == expected_sha, f"Tampered Stage 2 script {script_rel}")

    return freeze, seal, lock


def load_predictions_for_track(
    pred_rows: list[dict[str, str]],
    track: str,
    expected_row_ids: list[str],
    expected_models: tuple[str, ...],
) -> dict[str, np.ndarray]:
    expected_set = set(expected_row_ids)
    require(len(expected_set) == len(expected_row_ids), f"Duplicate row_ids in {track}")
    by_model: dict[str, dict[str, float]] = collections.defaultdict(dict)

    for r in pred_rows:
        if r["track"] != track:
            continue
        require(r["split"] == "test", f"Non-test row in sealed_predictions.csv: {r}")
        m = r["model"]
        rid = r["row_id"]
        require(m in expected_models, f"Unexpected model {m} in {track}")
        require(rid in expected_set, f"Unknown row_id {rid} in {track}")
        require(rid not in by_model[m], f"Duplicate prediction for ({track}, {m}, {rid})")
        sc = float(r["score"])
        require(math.isfinite(sc) and 0.0 <= sc <= 1.0, f"Invalid score {sc} for ({track}, {m}, {rid})")
        by_model[m][rid] = sc

    require(set(by_model) == set(expected_models), f"Missing models in {track}: {set(expected_models) - set(by_model)}")
    for m in expected_models:
        require(set(by_model[m]) == expected_set, f"Incomplete coverage for ({track}, {m})")

    return {m: np.array([by_model[m][rid] for rid in expected_row_ids], dtype=np.float64) for m in expected_models}


def evaluate_all(root: Path) -> tuple[dict[str, Any], dict[str, Any]]:
    freeze, seal, lock = verify_full_seal_chain(root)

    pred_rows = list(csv.DictReader((root / lock["sealed_predictions_path"]).open()))
    require(len(pred_rows) == lock["sealed_prediction_rows"], "Prediction row count mismatch")

    # =========================================================================
    # TRACK A EVALUATION
    # =========================================================================
    ta_struct = [r for r in csv.DictReader((root / "data/processed/track_a_structures.csv").open()) if r["split"] == "test"]
    ta_labels_map = {r["row_id"]: r for r in csv.DictReader((root / "data/processed/sealed_test_labels/track_a_test_labels.csv").open())}
    require(set(r["row_id"] for r in ta_struct) == set(ta_labels_map), "Track A test structure/label mismatch")

    ta_rids = [r["row_id"] for r in ta_struct]
    ta_cids = [r["compound_id"] for r in ta_struct]
    ta_scaff = np.array([r["scaffold_group"] for r in ta_struct])
    ta_fams = np.array([r["source_family"] for r in ta_struct])
    ta_y3 = np.array([int(ta_labels_map[rid]["y_coop_3x"]) for rid in ta_rids], dtype=int)
    ta_y10 = np.array([int(ta_labels_map[rid]["y_coop_10x"]) for rid in ta_rids], dtype=int)
    ta_dp = np.array([float(ta_labels_map[rid]["delta_p_activity"]) for rid in ta_rids], dtype=float)
    # Censoring relation for delta_p = log10(SAM / MTA): exact unless censored against direction
    ta_dp_rel = []
    for rid in ta_rids:
        lr = ta_labels_map[rid]
        rm, rs = lr["mta_rel"], lr["sam_rel"]
        if rm == "=" and rs == "=":
            ta_dp_rel.append("=")
        elif rm in ("=", "<", "<=") and rs in ("=", ">", ">="):
            ta_dp_rel.append("<")  # lower-bounded upward in delta_p (true delta_p >= observed)
        elif rm in ("=", ">", ">=") and rs in ("=", "<", "<="):
            ta_dp_rel.append(">")  # upper-bounded downward in delta_p (true delta_p <= observed)
        else:
            ta_dp_rel.append("ambig")
    ta_dp_rel_arr = np.array(ta_dp_rel)

    ta_models = ("B0", "P0", "B1", "M1", "M2")
    ta_scores = load_predictions_for_track(pred_rows, "track_a", ta_rids, ta_models)

    ta_metrics = {}
    for m in ta_models:
        s = ta_scores[m]
        thr = float(lock["track_a"][m]["validation_threshold"])
        rho_res = spearmanr(ta_dp, s)
        rho_val = float(rho_res.statistic) if np.isfinite(rho_res.statistic) else 0.0
        ci_pooled, npairs_pooled = censor_aware_concordance_index(ta_dp, ta_dp_rel_arr, s, min_diff=0.10)
        within_study = within_group_weighted_c_index(ta_fams, ta_dp, ta_dp_rel_arr, s)
        per_source_ap = {}
        for fam in sorted(set(ta_fams)):
            mask = (ta_fams == fam)
            if 0 < ta_y3[mask].sum() < mask.sum():
                per_source_ap[fam] = {
                    "n": int(mask.sum()),
                    "pos_3x": int(ta_y3[mask].sum()),
                    "ap_3x": float(average_precision_score(ta_y3[mask], s[mask])),
                    "roc_auc_3x": float(roc_auc_score(ta_y3[mask], s[mask])),
                }
        ta_metrics[m] = {
            "primary_y_coop_3x": ranking_summary(ta_y3, s, ta_cids, top_k=20),
            "stringent_y_coop_10x": ranking_summary(ta_y10, s, ta_cids, top_k=20),
            "threshold_y_coop_3x": binary_threshold_metrics(ta_y3, s, thr),
            "continuous_delta_p": {
                "spearman_rho_pooled": rho_val,
                "pooled_c_index": ci_pooled,
                "pooled_comparable_pairs": npairs_pooled,
                "within_study_c_index": within_study["weighted_c_index"],
                "within_study_breakdown": within_study["per_group"],
            },
            "per_source_y_coop_3x": per_source_ap,
        }

    ta_comparisons = [("M2", "B1"), ("M2", "M1"), ("M1", "B1"), ("P0", "B1"), ("M2", "P0")]
    ta_boot_ap3 = paired_scaffold_bootstrap(
        ta_scaff, ta_scores,
        lambda idx, s_sub: float(average_precision_score(ta_y3[idx], s_sub)) if len(np.unique(ta_y3[idx])) == 2 else None,
        ta_comparisons, n_draws=2000, seed=20261120,
    )
    ta_boot_ap10 = paired_scaffold_bootstrap(
        ta_scaff, ta_scores,
        lambda idx, s_sub: float(average_precision_score(ta_y10[idx], s_sub)) if len(np.unique(ta_y10[idx])) == 2 else None,
        ta_comparisons, n_draws=2000, seed=20261120,
    )
    ta_boot_auc3 = paired_scaffold_bootstrap(
        ta_scaff, {k: v for k, v in ta_scores.items() if k != "B0"},
        lambda idx, s_sub: float(roc_auc_score(ta_y3[idx], s_sub)) if len(np.unique(ta_y3[idx])) == 2 else None,
        ta_comparisons, n_draws=2000, seed=20261120,
    )
    ta_boot_cidx = paired_scaffold_bootstrap(
        ta_scaff, {k: v for k, v in ta_scores.items() if k != "B0"},
        lambda idx, s_sub: censor_aware_concordance_index(ta_dp[idx], ta_dp_rel_arr[idx], s_sub, min_diff=0.10)[0],
        ta_comparisons, n_draws=2000, seed=20261120,
    )

    # =========================================================================
    # TRACK B EVALUATION
    # =========================================================================
    tb_struct = [r for r in csv.DictReader((root / "data/processed/track_b_structures.csv").open()) if r["split"] == "test"]
    tb_labels_map = {r["row_id"]: r for r in csv.DictReader((root / "data/processed/sealed_test_labels/track_b_test_labels.csv").open())}
    require(set(r["row_id"] for r in tb_struct) == set(tb_labels_map), "Track B test structure/label mismatch")

    tb_rids = [r["row_id"] for r in tb_struct]
    tb_cids = [r["compound_id"] for r in tb_struct]
    tb_assays = np.array([r["assay_id"] for r in tb_struct])
    tb_scaff = np.array([r["scaffold_group"] for r in tb_struct])
    tb_chiral = np.array([r["has_defined_chirality"] == "True" for r in tb_struct])
    tb_y = np.array([int(tb_labels_map[rid]["binary_y_le100nM"]) for rid in tb_rids], dtype=int)
    tb_p = np.array([float(tb_labels_map[rid]["p_activity"]) for rid in tb_rids], dtype=float)
    tb_rel = np.array([tb_labels_map[rid]["relation"] for rid in tb_rids])
    tb_weak = np.array([int(tb_labels_map[rid]["is_weak_gt1000nM"]) for rid in tb_rids], dtype=int)

    strict_vs_mask = (tb_y == 1) | (tb_weak == 1)
    tb_models = ("B0", "B1", "M1", "M2")
    tb_scores = load_predictions_for_track(pred_rows, "track_b", tb_rids, tb_models)

    tb_metrics = {}
    for m in tb_models:
        s = tb_scores[m]
        thr = float(lock["track_b"][m]["validation_threshold"])
        within_c = within_group_weighted_c_index(tb_assays, tb_p, tb_rel, s)
        strict_vs = ranking_summary(
            tb_y[strict_vs_mask], s[strict_vs_mask],
            [tb_cids[i] for i in np.flatnonzero(strict_vs_mask)], top_k=25,
        )
        thr_m = binary_threshold_metrics(tb_y, s, thr)
        thr_m["predicted_positive_rate"] = float((s >= thr).mean())
        thr_m["true_positive_rate_prevalence"] = float(tb_y.mean())

        per_assay_bin = {}
        for aid in sorted(set(tb_assays)):
            am = (tb_assays == aid)
            pos_a = int(tb_y[am].sum())
            per_assay_bin[aid] = {
                "n": int(am.sum()),
                "positives_le100nM": pos_a,
                "weak_gt1000nM": int(tb_weak[am].sum()),
                "ap": float(average_precision_score(tb_y[am], s[am])) if 0 < pos_a < am.sum() else None,
                "roc_auc": float(roc_auc_score(tb_y[am], s[am])) if 0 < pos_a < am.sum() else None,
                "c_index": within_c["per_group"][aid]["c_index"],
                "spearman_rho_exact": within_c["per_group"][aid]["spearman_rho_exact"],
            }

        tb_metrics[m] = {
            "pooled_ranking": ranking_summary(tb_y, s, tb_cids, top_k=25),
            "within_assay_concordance": within_c,
            "strict_active_le100nM_vs_weak_gt1000nM": strict_vs,
            "threshold_calibration": thr_m,
            "chiral_subset_ranking": ranking_summary(
                tb_y[tb_chiral], s[tb_chiral],
                [tb_cids[i] for i in np.flatnonzero(tb_chiral)], top_k=20,
            ),
            "per_assay": per_assay_bin,
        }

    tb_comparisons = [("M2", "B1"), ("M2", "M1"), ("M1", "B1")]
    tb_boot_ap = paired_scaffold_bootstrap(
        tb_scaff, tb_scores,
        lambda idx, s_sub: float(average_precision_score(tb_y[idx], s_sub)) if len(np.unique(tb_y[idx])) == 2 else None,
        tb_comparisons, n_draws=2000, seed=20261120,
    )
    tb_boot_auc = paired_scaffold_bootstrap(
        tb_scaff, {k: v for k, v in tb_scores.items() if k != "B0"},
        lambda idx, s_sub: float(roc_auc_score(tb_y[idx], s_sub)) if len(np.unique(tb_y[idx])) == 2 else None,
        tb_comparisons, n_draws=2000, seed=20261120,
    )
    tb_boot_cidx = paired_scaffold_bootstrap(
        tb_scaff, {k: v for k, v in tb_scores.items() if k != "B0"},
        lambda idx, s_sub: within_group_weighted_c_index(tb_assays[idx], tb_p[idx], tb_rel[idx], s_sub)["weighted_c_index"],
        tb_comparisons, n_draws=2000, seed=20261120,
    )
    tb_boot_strict_ap = paired_scaffold_bootstrap(
        tb_scaff[strict_vs_mask], {k: v[strict_vs_mask] for k, v in tb_scores.items()},
        lambda idx, s_sub: float(average_precision_score(tb_y[strict_vs_mask][idx], s_sub)) if len(np.unique(tb_y[strict_vs_mask][idx])) == 2 else None,
        tb_comparisons, n_draws=2000, seed=20261120,
    )

    # Biophysical diagnostic summary on Track A test (3D steric clash vs empirical delta_pActivity)
    phys_rows = {r["compound_id"]: r for r in csv.DictReader((root / "data/processed/features/pocket_physics_3d_table.csv").open())}
    dE_vec = np.array([float(phys_rows[cid]["delta_vina_mta_minus_sam"]) for cid in ta_cids])
    steric_vec = np.array([float(phys_rows[cid]["steric_switch_differential"]) for cid in ta_cids])
    o3a_diff_vec = np.array([float(phys_rows[cid]["delta_o3a_mta_minus_sam"]) for cid in ta_cids])

    biophys_diag = {
        "track_a_test_n": len(ta_cids),
        "mean_delta_vina_mta_minus_sam_kcal": {
            "mta_cooperative_y3_eq_1": float(dE_vec[ta_y3 == 1].mean()),
            "non_cooperative_y3_eq_0": float(dE_vec[ta_y3 == 0].mean()),
        },
        "mean_steric_switch_differential": {
            "mta_cooperative_y3_eq_1": float(steric_vec[ta_y3 == 1].mean()),
            "non_cooperative_y3_eq_0": float(steric_vec[ta_y3 == 0].mean()),
        },
        "mean_delta_o3a_mta_minus_sam": {
            "mta_cooperative_y3_eq_1": float(o3a_diff_vec[ta_y3 == 1].mean()),
            "non_cooperative_y3_eq_0": float(o3a_diff_vec[ta_y3 == 0].mean()),
        },
        "spearman_steric_switch_vs_delta_p": float(spearmanr(steric_vec, ta_dp).statistic),
        "spearman_neg_delta_vina_vs_delta_p": float(spearmanr(-dE_vec, ta_dp).statistic),
    }

    report = {
        "status": "SCORED_LOCKED_PREDICTIONS_V2",
        "protocol_id": "PRMT5-V2-3D-COOPERATIVITY-2026",
        "protocol_freeze_sha256": sha256_file(root / "receipts/protocol_freeze.json"),
        "prediction_seal_sha256": sha256_file(root / "receipts/prediction_seal.json"),
        "model_selection_lock_sha256": sha256_file(root / "manifests/model_selection_lock.json"),
        "sealed_predictions_sha256": lock["sealed_predictions_sha256"],
        "track_a": {
            "n": int(len(ta_y3)),
            "scaffolds": int(len(set(ta_scaff))),
            "sources": int(len(set(ta_fams))),
            "positives_y_coop_3x": int(ta_y3.sum()),
            "positives_y_coop_10x": int(ta_y10.sum()),
            "metrics": ta_metrics,
            "bootstrap_2000": {
                "ap_y_coop_3x": ta_boot_ap3,
                "ap_y_coop_10x": ta_boot_ap10,
                "roc_auc_y_coop_3x": ta_boot_auc3,
                "pooled_c_index_delta_p": ta_boot_cidx,
            },
        },
        "track_b": {
            "n": int(len(tb_y)),
            "scaffolds": int(len(set(tb_scaff))),
            "assays": int(len(set(tb_assays))),
            "positives_le100nM": int(tb_y.sum()),
            "weak_gt1000nM": int(tb_weak.sum()),
            "chiral_compounds": int(tb_chiral.sum()),
            "metrics": tb_metrics,
            "bootstrap_2000": {
                "pooled_ap": tb_boot_ap,
                "pooled_roc_auc": tb_boot_auc,
                "within_assay_weighted_c_index": tb_boot_cidx,
                "strict_active_vs_weak_ap": tb_boot_strict_ap,
            },
        },
        "biophysical_diagnostics": biophys_diag,
    }
    boot_summary = {
        "seed": 20261120,
        "draws": 2000,
        "track_a_ap3": ta_boot_ap3,
        "track_b_ap": tb_boot_ap,
        "track_b_within_c_index": tb_boot_cidx,
    }
    return report, boot_summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--output-dir", type=Path, help="Optional empty directory for deterministic replay verification")
    args = parser.parse_args()
    root = args.root.resolve()

    if args.output_dir is None:
        metrics_path = root / "reports/results/final_metrics_v2.json"
        boot_path = root / "reports/results/bootstrap_summary_v2.json"
    else:
        metrics_path = args.output_dir.resolve() / "final_metrics_v2.json"
        boot_path = args.output_dir.resolve() / "bootstrap_summary_v2.json"

    require(
        not metrics_path.exists() and not boot_path.exists(),
        f"Scoring output already exists ({metrics_path}); refusing to overwrite. Use --output-dir for replay.",
    )

    report, boot_summary = evaluate_all(root)
    for path, payload in ((metrics_path, report), (boot_path, boot_summary)):
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("x", encoding="utf-8") as fh:
            fh.write(json.dumps(payload, indent=2, allow_nan=False) + "\n")

    print("=== STAGE 3 LOCKED TEST EVALUATION COMPLETE ===")
    print("Track A (y_coop_3x AP):", {m: round(report["track_a"]["metrics"][m]["primary_y_coop_3x"]["average_precision"], 4) for m in report["track_a"]["metrics"]})
    print("Track A (y_coop_10x AP):", {m: round(report["track_a"]["metrics"][m]["stringent_y_coop_10x"]["average_precision"], 4) for m in report["track_a"]["metrics"]})
    print("Track A (delta_p C-index):", {m: report["track_a"]["metrics"][m]["continuous_delta_p"]["pooled_c_index"] for m in report["track_a"]["metrics"]})
    print("Track B (Pooled AP):", {m: round(report["track_b"]["metrics"][m]["pooled_ranking"]["average_precision"], 4) for m in report["track_b"]["metrics"]})
    print("Track B (Within-Assay C-index):", {m: report["track_b"]["metrics"][m]["within_assay_concordance"]["weighted_c_index"] for m in report["track_b"]["metrics"]})
    print("Track B (Sensitivity at Val Threshold):", {m: report["track_b"]["metrics"][m]["threshold_calibration"]["sensitivity"] for m in report["track_b"]["metrics"]})


if __name__ == "__main__":
    main()
