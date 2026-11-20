#!/usr/bin/env python3
"""Stage 3B: Generate publication-grade biophysical & benchmark figures for PRMT5 v2."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


PALETTE = {
    "State1": "#1f77b4",
    "State2": "#d62728",
    "State3": "#2ca02c",
    "B0": "#8c8c8c",
    "P0": "#ff7f0e",
    "B1": "#7f7f7f",
    "M1": "#1f77b4",
    "M2": "#d62728",
}


def setup_style() -> None:
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 10,
        "axes.titlesize": 11,
        "axes.titleweight": "bold",
        "axes.labelsize": 10,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "figure.dpi": 220,
    })


def plot_fig1_steric_switch(root: Path, out_dir: Path, report: dict) -> None:
    pdb_manifest = json.loads((root / "manifests/pdb_receptors.json").read_text(encoding="utf-8"))
    crystals = pdb_manifest["crystallographic_steric_switch_audit"]

    ta_struct = [r for r in csv.DictReader((root / "data/processed/track_a_structures.csv").open()) if r["split"] == "test"]
    ta_labels = {r["row_id"]: r for r in csv.DictReader((root / "data/processed/sealed_test_labels/track_a_test_labels.csv").open())}
    phys_map = {r["compound_id"]: r for r in csv.DictReader((root / "data/processed/features/pocket_physics_3d_table.csv").open())}
    preds = [r for r in csv.DictReader((root / "data/predictions/sealed_predictions.csv").open()) if r["track"] == "track_a" and r["model"] == "P0"]
    p0_map = {r["row_id"]: float(r["score"]) for r in preds}

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13.2, 5.4))

    # Panel A: 11 Experimental Co-Crystal Structures in (min_d_SAM_pocket, min_d_Phe327_MTA) space
    ax1.axvspan(0.0, 2.8, color="#fde0dd", alpha=0.55, label="Steric clash with SAM:CE / Glu435 (<2.8 Å)")
    ax1.axhspan(0.0, 2.5, color="#deebf7", alpha=0.55, label="Steric clash with MTA-state Phe327 (<2.5 Å)")
    ax1.axvline(2.8, color="#c51b8a", linestyle="--", linewidth=1.0, alpha=0.8)
    ax1.axhline(2.5, color="#2b8cbe", linestyle="--", linewidth=1.0, alpha=0.8)

    state_labels_seen = set()
    label_offsets = {
        "4X61": (0.06, 0.12),
        "5EML": (0.06, -0.18),
        "6CKC": (0.06, 0.08),
        "6RLL": (0.06, -0.16),
        "7MX7": (0.06, 0.08),
        "7S1S": (0.06, -0.15),
        "8CSG": (-0.42, -0.22),
        "8VEY": (0.07, 0.04),
        "9C10": (-0.52, 0.12),
        "9N3R": (0.07, 0.08),
        "9T43": (-0.18, -0.24),
    }
    for c in crystals:
        x = min(c["min_dist_SAM_CE_A"], c["min_dist_Glu435_SAM_state_A"])
        y = c["min_dist_Phe327_MTA_state_A"]
        st = c["state"].split("_")[0]
        lbl_map = {
            "State1": "State 1: PRMT5•SAM + THIQ (pre-2023)",
            "State2": "State 2: PRMT5•MTA cooperative (★=post-2023)",
            "State3": "State 3: SAM-pocket nucleoside (pre-2023)",
        }
        lbl = lbl_map[st] if st not in state_labels_seen else None
        state_labels_seen.add(st)
        marker = "*" if c["quarantined_post2023"] else "o"
        size = 150 if c["quarantined_post2023"] else 85
        ax1.scatter(x, y, c=PALETTE[st], marker=marker, s=size, edgecolors="black", linewidths=0.7, zorder=4, label=lbl)
        dx, dy = label_offsets.get(c["pdb_id"], (0.06, 0.07))
        ax1.annotate(f"{c['pdb_id']} ({c['ligand_name']})", (x + dx, y + dy), fontsize=7.8, zorder=5)

    ax1.set_xlim(0.2, 4.35)
    ax1.set_ylim(0.15, 4.35)
    ax1.set_xlabel("Min Distance to SAM:CE or SAM-State Glu435 (Å)")
    ax1.set_ylabel("Min Distance to MTA-State Swung-In Phe327 (Å)")
    ax1.set_title("A. Experimental Co-Crystal 3D Steric-Switch Orthants (11 PDBs)")
    ax1.legend(loc="upper right", fontsize=7.4, framealpha=0.94)

    # Panel B: Track A Held-Out 2024-2025 Paired Test Cohort (n=98)
    fam_colors = {
        "MECHSWITCH_JMC2025": "#d62728",
        "TNG908_JMC2024": "#1f77b4",
        "AZ_JMC2024": "#2ca02c",
        "AM9747_JMC2025": "#9467bd",
    }
    for fam, col in fam_colors.items():
        xs, ys = [], []
        for r in ta_struct:
            if r["source_family"] == fam:
                xs.append(p0_map[r["row_id"]])
                ys.append(float(ta_labels[r["row_id"]]["delta_p_activity"]))
        ax2.scatter(xs, ys, c=col, s=48, alpha=0.82, edgecolors="white", linewidths=0.5, label=f"{fam} (n={len(xs)})")

    ax2.axhline(0.5, color="#444444", linestyle="--", linewidth=1.0, label="3.16x MTA-cooperative boundary (Δp=0.5)")
    ax2.axhline(0.0, color="#999999", linestyle=":", linewidth=0.9)
    rho_p0 = report["track_a"]["metrics"]["P0"]["continuous_delta_p"]["spearman_rho_pooled"]
    cidx_p0 = report["track_a"]["metrics"]["P0"]["continuous_delta_p"]["pooled_c_index"]
    ax2.set_xlabel("Zero-Shot 3D Cofactor-Pocket Physics Score P0 (No Labels Fitted)")
    ax2.set_ylabel("Experimental ΔpActivity = pActivity(MTA) − pActivity(SAM)")
    ax2.set_title(f"B. Zero-Shot 3D Physics on Untouched 2024–2025 Series (ρ={rho_p0:+.2f}, C={cidx_p0:.2f})")
    ax2.legend(loc="lower right", fontsize=8.0, framealpha=0.92)

    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(out_dir / f"fig1_cofactor_steric_switch.{ext}", bbox_inches="tight")
    plt.close(fig)


def plot_fig2_benchmark_ci(out_dir: Path, report: dict) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(15.0, 5.0))

    # Panel A: Track A Held-Out AP (y_coop_3x and y_coop_10x) with 95% CIs
    ax = axes[0]
    ta_models = ["B0", "P0", "B1", "M1", "M2"]
    x = np.arange(len(ta_models))
    w = 0.36

    b3 = report["track_a"]["bootstrap_2000"]["ap_y_coop_3x"]["models"]
    b10 = report["track_a"]["bootstrap_2000"]["ap_y_coop_10x"]["models"]

    pts3 = [b3[m]["point"] for m in ta_models]
    err3_lo = [b3[m]["point"] - b3[m]["ci95"][0] for m in ta_models]
    err3_hi = [b3[m]["ci95"][1] - b3[m]["point"] for m in ta_models]

    pts10 = [b10[m]["point"] for m in ta_models]
    err10_lo = [b10[m]["point"] - b10[m]["ci95"][0] for m in ta_models]
    err10_hi = [b10[m]["ci95"][1] - b10[m]["point"] for m in ta_models]

    ax.bar(x - w / 2, pts3, width=w, color="#2b8cbe", alpha=0.88, yerr=[err3_lo, err3_hi], capsize=3.5, label="Primary ≥3.16x MTA-Coop AP")
    ax.bar(x + w / 2, pts10, width=w, color="#e34a33", alpha=0.88, yerr=[err10_lo, err10_hi], capsize=3.5, label="Stringent ≥10x MTA-Selectivity AP")
    ax.set_xticks(x, ta_models)
    ax.set_ylim(0.0, 1.05)
    ax.set_ylabel("Held-Out Average Precision (95% Scaffold CI)")
    ax.set_title("A. Track A: Cofactor Cooperativity (n=98, 2024–2025)")
    ax.legend(loc="upper left", fontsize=8.0)

    # Panel B: Track B Held-Out Within-Assay C-index & Pooled AP
    ax = axes[1]
    tb_models = ["B0", "B1", "M1", "M2"]
    xb = np.arange(len(tb_models))
    b_ap = report["track_b"]["bootstrap_2000"]["pooled_ap"]["models"]
    b_ci = report["track_b"]["bootstrap_2000"]["within_assay_weighted_c_index"]["models"]

    pts_ap = [b_ap[m]["point"] for m in tb_models]
    ap_lo = [b_ap[m]["point"] - b_ap[m]["ci95"][0] for m in tb_models]
    ap_hi = [b_ap[m]["ci95"][1] - b_ap[m]["point"] for m in tb_models]

    pts_ci = [0.5 if m == "B0" else b_ci[m]["point"] for m in tb_models]
    ci_lo = [0.0 if m == "B0" else b_ci[m]["point"] - b_ci[m]["ci95"][0] for m in tb_models]
    ci_hi = [0.0 if m == "B0" else b_ci[m]["ci95"][1] - b_ci[m]["point"] for m in tb_models]

    ax.bar(xb - w / 2, pts_ci, width=w, color="#31a354", alpha=0.88, yerr=[ci_lo, ci_hi], capsize=3.5, label="Within-Assay Weighted C-index")
    ax.bar(xb + w / 2, pts_ap, width=w, color="#756bb1", alpha=0.88, yerr=[ap_lo, ap_hi], capsize=3.5, label="Pooled ≤100 nM AP (prev=0.313)")
    ax.axhline(0.5, color="#888888", linestyle=":", linewidth=0.9)
    ax.set_xticks(xb, tb_models)
    ax.set_ylim(0.0, 1.05)
    ax.set_ylabel("Score (95% Paired Scaffold Bootstrap CI)")
    ax.set_title("B. Track B: Untouched 2024–2025 Potency (n=182)")
    ax.legend(loc="upper left", fontsize=8.0)

    # Panel C: Paired Scaffold Bootstrap Differences vs v1-Style B1
    ax = axes[2]
    diff_items = [
        ("Track A AP(3x): M2 − B1", report["track_a"]["bootstrap_2000"]["ap_y_coop_3x"]["paired_differences"]["M2_minus_B1"]),
        ("Track A AP(3x): M2 − M1", report["track_a"]["bootstrap_2000"]["ap_y_coop_3x"]["paired_differences"]["M2_minus_M1"]),
        ("Track A AP(3x): P0 − B1", report["track_a"]["bootstrap_2000"]["ap_y_coop_3x"]["paired_differences"]["P0_minus_B1"]),
        ("Track B AP: M2 − B1", report["track_b"]["bootstrap_2000"]["pooled_ap"]["paired_differences"]["M2_minus_B1"]),
        ("Track B AP: M1 − B1", report["track_b"]["bootstrap_2000"]["pooled_ap"]["paired_differences"]["M1_minus_B1"]),
        ("Track B C-idx: M2 − B1", report["track_b"]["bootstrap_2000"]["within_assay_weighted_c_index"]["paired_differences"]["M2_minus_B1"]),
    ]
    yd = np.arange(len(diff_items))[::-1]
    for yi, (lbl, d) in zip(yd, diff_items):
        pt = d["point_diff"]
        lo, hi = d["ci95"]
        col = "#d62728" if lo > 0 else ("#1f77b4" if pt > 0 else "#666666")
        ax.errorbar(pt, yi, xerr=[[pt - lo], [hi - pt]], fmt="o", color=col, capsize=4, markersize=6, linewidth=1.6)
    ax.axvline(0.0, color="black", linestyle="--", linewidth=1.0)
    ax.set_yticks(yd, [x[0] for x in diff_items], fontsize=8.5)
    ax.set_xlabel("Paired Metric Gain (2,000 Scaffold Resamples, 95% CI)")
    ax.set_title("C. Paired Scaffold-Bootstrap Contrasts")

    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(out_dir / f"fig2_track_a_b_benchmark_ci.{ext}", bbox_inches="tight")
    plt.close(fig)


def plot_fig3_calibration_and_enrichment(root: Path, out_dir: Path, report: dict) -> None:
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.8, 5.0))

    # Panel A: v1 SK Threshold Collapse Context vs v2 Track B Sensitivity & Balanced Accuracy
    categories = [
        "v1 SK (MoLFormer)\n[Historical v1 Ref]",
        "v2 Track B: B1\n(v1-Style Achiral 2D)",
        "v2 Track B: M1\n(Chiral + Assay Ranker)",
        "v2 Track B: M2\n(Chiral + 3D + Ranker)",
    ]
    tb_m = report["track_b"]["metrics"]
    sens_vals = [
        2.0 / 69.0,  # v1 published 2/69 positives recovered on SK at Pharmablock threshold!
        tb_m["B1"]["threshold_calibration"]["sensitivity"],
        tb_m["M1"]["threshold_calibration"]["sensitivity"],
        tb_m["M2"]["threshold_calibration"]["sensitivity"],
    ]
    bal_vals = [
        0.5 * (2.0 / 69.0 + 175.0 / 177.0),  # v1 SK balanced accuracy
        tb_m["B1"]["threshold_calibration"]["balanced_accuracy"],
        tb_m["M1"]["threshold_calibration"]["balanced_accuracy"],
        tb_m["M2"]["threshold_calibration"]["balanced_accuracy"],
    ]
    x = np.arange(len(categories))
    w = 0.35
    ax1.bar(x - w / 2, sens_vals, width=w, color="#fb6a4a", label="Sensitivity (Recall of ≤100 nM Actives)")
    ax1.bar(x + w / 2, bal_vals, width=w, color="#2b8cbe", label="Balanced Accuracy at Locked Threshold")
    for xi, (sv, bv) in enumerate(zip(sens_vals, bal_vals)):
        ax1.text(xi - w / 2, sv + 0.02, f"{sv:.2f}", ha="center", fontsize=8.2)
        ax1.text(xi + w / 2, bv + 0.02, f"{bv:.2f}", ha="center", fontsize=8.2)
    ax1.set_xticks(x, categories, fontsize=8.5)
    ax1.set_ylim(0.0, 1.05)
    ax1.set_ylabel("Operating-Point Performance")
    ax1.set_title("A. Resolving Cross-Assay [SAM] Threshold Collapse")
    ax1.legend(loc="upper left", fontsize=8.2)

    # Panel B: Cumulative Screening Enrichment Curve on Track B Test (k = 1..60 out of 182)
    tb_labels = {r["row_id"]: int(r["binary_y_le100nM"]) for r in csv.DictReader((root / "data/processed/sealed_test_labels/track_b_test_labels.csv").open())}
    preds = [r for r in csv.DictReader((root / "data/predictions/sealed_predictions.csv").open()) if r["track"] == "track_b"]
    ks = np.arange(1, 61)
    n_tot = len(tb_labels)
    pos_tot = sum(tb_labels.values())
    ax2.plot(ks, ks * (pos_tot / n_tot), color="#888888", linestyle="--", linewidth=1.4, label="Random Selection Expectation")

    for m_name, col, lw in [("B1", "#636363", 1.6), ("M1", "#1f77b4", 2.0), ("M2", "#d62728", 2.3)]:
        m_rows = sorted(
            [r for r in preds if r["model"] == m_name],
            key=lambda r: (-float(r["score"]), r["compound_id"]),
        )
        cum_hits = np.cumsum([tb_labels[r["row_id"]] for r in m_rows])[:60]
        ax2.plot(ks, cum_hits, color=col, linewidth=lw, label=f"{m_name} (Top-25 hits = {int(cum_hits[24])})")

    ax2.axvline(25, color="#444444", linestyle=":", linewidth=1.0)
    ax2.set_xlabel("Top-K Screened Compounds (out of 182 Held-Out 2024–2025 Parents)")
    ax2.set_ylabel("Recovered Active Inhibitors (IC50 / Ki ≤ 100 nM)")
    ax2.set_title("B. Track B Virtual Screening Hit Recovery Curve")
    ax2.legend(loc="upper left", fontsize=8.5)

    fig.tight_layout()
    for ext in ("png", "svg"):
        fig.savefig(out_dir / f"fig3_v1_threshold_collapse_vs_v2_calibration.{ext}", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    root = args.root.resolve()

    setup_style()
    report = json.loads((root / "reports/results/final_metrics_v2.json").read_text(encoding="utf-8"))
    out_dir = root / "reports/figures"
    out_dir.mkdir(parents=True, exist_ok=True)

    plot_fig1_steric_switch(root, out_dir, report)
    plot_fig2_benchmark_ci(out_dir, report)
    plot_fig3_calibration_and_enrichment(root, out_dir, report)
    print("Generated publication figures in", out_dir)


if __name__ == "__main__":
    main()
