"""Automated zero-leakage, SHA-256 seal-chain, replay, and report-consistency contract tests for PRMT5 v2."""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from chemistry_v2 import standardize_structure
from evaluate_locked import evaluate_all, verify_full_seal_chain


def sha256_str(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_csv(rel: str) -> list[dict[str, str]]:
    with (ROOT / rel).open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


class TestZeroLeakageAndQuarantine:
    @pytest.mark.parametrize("track_rel", [
        "data/processed/track_a_structures.csv",
        "data/processed/track_b_structures.csv",
    ])
    def test_zero_cross_split_leakage_on_all_four_axes(self, track_rel: str) -> None:
        rows = load_csv(track_rel)
        by_split: dict[str, list[dict[str, str]]] = {"train": [], "validation": [], "test": []}
        for r in rows:
            by_split[r["split"]].append(r)
        assert all(len(v) > 0 for v in by_split.values())

        for s1, s2 in [("train", "validation"), ("train", "test"), ("validation", "test")]:
            for col in ("canonical_smiles", "nonisomeric_smiles", "scaffold_group", "source_family"):
                v1 = {r[col] for r in by_split[s1]}
                v2 = {r[col] for r in by_split[s2]}
                assert len(v1 & v2) == 0, f"Leakage in {track_rel} ({s1} vs {s2}) on {col}: {v1 & v2}"

    @pytest.mark.parametrize("track_rel", [
        "data/processed/track_a_structures.csv",
        "data/processed/track_b_structures.csv",
    ])
    def test_zero_overlap_with_sealed_v1_sk_246_test_set(self, track_rel: str) -> None:
        sk = json.loads((ROOT / "manifests/v1_sk_quarantine.json").read_text(encoding="utf-8"))
        assert sk["count"] == 246
        sk_iso = set(sk["iso_sha256"])
        sk_non = set(sk["noniso_sha256"])
        sk_sc = set(sk["scaffold_sha256"])

        for r in load_csv(track_rel):
            assert sha256_str(r["canonical_smiles"]) not in sk_iso
            assert sha256_str(r["nonisomeric_smiles"]) not in sk_non
            assert sha256_str(r["scaffold_group"]) not in sk_sc
            assert r["source_family"] != sk["source_family"]

    def test_test_labels_physically_excluded_from_structures_and_dev_labels(self) -> None:
        ta_s = load_csv("data/processed/track_a_structures.csv")
        tb_s = load_csv("data/processed/track_b_structures.csv")
        assert "y_coop_3x" not in ta_s[0] and "delta_p_activity" not in ta_s[0]
        assert "binary_y_le100nM" not in tb_s[0] and "p_activity" not in tb_s[0]

        ta_dev = load_csv("data/processed/track_a_dev_labels.csv")
        tb_dev = load_csv("data/processed/track_b_dev_labels.csv")
        assert {r["split"] for r in ta_dev} == {"train", "validation"}
        assert {r["split"] for r in tb_dev} == {"train", "validation"}

        ta_te = load_csv("data/processed/sealed_test_labels/track_a_test_labels.csv")
        tb_te = load_csv("data/processed/sealed_test_labels/track_b_test_labels.csv")
        assert {r["split"] for r in ta_te} == {"test"}
        assert {r["split"] for r in tb_te} == {"test"}

    def test_post2023_crystals_strictly_quarantined_from_vina_receptors(self) -> None:
        pdb_m = json.loads((ROOT / "manifests/pdb_receptors.json").read_text(encoding="utf-8"))
        for rec in pdb_m["receptors"]:
            year = int(rec["release_date"][:4])
            if year >= 2023:
                assert rec["quarantined_post2023_diagnostic_only"] is True
                assert "receptor_pdbqt" not in rec
            else:
                assert rec["quarantined_post2023_diagnostic_only"] is False

    def test_chemistry_standardization_idempotence(self) -> None:
        sample = load_csv("data/processed/track_a_structures.csv")[:25] + load_csv("data/processed/track_b_structures.csv")[:25]
        for r in sample:
            st = standardize_structure(r["canonical_smiles"])
            assert st["canonical_smiles"] == r["canonical_smiles"]
            assert st["nonisomeric_smiles"] == r["nonisomeric_smiles"]
            assert st["scaffold_group"] == r["scaffold_group"]


class TestSealChainAndTamperRejection:
    def test_live_seal_chain_passes(self) -> None:
        freeze, seal, lock = verify_full_seal_chain(ROOT)
        assert freeze["models_fitted"] is False
        assert seal["test_labels_accessed"] is False
        assert lock["test_metrics_computed"] is False

    def _copy_minimal_tree(self, tmp_path: Path) -> Path:
        sandbox = tmp_path / "repo"
        for rel_dir in ("configs", "manifests", "receipts", "scripts", "data/processed", "data/raw_curated", "data/predictions"):
            shutil.copytree(ROOT / rel_dir, sandbox / rel_dir)
        return sandbox

    def test_altered_prediction_score_rejected(self, tmp_path: Path) -> None:
        sb = self._copy_minimal_tree(tmp_path)
        p = sb / "data/predictions/sealed_predictions.csv"
        txt = p.read_text(encoding="utf-8")
        lines = txt.splitlines()
        parts = lines[1].split(",")
        parts[-1] = "0.99999999" if parts[-1] != "0.99999999" else "0.11111111"
        lines[1] = ",".join(parts)
        p.write_text("\n".join(lines) + "\n", encoding="utf-8")
        with pytest.raises(ValueError, match="Tampered sealed_predictions.csv"):
            verify_full_seal_chain(sb)

    def test_altered_protocol_yaml_rejected(self, tmp_path: Path) -> None:
        sb = self._copy_minimal_tree(tmp_path)
        p = sb / "configs/protocol_v2.yaml"
        p.write_text(p.read_text(encoding="utf-8") + "\n# unauthorized edit\n", encoding="utf-8")
        with pytest.raises(ValueError, match="Tampered frozen artifact configs/protocol_v2.yaml"):
            verify_full_seal_chain(sb)

    def test_prediction_seal_leak_flag_violation_rejected(self, tmp_path: Path) -> None:
        sb = self._copy_minimal_tree(tmp_path)
        sp = sb / "receipts/prediction_seal.json"
        d = json.loads(sp.read_text(encoding="utf-8"))
        d["test_labels_accessed"] = True
        sp.write_text(json.dumps(d, indent=2) + "\n", encoding="utf-8")
        with pytest.raises(ValueError, match="test_labels_accessed=False"):
            verify_full_seal_chain(sb)

    def test_deterministic_replay_matches_published_json_and_refuses_overwrite(self, tmp_path: Path) -> None:
        replay_dir = tmp_path / "replay_out"
        cmd = [sys.executable, str(ROOT / "scripts/evaluate_locked.py"), "--root", str(ROOT), "--output-dir", str(replay_dir)]
        res = subprocess.run(cmd, capture_output=True, text=True, check=True)
        assert "STAGE 3 LOCKED TEST EVALUATION COMPLETE" in res.stdout

        pub_metrics = json.loads((ROOT / "reports/results/final_metrics_v2.json").read_text(encoding="utf-8"))
        rep_metrics = json.loads((replay_dir / "final_metrics_v2.json").read_text(encoding="utf-8"))
        assert pub_metrics == rep_metrics

        # Second invocation to same non-empty directory must fail rather than overwrite!
        res_dup = subprocess.run(cmd, capture_output=True, text=True)
        assert res_dup.returncode != 0
        assert "refusing to overwrite" in res_dup.stderr


class TestMarkdownReportExactConsistency:
    def test_saved_metrics_and_cohort_counts_appear_accurately_in_reports(self) -> None:
        m0 = json.loads((ROOT / "manifests/m0_source_audit.json").read_text(encoding="utf-8"))
        fm = json.loads((ROOT / "reports/results/final_metrics_v2.json").read_text(encoding="utf-8"))

        ta_m2_ap3 = f"{fm['track_a']['metrics']['M2']['primary_y_coop_3x']['average_precision']:.4f}"
        ta_b1_ap3 = f"{fm['track_a']['metrics']['B1']['primary_y_coop_3x']['average_precision']:.4f}"
        ta_p0_ap3 = f"{fm['track_a']['metrics']['P0']['primary_y_coop_3x']['average_precision']:.4f}"
        tb_m2_ap = f"{fm['track_b']['metrics']['M2']['pooled_ranking']['average_precision']:.4f}"
        tb_b1_ap = f"{fm['track_b']['metrics']['B1']['pooled_ranking']['average_precision']:.4f}"
        tb_m2_ci = f"{fm['track_b']['metrics']['M2']['within_assay_concordance']['weighted_c_index']:.4f}"

        for doc_rel in ("README.md", "README.zh-CN.md", "reports/RESULTS.md", "reports/RESULTS.zh-CN.md"):
            doc_path = ROOT / doc_rel
            assert doc_path.exists(), f"Missing required report: {doc_rel}"
            txt = doc_path.read_text(encoding="utf-8")
            assert str(m0["track_a"]["splits"]["test"]["rows"]) in txt
            assert str(m0["track_b"]["splits"]["test"]["rows"]) in txt
            for metric_str in (ta_m2_ap3, ta_b1_ap3, ta_p0_ap3, tb_m2_ap, tb_b1_ap, tb_m2_ci):
                assert metric_str in txt, f"Exact metric {metric_str} missing from {doc_rel}"
