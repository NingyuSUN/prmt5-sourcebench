#!/usr/bin/env python3
"""Stage 1 (M0): Build zero-leakage Track A & Track B cohorts, verify invariants, and seal protocol_freeze.json.

No model fitting and no molecular/3D feature computation are performed in this script.
Test target labels are isolated into `data/processed/sealed_test_labels/` so Stage 2
training and featurization scripts never read test labels.
"""
from __future__ import annotations

import argparse
import collections
import csv
import datetime
import hashlib
import json
import math
from pathlib import Path
import sys
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
from chemistry_v2 import standardize_structure

KNOWN_REFERENCE_CONTROLS = {
    "O=C(NC[C@H](O)CN1CCc2ccccc2C1)c1cc(NC2COC2)ncn1": "EPZ015666",
    "CC(=O)N1CCC(Nc2cc(C(=O)NC[C@H](O)CN3CCc4ccccc4C3)ncn2)CC1": "GSK3326595",
    "Nc1nc2cc(CC[C@H]3C[C@@H](n4ccc5c(N)ncnc54)[C@H](O)[C@@H]3O)ccc2cc1Br": "JNJ-64619178",
    "Nc1ncnc2c1ccn2[C@@H]1O[C@H]([C@H](O)c2ccccc2)[C@@H](O)[C@H]1O": "LLY-283",
    "Cn1ncc(-c2ccc3c(=O)[nH]nc(CN)c3c2)c1-c1c(F)c(Cl)cc(OC2CC2)c1C#N": "MRTX1719",
    "Cn1ncc(-c2ccc3c(=O)[nH]nc(CN)c3c2)c1-c1sc2ccccc2c1C#N": "Mirati_ref_analog",
    "Cc1cc(NC(=O)C(=O)N2C[C@@H](C)CC[C@@H]2c2ccc3scnc3c2)cnc1N": "TNG908",
    "Cc1cc2cc(C(=O)N(Cc3ccc(C(F)(F)F)cn3)CC(C)C)ccc2nc1N": "AM-9747",
}

# Audited assay protocol covariates (SAM concentration in uM, readout family, displacement flag)
ASSAY_COVARIATES = {
    "EPIZYME_FLASHPLATE": {"sam_uM": 0.075, "readout": "FlashPlate_radiometric", "is_displacement": 0},
    "CHEMBL5733161": {"sam_uM": 1.0, "readout": "homogeneous_biochemical", "is_displacement": 0},
    "CHEMBL5734993": {"sam_uM": 1.0, "readout": "Tris_biochemical_nucleoside", "is_displacement": 0},
    "CHEMBL5737088": {"sam_uM": 1.0, "readout": "Tris_biochemical_nucleoside", "is_displacement": 0},
    "CHEMBL5737758": {"sam_uM": 1.0, "readout": "Tris_biochemical_nucleoside", "is_displacement": 0},
    "CHEMBL5739093": {"sam_uM": 10.0, "readout": "AlphaLISA", "is_displacement": 0},
    "CHEMBL5167913": {"sam_uM": 0.5, "readout": "3H_SAM_displacement", "is_displacement": 1},
    "CHEMBL5738764": {"sam_uM": 1.0, "readout": "FlashPlate_radiometric", "is_displacement": 0},
    "CHEMBL5738200": {"sam_uM": 1.0, "readout": "TR_FRET", "is_displacement": 0},
    "CHEMBL5623742": {"sam_uM": 1.0, "readout": "enzymatic_biochemical", "is_displacement": 0},
    "CHEMBL5544428": {"sam_uM": 1.0, "readout": "AlphaLISA", "is_displacement": 0},
    "CHEMBL5509028": {"sam_uM": 1.0, "readout": "SPA_radiometric", "is_displacement": 0},
    "CHEMBL6088129": {"sam_uM": 1.0, "readout": "AlphaLISA", "is_displacement": 0},
    "CHEMBL6090539": {"sam_uM": 1.0, "readout": "3H_SAM_enzymatic", "is_displacement": 0},
    "CHEMBL6157149": {"sam_uM": 10.0, "readout": "FP_peptide_displacement_SAM", "is_displacement": 1},
}


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sha256_str(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"Refusing to write empty CSV: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def verify_zero_leakage(track_name: str, split_rows: dict[str, list[dict]], sk_quarantine: dict) -> dict:
    """Hard check for zero cross-split overlap across iso, noniso, Murcko scaffold, family, and v1 SK."""
    sk_iso = set(sk_quarantine["iso_sha256"])
    sk_non = set(sk_quarantine["noniso_sha256"])
    sk_sc = set(sk_quarantine["scaffold_sha256"])

    sets_by_split = {}
    for sp, rows in split_rows.items():
        iso = {r["canonical_smiles"] for r in rows}
        non = {r["nonisomeric_smiles"] for r in rows}
        sc = {r["scaffold_group"] for r in rows}
        fam = {r["source_family"] for r in rows}
        sets_by_split[sp] = {"iso": iso, "noniso": non, "scaffold": sc, "family": fam}

        # Verify 0 overlap with v1 SK test set!
        ov_sk_iso = sum(1 for x in iso if sha256_str(x) in sk_iso)
        ov_sk_non = sum(1 for x in non if sha256_str(x) in sk_non)
        ov_sk_sc = sum(1 for x in sc if sha256_str(x) in sk_sc)
        ov_sk_fam = int(sk_quarantine["source_family"] in fam)
        if any((ov_sk_iso, ov_sk_non, ov_sk_sc, ov_sk_fam)):
            raise RuntimeError(
                f"{track_name} split {sp} leaked v1 SK test set: "
                f"iso={ov_sk_iso}, noniso={ov_sk_non}, scaff={ov_sk_sc}, fam={ov_sk_fam}"
            )

    pair_report = {}
    for s1, s2 in [("train", "validation"), ("train", "test"), ("validation", "test")]:
        a, b = sets_by_split[s1], sets_by_split[s2]
        ov = {k: len(a[k] & b[k]) for k in ("iso", "noniso", "scaffold", "family")}
        if any(ov.values()):
            raise RuntimeError(f"{track_name} leakage between {s1} and {s2}: {ov}")
        pair_report[f"{s1}_vs_{s2}"] = ov

    return {
        "splits": {
            sp: {
                "rows": len(split_rows[sp]),
                "unique_isomeric_parents": len(sets_by_split[sp]["iso"]),
                "unique_nonisomeric_parents": len(sets_by_split[sp]["noniso"]),
                "unique_murcko_scaffolds": len(sets_by_split[sp]["scaffold"]),
                "unique_source_families": len(sets_by_split[sp]["family"]),
                "v1_sk_overlap": 0,
            }
            for sp in ("train", "validation", "test")
        },
        "pairwise_overlaps": pair_report,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    root = args.root.resolve()

    protocol_path = root / "configs/protocol_v2.yaml"
    protocol = yaml.safe_load(protocol_path.read_text(encoding="utf-8"))
    sk_quarantine = json.loads((root / "manifests/v1_sk_quarantine.json").read_text(encoding="utf-8"))

    chembl_rows = list(csv.DictReader((root / "data/raw_curated/chembl_selected_activities.csv").open()))
    mir_rows = list(csv.DictReader((root / "data/raw_curated/mirati_table9a_unique.csv").open()))
    v1_rows = list(csv.DictReader((root / "data/raw_curated/v1_pre2020_cohorts.csv").open()))

    exclusions: list[dict] = []
    std_cache: dict[str, dict[str, str]] = {}

    def get_std(raw_smi: str) -> dict[str, str]:
        if raw_smi not in std_cache:
            std_cache[raw_smi] = standardize_structure(raw_smi)
        return std_cache[raw_smi]

    # Helper: aggregate ChEMBL single-assay parent measurements
    def aggregate_chembl_assay(assay_id: str, quarantine_controls: bool) -> dict[str, dict]:
        grouped = collections.defaultdict(list)
        for r in chembl_rows:
            if r["assay_chembl_id"] != assay_id:
                continue
            if r["data_validity_comment"]:
                exclusions.append({"track": "both", "assay_id": assay_id, "smiles": r["raw_canonical_smiles"], "reason": "data_validity_comment"})
                continue
            try:
                val = float(r["standard_value"])
            except (TypeError, ValueError):
                val = -1.0
            if val <= 0 or r["standard_units"] != "nM":
                exclusions.append({"track": "both", "assay_id": assay_id, "smiles": r["raw_canonical_smiles"], "reason": "nonpositive_or_non_nM"})
                continue
            st = get_std(r["raw_canonical_smiles"])
            if quarantine_controls and st["canonical_smiles"] in KNOWN_REFERENCE_CONTROLS:
                exclusions.append({
                    "track": "novel_series",
                    "assay_id": assay_id,
                    "smiles": st["canonical_smiles"],
                    "reason": f"historical_reference_control:{KNOWN_REFERENCE_CONTROLS[st['canonical_smiles']]}",
                })
                continue
            grouped[st["canonical_smiles"]].append((r, st, val))

        out = {}
        for iso_smi, items in sorted(grouped.items()):
            rels = {x[0]["standard_relation"] for x in items}
            vals = [x[2] for x in items]
            st = items[0][1]
            r0 = items[0][0]
            if rels == {"="}:
                if max(vals) / min(vals) > 3.0:
                    exclusions.append({"track": "both", "assay_id": assay_id, "smiles": iso_smi, "reason": "replicate_spread_gt_3x"})
                    continue
                val_agg = math.exp(sum(math.log(v) for v in vals) / len(vals))
                rel_agg = "="
            elif all(v > 100.0 for v in vals) and rels <= {"=", ">", ">="}:
                val_agg = max(vals)
                rel_agg = ">"
            elif all(v <= 100.0 for v in vals) and rels <= {"=", "<", "<="}:
                val_agg = min(vals)
                rel_agg = "<"
            else:
                exclusions.append({"track": "both", "assay_id": assay_id, "smiles": iso_smi, "reason": "ambiguous_censored_replicates"})
                continue
            out[iso_smi] = {
                **st,
                "relation": rel_agg,
                "value_nM": val_agg,
                "document_id": r0["document_chembl_id"],
                "publication_year": int(r0["document_year"]),
            }
        return out

    # =========================================================================
    # BUILD TRACK A: Cofactor-Cooperativity (PRMT5•MTA vs PRMT5•SAM)
    # =========================================================================
    track_a = {"train": [], "validation": [], "test": []}
    ta_seen_iso: set[str] = set()

    # A-Train 1: Mirati Table 9a direct paired FlashPlate H4 (2021-03-18)
    mir_grouped = collections.defaultdict(lambda: collections.defaultdict(list))
    for mr in mir_rows:
        val = float(mr["value_nM"])
        if val <= 0:
            exclusions.append({"track": "track_a", "assay_id": "MIRATI_T9A", "smiles": mr["raw_canonical_smiles"], "reason": "zero_or_negative_source_IC50"})
            continue
        st = get_std(mr["raw_canonical_smiles"])
        mir_grouped[st["canonical_smiles"]][mr["mta_condition"]].append((mr["relation"], val, st))

    for iso_smi, conds in sorted(mir_grouped.items()):
        if "MTA_2uM" in conds and "MTA_absent" in conds:
            m_set = {(r, v) for r, v, _ in conds["MTA_2uM"]}
            s_set = {(r, v) for r, v, _ in conds["MTA_absent"]}
            if len(m_set) != 1 or len(s_set) != 1:
                exclusions.append({"track": "track_a", "assay_id": "MIRATI_T9A", "smiles": iso_smi, "reason": "within_table_duplicate_disagreement"})
                continue
            (rm, vm), (rs, vs) = next(iter(m_set)), next(iter(s_set))
            st = conds["MTA_2uM"][0][2]
            dp = math.log10(vs / vm)
            y3 = int(dp >= 0.5 and vm < 10000.0 and rm in ("=", "<", "<="))
            y10 = int(dp >= 1.0 and vm <= 1000.0 and rm in ("=", "<", "<="))
            track_a["train"].append({
                "row_id": f"TA_TR_{st['compound_id']}",
                **st,
                "split": "train",
                "source_family": "MIRATI_72659905",
                "document_id": "US20210078994A1",
                "publication_year": 2021,
                "assay_pair_type": "FlashPlate_H4_MTA2uM_vs_noMTA",
                "is_direct_paired": 1,
                "y_coop_3x": y3,
                "y_coop_10x": y10,
                "delta_p_activity": round(dp, 5),
                "mta_rel": rm,
                "mta_value_nM": round(vm, 4),
                "sam_rel": rs,
                "sam_value_nM": round(vs, 4),
            })
            ta_seen_iso.add(iso_smi)

    # A-Train 2: Active (<=100 nM) 1st-gen SAM-cooperative Epizyme & SAM-pocket Nucleoside controls (y_coop=0)
    for r in v1_rows:
        if r["v1_split"] == "train" and int(r["binary_y_le100nM"]) == 1:
            st = get_std(r["raw_canonical_smiles"])
            if st["canonical_smiles"] not in ta_seen_iso:
                track_a["train"].append({
                    "row_id": f"TA_TR_{st['compound_id']}",
                    **st,
                    "split": "train",
                    "source_family": r["source_family"].split(";")[0],
                    "document_id": f"PATFAM_{r['source_family'].split(';')[0]}",
                    "publication_year": int(r["source_date_bound"][:4]),
                    "assay_pair_type": "State1_SAM_uncompetitive_active_control",
                    "is_direct_paired": 0,
                    "y_coop_3x": 0,
                    "y_coop_10x": 0,
                    "delta_p_activity": -1.0,
                    "mta_rel": ">",
                    "mta_value_nM": 1000.0,
                    "sam_rel": "<=",
                    "sam_value_nM": 100.0,
                })
                ta_seen_iso.add(st["canonical_smiles"])

    for aid in ("CHEMBL5734993", "CHEMBL5737088", "CHEMBL5737758"):
        for iso_smi, rec in aggregate_chembl_assay(aid, quarantine_controls=False).items():
            if rec["relation"] in ("=", "<") and rec["value_nM"] <= 100.0 and iso_smi not in ta_seen_iso:
                track_a["train"].append({
                    "row_id": f"TA_TR_{rec['compound_id']}",
                    "compound_id": rec["compound_id"],
                    "canonical_smiles": rec["canonical_smiles"],
                    "nonisomeric_smiles": rec["nonisomeric_smiles"],
                    "scaffold_group": rec["scaffold_group"],
                    "has_defined_chirality": rec["has_defined_chirality"],
                    "split": "train",
                    "source_family": "NUCPAT_FAM",
                    "document_id": rec["document_id"],
                    "publication_year": min(rec["publication_year"], 2021),
                    "assay_pair_type": "State3_SAM_pocket_nucleoside_active_control",
                    "is_direct_paired": 0,
                    "y_coop_3x": 0,
                    "y_coop_10x": 0,
                    "delta_p_activity": -1.0,
                    "mta_rel": ">",
                    "mta_value_nM": 1000.0,
                    "sam_rel": rec["relation"],
                    "sam_value_nM": round(rec["value_nM"], 4),
                })
                ta_seen_iso.add(iso_smi)

    ta_train_noniso = {r["nonisomeric_smiles"] for r in track_a["train"]}
    ta_train_scaff = {r["scaffold_group"] for r in track_a["train"]}

    def add_paired_assay_to_track_a(split: str, prefix: str, a_mta: str, a_sam: str, source_family: str, pair_desc: str, disallow_noniso: set[str], disallow_scaff: set[str]) -> None:
        m_dict = aggregate_chembl_assay(a_mta, quarantine_controls=True)
        s_dict = aggregate_chembl_assay(a_sam, quarantine_controls=True)
        for iso_smi in sorted(set(m_dict) & set(s_dict)):
            mrec, srec = m_dict[iso_smi], s_dict[iso_smi]
            if mrec["nonisomeric_smiles"] in disallow_noniso or mrec["scaffold_group"] in disallow_scaff:
                exclusions.append({"track": "track_a", "assay_id": f"{a_mta}+{a_sam}", "smiles": iso_smi, "reason": "cross_split_scaffold_or_noniso_quarantine"})
                continue
            vm, vs = mrec["value_nM"], srec["value_nM"]
            dp = math.log10(vs / vm)
            y3 = int(dp >= 0.5 and vm < 10000.0 and mrec["relation"] in ("=", "<"))
            y10 = int(dp >= 1.0 and vm <= 1000.0 and mrec["relation"] in ("=", "<"))
            track_a[split].append({
                "row_id": f"{prefix}_{mrec['compound_id']}",
                "compound_id": mrec["compound_id"],
                "canonical_smiles": mrec["canonical_smiles"],
                "nonisomeric_smiles": mrec["nonisomeric_smiles"],
                "scaffold_group": mrec["scaffold_group"],
                "has_defined_chirality": mrec["has_defined_chirality"],
                "split": split,
                "source_family": source_family,
                "document_id": mrec["document_id"],
                "publication_year": mrec["publication_year"],
                "assay_pair_type": pair_desc,
                "is_direct_paired": 1,
                "y_coop_3x": y3,
                "y_coop_10x": y10,
                "delta_p_activity": round(dp, 5),
                "mta_rel": mrec["relation"],
                "mta_value_nM": round(vm, 4),
                "sam_rel": srec["relation"],
                "sam_value_nM": round(vs, 4),
            })

    # A-Validation (2019-2022)
    add_paired_assay_to_track_a(
        "validation", "TA_VAL", "CHEMBL5109533", "CHEMBL5109534",
        "TNG_BMC2022", "HCT116_MTAP_KO_vs_WT_SDMA", ta_train_noniso, ta_train_scaff
    )
    for iso_smi, rec in aggregate_chembl_assay("CHEMBL5130634", quarantine_controls=True).items():
        if rec["value_nM"] <= 200.0 and rec["nonisomeric_smiles"] not in ta_train_noniso and rec["scaffold_group"] not in ta_train_scaff:
            track_a["validation"].append({
                "row_id": f"TA_VAL_{rec['compound_id']}",
                "compound_id": rec["compound_id"],
                "canonical_smiles": rec["canonical_smiles"],
                "nonisomeric_smiles": rec["nonisomeric_smiles"],
                "scaffold_group": rec["scaffold_group"],
                "has_defined_chirality": rec["has_defined_chirality"],
                "split": "validation",
                "source_family": "CALU6_JMC2022",
                "document_id": rec["document_id"],
                "publication_year": rec["publication_year"],
                "assay_pair_type": "Calu6_MTAP_null_SDMA_selective",
                "is_direct_paired": 0,
                "y_coop_3x": 1,
                "y_coop_10x": 1,
                "delta_p_activity": 1.2,
                "mta_rel": rec["relation"],
                "mta_value_nM": round(rec["value_nM"], 4),
                "sam_rel": ">",
                "sam_value_nM": round(rec["value_nM"] * 15.85, 4),
            })
    for r in v1_rows:
        if r["v1_split"] == "validation" and int(r["binary_y_le100nM"]) == 1:
            st = get_std(r["raw_canonical_smiles"])
            if st["nonisomeric_smiles"] not in ta_train_noniso and st["scaffold_group"] not in ta_train_scaff:
                ev = json.loads(r["source_evidence"])[0]
                vals = [float(v) for v in ev["values"].values()]
                v_sam = math.exp(sum(math.log(x) for x in vals) / len(vals))
                track_a["validation"].append({
                    "row_id": f"TA_VAL_{st['compound_id']}",
                    **st,
                    "split": "validation",
                    "source_family": "PHARMABLOCK_2019",
                    "document_id": "WO2019173804A1",
                    "publication_year": 2019,
                    "assay_pair_type": "State1_SAM10uM_AlphaLISA_active_control",
                    "is_direct_paired": 0,
                    "y_coop_3x": 0,
                    "y_coop_10x": 0,
                    "delta_p_activity": -1.0,
                    "mta_rel": ">",
                    "mta_value_nM": 1000.0,
                    "sam_rel": "=",
                    "sam_value_nM": round(v_sam, 4),
                })
    val_seen_iso = {r["canonical_smiles"] for r in track_a["validation"]}
    for iso_smi, rec in aggregate_chembl_assay("CHEMBL5167913", quarantine_controls=True).items():
        if rec["value_nM"] <= 100.0 and iso_smi not in val_seen_iso and rec["nonisomeric_smiles"] not in ta_train_noniso and rec["scaffold_group"] not in ta_train_scaff:
            track_a["validation"].append({
                "row_id": f"TA_VAL_{rec['compound_id']}",
                "compound_id": rec["compound_id"],
                "canonical_smiles": rec["canonical_smiles"],
                "nonisomeric_smiles": rec["nonisomeric_smiles"],
                "scaffold_group": rec["scaffold_group"],
                "has_defined_chirality": rec["has_defined_chirality"],
                "split": "validation",
                "source_family": "JMC2022_5154928",
                "document_id": rec["document_id"],
                "publication_year": rec["publication_year"],
                "assay_pair_type": "3H_SAM_displacement_active_control",
                "is_direct_paired": 0,
                "y_coop_3x": 0,
                "y_coop_10x": 0,
                "delta_p_activity": -1.0,
                "mta_rel": ">",
                "mta_value_nM": 1000.0,
                "sam_rel": rec["relation"],
                "sam_value_nM": round(rec["value_nM"], 4),
            })

    ta_dev_noniso = ta_train_noniso | {r["nonisomeric_smiles"] for r in track_a["validation"]}
    ta_dev_scaff = ta_train_scaff | {r["scaffold_group"] for r in track_a["validation"]}

    # A-Test (2024-2025): 100% Direct Experimental Paired Biochemical Series!
    for a_mta, a_sam, fam, desc in [
        ("CHEMBL5547645", "CHEMBL5547646", "TNG908_JMC2024", "FA_peptide_displacement_MTA_vs_SAM"),
        ("CHEMBL5577992", "CHEMBL5577993", "AZ_JMC2024", "MTaseGlo_H4_SAM_plus_MTA_vs_minus_MTA"),
        ("CHEMBL6157148", "CHEMBL6157149", "MECHSWITCH_JMC2025", "FA_Ki_MTA_vs_SAM"),
        ("CHEMBL6119736", "CHEMBL6119737", "AM9747_JMC2025", "MTaseGlo_H2A_SAM_MTA_vs_SAM"),
    ]:
        add_paired_assay_to_track_a("test", "TA_TE", a_mta, a_sam, fam, desc, ta_dev_noniso, ta_dev_scaff)

    if not all(r["is_direct_paired"] == 1 for r in track_a["test"]):
        raise RuntimeError("Track A test must consist 100% of direct experimentally paired compounds")

    # =========================================================================
    # BUILD TRACK B: Cross-Document Potency Ranking on 2023-2025 PRMT5 Series
    # =========================================================================
    track_b = {"train": [], "validation": [], "test": []}
    tb_seen_train: set[str] = set()

    # B-Train 1: v1 Epizyme parents (75 nM SAM), deduplicated under v2 tautomer standardization
    epz_grouped = collections.defaultdict(list)
    for r in v1_rows:
        if r["v1_split"] == "train":
            st = get_std(r["raw_canonical_smiles"])
            epz_grouped[st["canonical_smiles"]].append((r, st))
    for iso_smi, items in sorted(epz_grouped.items()):
        ys = {int(x[0]["binary_y_le100nM"]) for x in items}
        if len(ys) != 1:
            exclusions.append({"track": "track_b", "assay_id": "EPIZYME_FLASHPLATE", "smiles": iso_smi, "reason": "tautomer_merge_label_conflict"})
            continue
        r, st = items[0]
        ev = json.loads(r["source_evidence"])[0]
        lo = float(ev["lower_nM"]) if ev["lower_nM"] else 10.0
        up = float(ev["upper_nM"]) if ev["upper_nM"] else 15000.0
        val_surr = math.sqrt(max(lo, 1.0) * up) if ev["upper_nM"] else lo * 2.0
        cov = ASSAY_COVARIATES["EPIZYME_FLASHPLATE"]
        track_b["train"].append({
            "row_id": f"TB_TR_EPZ_{st['compound_id']}",
            **st,
            "assay_id": "EPIZYME_FLASHPLATE",
            "split": "train",
            "source_family": r["source_family"].split(";")[0],
            "document_id": f"PATFAM_{r['source_family'].split(';')[0]}",
            "publication_year": int(r["source_date_bound"][:4]),
            "sam_uM": cov["sam_uM"],
            "log10_sam_uM": round(math.log10(cov["sam_uM"]), 4),
            "assay_readout": cov["readout"],
            "is_displacement": cov["is_displacement"],
            "binary_y_le100nM": int(r["binary_y_le100nM"]),
            "p_activity": round(9.0 - math.log10(val_surr), 5),
            "relation": "interval",
            "value_nM": round(val_surr, 4),
            "is_exact_numeric": 0,
            "is_weak_gt1000nM": int(lo >= 1000.0),
        })
        tb_seen_train.add(iso_smi)

    # B-Train 2: US9856252_2018 + NUCPAT_FAM
    for aid, fam in [
        ("CHEMBL5733161", "US9856252_2018"),
        ("CHEMBL5734993", "NUCPAT_FAM"),
        ("CHEMBL5737088", "NUCPAT_FAM"),
        ("CHEMBL5737758", "NUCPAT_FAM"),
    ]:
        cov = ASSAY_COVARIATES[aid]
        for iso_smi, rec in aggregate_chembl_assay(aid, quarantine_controls=False).items():
            if iso_smi in tb_seen_train:
                continue
            y = 1 if (rec["relation"] in ("=", "<") and rec["value_nM"] <= 100.0) else (0 if (rec["relation"] in ("=", ">") and rec["value_nM"] > 100.0) else None)
            if y is None:
                continue
            track_b["train"].append({
                "row_id": f"TB_TR_{aid}_{rec['compound_id']}",
                "compound_id": rec["compound_id"],
                "canonical_smiles": rec["canonical_smiles"],
                "nonisomeric_smiles": rec["nonisomeric_smiles"],
                "scaffold_group": rec["scaffold_group"],
                "has_defined_chirality": rec["has_defined_chirality"],
                "assay_id": aid,
                "split": "train",
                "source_family": fam,
                "document_id": rec["document_id"],
                "publication_year": min(rec["publication_year"], 2021),
                "sam_uM": cov["sam_uM"],
                "log10_sam_uM": round(math.log10(cov["sam_uM"]), 4),
                "assay_readout": cov["readout"],
                "is_displacement": cov["is_displacement"],
                "binary_y_le100nM": y,
                "p_activity": round(9.0 - math.log10(rec["value_nM"]), 5),
                "relation": rec["relation"],
                "value_nM": round(rec["value_nM"], 4),
                "is_exact_numeric": int(rec["relation"] == "="),
                "is_weak_gt1000nM": int(rec["value_nM"] > 1000.0 and rec["relation"] in ("=", ">")),
            })
            tb_seen_train.add(iso_smi)

    tb_train_noniso = {r["nonisomeric_smiles"] for r in track_b["train"]}
    tb_train_scaff = {r["scaffold_group"] for r in track_b["train"]}

    # B-Validation (2019-2022)
    cov_pb = ASSAY_COVARIATES["CHEMBL5739093"]
    for r in v1_rows:
        if r["v1_split"] == "validation":
            st = get_std(r["raw_canonical_smiles"])
            ev = json.loads(r["source_evidence"])[0]
            vals = [float(v) for v in ev["values"].values()]
            val_pb = math.exp(sum(math.log(x) for x in vals) / len(vals))
            track_b["validation"].append({
                "row_id": f"TB_VAL_CHEMBL5739093_{st['compound_id']}",
                **st,
                "assay_id": "CHEMBL5739093",
                "split": "validation",
                "source_family": "PHARMABLOCK_2019",
                "document_id": "WO2019173804A1",
                "publication_year": 2019,
                "sam_uM": cov_pb["sam_uM"],
                "log10_sam_uM": round(math.log10(cov_pb["sam_uM"]), 4),
                "assay_readout": cov_pb["readout"],
                "is_displacement": cov_pb["is_displacement"],
                "binary_y_le100nM": int(r["binary_y_le100nM"]),
                "p_activity": round(9.0 - math.log10(val_pb), 5),
                "relation": "=",
                "value_nM": round(val_pb, 4),
                "is_exact_numeric": 1,
                "is_weak_gt1000nM": int(val_pb > 1000.0),
            })

    tb_val_seen = {r["canonical_smiles"] for r in track_b["validation"]}
    for aid, fam in [
        ("CHEMBL5167913", "JMC2022_5154928"),
        ("CHEMBL5738764", "US11524962_2022"),
        ("CHEMBL5738200", "US11365205_2022"),
    ]:
        cov = ASSAY_COVARIATES[aid]
        for iso_smi, rec in aggregate_chembl_assay(aid, quarantine_controls=True).items():
            if iso_smi in tb_val_seen or rec["nonisomeric_smiles"] in tb_train_noniso or rec["scaffold_group"] in tb_train_scaff:
                exclusions.append({"track": "track_b", "assay_id": aid, "smiles": iso_smi, "reason": "cross_split_quarantine"})
                continue
            y = 1 if (rec["relation"] in ("=", "<") and rec["value_nM"] <= 100.0) else (0 if (rec["relation"] in ("=", ">") and rec["value_nM"] > 100.0) else None)
            if y is None:
                continue
            track_b["validation"].append({
                "row_id": f"TB_VAL_{aid}_{rec['compound_id']}",
                "compound_id": rec["compound_id"],
                "canonical_smiles": rec["canonical_smiles"],
                "nonisomeric_smiles": rec["nonisomeric_smiles"],
                "scaffold_group": rec["scaffold_group"],
                "has_defined_chirality": rec["has_defined_chirality"],
                "assay_id": aid,
                "split": "validation",
                "source_family": fam,
                "document_id": rec["document_id"],
                "publication_year": rec["publication_year"],
                "sam_uM": cov["sam_uM"],
                "log10_sam_uM": round(math.log10(cov["sam_uM"]), 4),
                "assay_readout": cov["readout"],
                "is_displacement": cov["is_displacement"],
                "binary_y_le100nM": y,
                "p_activity": round(9.0 - math.log10(rec["value_nM"]), 5),
                "relation": rec["relation"],
                "value_nM": round(rec["value_nM"], 4),
                "is_exact_numeric": int(rec["relation"] == "="),
                "is_weak_gt1000nM": int(rec["value_nM"] > 1000.0 and rec["relation"] in ("=", ">")),
            })
            tb_val_seen.add(iso_smi)

    tb_dev_noniso = tb_train_noniso | {r["nonisomeric_smiles"] for r in track_b["validation"]}
    tb_dev_scaff = tb_train_scaff | {r["scaffold_group"] for r in track_b["validation"]}

    # B-Test (2024-2025 untouched peer-reviewed journal biochemical cohorts with >=20 novel parents)
    tb_test_seen: set[str] = set()
    for aid, fam in [
        ("CHEMBL5544428", "EJMC2024_5543090"),
        ("CHEMBL5509028", "EJMC2024_5506008"),
        ("CHEMBL6088129", "EJMC2025_6087401"),
        ("CHEMBL6090539", "JMC2025_6087464"),
        ("CHEMBL6157149", "JMC2025_6153116"),
    ]:
        cov = ASSAY_COVARIATES[aid]
        for iso_smi, rec in aggregate_chembl_assay(aid, quarantine_controls=True).items():
            if iso_smi in tb_test_seen or rec["nonisomeric_smiles"] in tb_dev_noniso or rec["scaffold_group"] in tb_dev_scaff:
                exclusions.append({"track": "track_b", "assay_id": aid, "smiles": iso_smi, "reason": "cross_split_scaffold_or_control_quarantine"})
                continue
            y = 1 if (rec["relation"] in ("=", "<") and rec["value_nM"] <= 100.0) else (0 if (rec["relation"] in ("=", ">") and rec["value_nM"] > 100.0) else None)
            if y is None:
                continue
            track_b["test"].append({
                "row_id": f"TB_TE_{aid}_{rec['compound_id']}",
                "compound_id": rec["compound_id"],
                "canonical_smiles": rec["canonical_smiles"],
                "nonisomeric_smiles": rec["nonisomeric_smiles"],
                "scaffold_group": rec["scaffold_group"],
                "has_defined_chirality": rec["has_defined_chirality"],
                "assay_id": aid,
                "split": "test",
                "source_family": fam,
                "document_id": rec["document_id"],
                "publication_year": rec["publication_year"],
                "sam_uM": cov["sam_uM"],
                "log10_sam_uM": round(math.log10(cov["sam_uM"]), 4),
                "assay_readout": cov["readout"],
                "is_displacement": cov["is_displacement"],
                "binary_y_le100nM": y,
                "p_activity": round(9.0 - math.log10(rec["value_nM"]), 5),
                "relation": rec["relation"],
                "value_nM": round(rec["value_nM"], 4),
                "is_exact_numeric": int(rec["relation"] == "="),
                "is_weak_gt1000nM": int(rec["value_nM"] > 1000.0 and rec["relation"] in ("=", ">")),
            })
            tb_test_seen.add(iso_smi)

    # =========================================================================
    # VERIFY ZERO LEAKAGE ON BOTH TRACKS
    # =========================================================================
    audit_a = verify_zero_leakage("Track_A", track_a, sk_quarantine)
    audit_b = verify_zero_leakage("Track_B", track_b, sk_quarantine)

    # Separate inputs (`structures.csv`), development labels (`dev_labels.csv`), and sealed test labels (`sealed_test_labels/`)
    ta_struct_cols = [
        "row_id", "compound_id", "split", "canonical_smiles", "nonisomeric_smiles",
        "scaffold_group", "has_defined_chirality", "source_family", "document_id",
        "publication_year", "assay_pair_type", "is_direct_paired",
    ]
    ta_label_cols = [
        "row_id", "compound_id", "split", "source_family", "document_id",
        "y_coop_3x", "y_coop_10x", "delta_p_activity", "mta_rel", "mta_value_nM", "sam_rel", "sam_value_nM",
    ]
    tb_struct_cols = [
        "row_id", "compound_id", "assay_id", "split", "canonical_smiles", "nonisomeric_smiles",
        "scaffold_group", "has_defined_chirality", "source_family", "document_id",
        "publication_year", "sam_uM", "log10_sam_uM", "assay_readout", "is_displacement",
    ]
    tb_label_cols = [
        "row_id", "compound_id", "assay_id", "split", "source_family", "document_id",
        "binary_y_le100nM", "p_activity", "relation", "value_nM", "is_exact_numeric", "is_weak_gt1000nM",
    ]

    ta_all = track_a["train"] + track_a["validation"] + track_a["test"]
    tb_all = track_b["train"] + track_b["validation"] + track_b["test"]

    write_csv(root / "data/processed/track_a_structures.csv", [{k: r[k] for k in ta_struct_cols} for r in ta_all])
    write_csv(root / "data/processed/track_a_dev_labels.csv", [{k: r[k] for k in ta_label_cols} for r in track_a["train"] + track_a["validation"]])
    write_csv(root / "data/processed/sealed_test_labels/track_a_test_labels.csv", [{k: r[k] for k in ta_label_cols} for r in track_a["test"]])

    write_csv(root / "data/processed/track_b_structures.csv", [{k: r[k] for k in tb_struct_cols} for r in tb_all])
    write_csv(root / "data/processed/track_b_dev_labels.csv", [{k: r[k] for k in tb_label_cols} for r in track_b["train"] + track_b["validation"]])
    write_csv(root / "data/processed/sealed_test_labels/track_b_test_labels.csv", [{k: r[k] for k in tb_label_cols} for r in track_b["test"]])

    write_csv(root / "data/processed/exclusions_v2.csv", exclusions)

    # Summarize counts & class prevalences for M0 manifests
    m0_summary = {
        "protocol_id": protocol["protocol_id"],
        "chembl_release": "ChEMBL_37 (2026-05-01)",
        "track_a": {
            **audit_a,
            "class_counts": {
                sp: {
                    "n": len(track_a[sp]),
                    "direct_paired": sum(r["is_direct_paired"] for r in track_a[sp]),
                    "y_coop_3x_pos": sum(r["y_coop_3x"] for r in track_a[sp]),
                    "y_coop_10x_pos": sum(r["y_coop_10x"] for r in track_a[sp]),
                    "chiral_parents": sum(r["has_defined_chirality"] == "True" for r in track_a[sp]),
                }
                for sp in ("train", "validation", "test")
            },
            "test_by_source": {
                fam: {
                    "n": sum(1 for r in track_a["test"] if r["source_family"] == fam),
                    "scaffolds": len({r["scaffold_group"] for r in track_a["test"] if r["source_family"] == fam}),
                    "y_coop_3x_pos": sum(r["y_coop_3x"] for r in track_a["test"] if r["source_family"] == fam),
                    "y_coop_10x_pos": sum(r["y_coop_10x"] for r in track_a["test"] if r["source_family"] == fam),
                }
                for fam in sorted({r["source_family"] for r in track_a["test"]})
            },
        },
        "track_b": {
            **audit_b,
            "class_counts": {
                sp: {
                    "n": len(track_b[sp]),
                    "pos_le100nM": sum(r["binary_y_le100nM"] for r in track_b[sp]),
                    "weak_gt1000nM": sum(r["is_weak_gt1000nM"] for r in track_b[sp]),
                    "exact_numeric": sum(r["is_exact_numeric"] for r in track_b[sp]),
                    "chiral_parents": sum(r["has_defined_chirality"] == "True" for r in track_b[sp]),
                }
                for sp in ("train", "validation", "test")
            },
            "test_by_assay": {
                aid: {
                    "source_family": next(r["source_family"] for r in track_b["test"] if r["assay_id"] == aid),
                    "n": sum(1 for r in track_b["test"] if r["assay_id"] == aid),
                    "scaffolds": len({r["scaffold_group"] for r in track_b["test"] if r["assay_id"] == aid}),
                    "pos_le100nM": sum(r["binary_y_le100nM"] for r in track_b["test"] if r["assay_id"] == aid),
                    "weak_gt1000nM": sum(r["is_weak_gt1000nM"] for r in track_b["test"] if r["assay_id"] == aid),
                }
                for aid in sorted({r["assay_id"] for r in track_b["test"]})
            },
        },
        "unique_molecules_across_both_tracks": len({r["canonical_smiles"] for r in ta_all} | {r["canonical_smiles"] for r in tb_all}),
        "exclusions_logged": len(exclusions),
    }

    m0_path = root / "manifests/m0_source_audit.json"
    m0_path.write_text(json.dumps(m0_summary, indent=2) + "\n", encoding="utf-8")

    # Write cohort_split_manifest.json binding every processed file hash
    bound_files = [
        "configs/protocol_v2.yaml",
        "manifests/pdb_receptors.json",
        "manifests/v1_sk_quarantine.json",
        "manifests/m0_source_audit.json",
        "data/raw_curated/chembl_selected_activities.csv",
        "data/raw_curated/mirati_table9a_unique.csv",
        "data/raw_curated/v1_pre2020_cohorts.csv",
        "data/processed/track_a_structures.csv",
        "data/processed/track_a_dev_labels.csv",
        "data/processed/sealed_test_labels/track_a_test_labels.csv",
        "data/processed/track_b_structures.csv",
        "data/processed/track_b_dev_labels.csv",
        "data/processed/sealed_test_labels/track_b_test_labels.csv",
        "data/processed/exclusions_v2.csv",
        "scripts/chemistry_v2.py",
        "scripts/prepare_receptors.py",
        "scripts/build_cohorts_m0.py",
    ]
    split_manifest = {
        "protocol_id": protocol["protocol_id"],
        "bound_sha256": {rel: sha256_file(root / rel) for rel in bound_files},
    }
    split_manifest_path = root / "manifests/cohort_split_manifest.json"
    split_manifest_path.write_text(json.dumps(split_manifest, indent=2) + "\n", encoding="utf-8")

    # Write Stage 1 Freeze Receipt (`receipts/protocol_freeze.json`)
    freeze_receipt = {
        "stage": "STAGE_1_M0_PROTOCOL_FREEZE",
        "status": "FROZEN_BEFORE_FEATURIZATION_AND_MODEL_FITTING",
        "frozen_at_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "models_fitted": False,
        "molecular_features_computed": False,
        "test_labels_accessed_by_models": False,
        "v1_sk_test_overlap": 0,
        "cross_split_leakage": 0,
        "protocol_sha256": sha256_file(protocol_path),
        "pdb_receptors_manifest_sha256": sha256_file(root / "manifests/pdb_receptors.json"),
        "cohort_split_manifest_sha256": sha256_file(split_manifest_path),
        "bound_sha256": split_manifest["bound_sha256"],
    }
    receipt_path = root / "receipts/protocol_freeze.json"
    receipt_path.parent.mkdir(parents=True, exist_ok=True)
    receipt_path.write_text(json.dumps(freeze_receipt, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(m0_summary, indent=2))
    print(f"Sealed Stage 1 protocol freeze receipt at {receipt_path}")


if __name__ == "__main__":
    main()
