#!/usr/bin/env python3
"""Extract minimal, auditable raw source tables into data/raw_curated/ for self-contained M0 cohort building."""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from chemistry_v2 import standardize_structure

RECON = Path("/home/seedna2/.gemini/antigravity/brain/5943980c-5b00-4b66-a5a5-88c5cedf77cb/scratch/m0_recon")
EXT = RECON / "v1_extract/prmt5-discovery"
V1_COHORT = Path("/home/seedna2/.gemini/antigravity/brain/c94c7cdc-dbd8-40be-8165-c1c697fb97c6/scratch/repos/prmt5-sourcebench/data/processed/cohort.csv")

INCLUDED_CHEMBL_ASSAYS = {
    # Track B Train & Track A Negatives (pre-2022)
    "CHEMBL5733161",  # US-9856252-B2 (2018, 1 uM SAM)
    "CHEMBL5734993",  # US-10570140-B2 (2020, nucleoside SAM-pocket)
    "CHEMBL5737088",  # US-11078205-B2 (2021, nucleoside SAM-pocket)
    "CHEMBL5737758",  # US-11220524-B2 (continuation family of US-10570140 nucleoside series)
    # Track A & Track B Validation (2022)
    "CHEMBL5109533",  # BMC 2022 TNG phthalazinones MTAP-KO
    "CHEMBL5109534",  # BMC 2022 TNG phthalazinones MTAP-WT
    "CHEMBL5130634",  # JMC 2022 Calu-6 MTAP-/-
    "CHEMBL5167913",  # JMC 2022 [3H]-SAM displacement biochemical
    "CHEMBL5738764",  # US-11524962-B2 (2022 FlashPlate biochemical)
    "CHEMBL5738200",  # US-11365205-B2 (2022 TR-FRET biochemical)
    # Track A Test (2024-2025 paired PRMT5•MTA vs PRMT5•SAM biochemical JMC cohorts)
    "CHEMBL5547645",  # JMC 2024 TNG908 series +MTA
    "CHEMBL5547646",  # JMC 2024 TNG908 series +SAM
    "CHEMBL5577992",  # JMC 2024 AZ-PRMT5i-1 series +MTA
    "CHEMBL5577993",  # JMC 2024 AZ-PRMT5i-1 series -MTA
    "CHEMBL6157148",  # JMC 2025 Mechanism Switching series +MTA Ki
    "CHEMBL6157149",  # JMC 2025 Mechanism Switching series +SAM Ki (also Track B Test)
    "CHEMBL6119736",  # JMC 2025 Amgen AM-9747 DEL series +SAM/MTA
    "CHEMBL6119737",  # JMC 2025 Amgen AM-9747 DEL series +SAM
    # Track B Test (2023-2025 peer-reviewed journal biochemical cohorts)
    "CHEMBL5623742",  # EJMC 2023 THIQ-indole biochemical
    "CHEMBL5544428",  # EJMC 2024 ring-opening nucleoside AlphaLISA
    "CHEMBL5509028",  # EJMC 2024 tetrahydro-beta-carboline SPA
    "CHEMBL6088129",  # EJMC 2025 piperazine nucleoside AlphaLISA
    "CHEMBL6090539",  # JMC 2025 dihydroisoquinolinone biochemical
}


def sha_str(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    out_dir = root / "data/raw_curated"
    out_dir.mkdir(parents=True, exist_ok=True)

    docs = {}
    for p in sorted((EXT / "data/raw").glob("documents_*.json")):
        for d in json.loads(p.read_text()).get("documents", []):
            docs[d["document_chembl_id"]] = d

    ledger = list(csv.DictReader((RECON / "v1_activity_ledger.csv").open()))
    ledger_by_aid = {r["activity_id"]: r for r in ledger}

    # 1. Curated ChEMBL activity slice
    chembl_rows = []
    for r in ledger:
        if r["assay_chembl_id"] in INCLUDED_CHEMBL_ASSAYS and r["canonical_smiles"]:
            d = docs.get(r["document_chembl_id"], {})
            chembl_rows.append({
                "activity_id": r["activity_id"],
                "molecule_chembl_id": r["molecule_chembl_id"],
                "assay_chembl_id": r["assay_chembl_id"],
                "document_chembl_id": r["document_chembl_id"],
                "document_year": d.get("year") or r["document_year"],
                "doc_type": d.get("doc_type") or "",
                "journal_or_patent": d.get("patent_id") or d.get("journal") or "",
                "pubmed_id": str(d.get("pubmed_id") or ""),
                "standard_type": r["standard_type"],
                "standard_relation": r["standard_relation"],
                "standard_value": r["standard_value"],
                "standard_units": r["standard_units"],
                "data_validity_comment": r["data_validity_comment"],
                "raw_canonical_smiles": r["canonical_smiles"],
                "assay_format": r["assay_format"],
                "confidence_score": r["confidence"],
                "assay_description": r["assay_description"],
            })
    chembl_rows.sort(key=lambda x: int(x["activity_id"]))
    with (out_dir / "chembl_selected_activities.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(chembl_rows[0]))
        w.writeheader()
        w.writerows(chembl_rows)

    # 2. Mirati Table 9a mapped rows (verified 2021-03-18 disclosure)
    mir_rows = list(csv.DictReader((EXT / "data/curation/mirati/mirati_activity_mapping.csv").open()))
    mir_out = []
    for mr in mir_rows:
        if mr["mapping_status"] == "unique" and mr["protocol"] == "FlashPlate_H4":
            lr = ledger_by_aid.get(mr["activity_id"])
            if lr and lr["canonical_smiles"]:
                mir_out.append({
                    "activity_id": mr["activity_id"],
                    "molecule_chembl_id": mr["molecule_chembl_id"],
                    "patent": mr["patent"],
                    "family_id": mr["family_id"],
                    "earliest_verified_table_publication": mr["earliest_verified_table_publication"],
                    "table": mr["table"],
                    "protocol": mr["protocol"],
                    "source_example": mr["source_example"],
                    "mta_condition": mr["mta_condition"],
                    "relation": mr["relation"],
                    "value_nM": mr["value_nM"],
                    "raw_canonical_smiles": lr["canonical_smiles"],
                })
    mir_out.sort(key=lambda x: int(x["activity_id"]))
    with (out_dir / "mirati_table9a_unique.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(mir_out[0]))
        w.writeheader()
        w.writerows(mir_out)

    # 3. v1 pre-2020 cohorts (Epizyme train 604 + Pharmablock validation 70) AND sealed SK 246 quarantine hashes
    v1_rows = list(csv.DictReader(V1_COHORT.open()))
    v1_pre2020 = []
    sk_quarantine = {"iso_sha256": [], "noniso_sha256": [], "scaffold_sha256": [], "count": 0, "source_family": "75338422", "assay_id": "CHEMBL5737316"}
    for r in v1_rows:
        st = standardize_structure(r["canonical_smiles"])
        if r["split"] in ("train", "validation"):
            v1_pre2020.append({
                "v1_compound_id": r["compound_id"],
                "v1_split": r["split"],
                "raw_canonical_smiles": r["canonical_smiles"],
                "binary_y_le100nM": r["binary_y"],
                "source_family": r["source_family"],
                "assay_ids": r["assay_ids"],
                "source_date_bound": r["source_date_bound"],
                "source_evidence": r["source_evidence"],
            })
        elif r["split"] == "test":
            sk_quarantine["count"] += 1
            sk_quarantine["iso_sha256"].append(sha_str(st["canonical_smiles"]))
            sk_quarantine["noniso_sha256"].append(sha_str(st["nonisomeric_smiles"]))
            sk_quarantine["scaffold_sha256"].append(sha_str(st["scaffold_group"]))

    sk_quarantine["iso_sha256"] = sorted(set(sk_quarantine["iso_sha256"]))
    sk_quarantine["noniso_sha256"] = sorted(set(sk_quarantine["noniso_sha256"]))
    sk_quarantine["scaffold_sha256"] = sorted(set(sk_quarantine["scaffold_sha256"]))

    with (out_dir / "v1_pre2020_cohorts.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(v1_pre2020[0]))
        w.writeheader()
        w.writerows(v1_pre2020)

    (root / "manifests/v1_sk_quarantine.json").write_text(json.dumps(sk_quarantine, indent=2) + "\n", encoding="utf-8")
    print(f"Extracted {len(chembl_rows)} ChEMBL rows, {len(mir_out)} Mirati T9a rows, {len(v1_pre2020)} v1 pre-2020 rows, and {sk_quarantine['count']} SK quarantine hashes.")


if __name__ == "__main__":
    main()
