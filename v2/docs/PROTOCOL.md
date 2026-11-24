# PRMT5 v2 Scientific Protocol (`PRMT5-V2-3D-COOPERATIVITY-2026`)

**Machine-readable protocol freeze:** [`configs/protocol_v2.yaml`](../configs/protocol_v2.yaml)  
**Stage 1 cryptographic freeze receipt:** [`receipts/protocol_freeze.json`](../receipts/protocol_freeze.json)  
**Stage 2 prediction seal receipt:** [`receipts/prediction_seal.json`](../receipts/prediction_seal.json)

---

## 1. Motivation: Why PRMT5 SourceBench v1 Failed

In PRMT5 SourceBench v1 (`prmt5-sourcebench`), a partially fine-tuned non-chiral 2D `MoLFormer` model trained on 3 early Epizyme patent families (`604` parents, `75 nM` SAM) and tuned on Pharmablock (`70` parents, `10 µM` SAM, `55.7%` positive) failed to transfer to SK (`246` parents, `1 µM` SAM), scoring **AP = 0.2553** below the constant-prevalence baseline (**0.2805**) and recovering only **2 of 69** positives at the validation-locked decision threshold (`2.9%` sensitivity). Three root causes drove that collapse:

1. **Lead-Optimization Steric Cliffs & Chirality Loss:** 243 of 246 SK molecules already shared the active `(S)`-hydroxypropyl-tetrahydroisoquinoline (`THIQ`) warhead (`69 +++` `≤100 nM` vs `174 ++` `100–1000 nM`). Stripping stereochemistry (`includeChirality=False` / nonisomeric SMILES) and omitting 3D active-site pocket geometry left 2D models blind to chiral orientation and R-group steric clashes.
2. **133-fold `[SAM]` Shift × Cofactor Kinetics:** Locking a single raw probability cutoff on `10 µM`-SAM Pharmablock (`55.7%` positive) caused severe miscalibration when deployed on lower-prevalence, shifted-`[SAM]` assays (`75 nM` to `10 µM` SAM).
3. **Clinical Paradigm Shift to 2nd-Generation `PRMT5•MTA` Inhibitors:** First-generation `PRMT5•SAM` inhibitors caused clinical dose-limiting hematologic toxicity because healthy bone marrow maintains high intracellular `[SAM]`. Modern precision oncology pivoted to **MTA-cooperative `PRMT5•MTA` selective inhibitors** (`MRTX1719` / `BMS-986504`, `TNG908`, `TNG462`, `AMG 193`, `AZD3470`) targeting `MTAP`-deleted tumors where `5'-methylthioadenosine (MTA)` accumulates (`2–10 µM`). Selectivity is governed by a **3D dual steric switch** inside the PRMT5 catalytic pocket.

---

## 2. Sacred Non-Leakage & Historical Quarantine Rules

1. **Strict Quarantine of v1's 246 SK Test Compounds:** v1's 246 SK test molecules (`CHEMBL5737316`, family `75338422`, `US20210101892A1`) are **never** re-evaluated, tuned against, or included in any PRMT5 v2 split. SHA-256 digests of all 246 SK isomeric SMILES, nonisomeric SMILES, and Bemis-Murcko scaffolds are recorded in [`manifests/v1_sk_quarantine.json`](../manifests/v1_sk_quarantine.json), with `0` overlap enforced across all v2 splits.
2. **Four-Axis Zero Cross-Split Leakage:** Across `(train, validation)`, `(train, test)`, and `(validation, test)` in both benchmark tracks, zero overlap is enforced simultaneously on:
   - `canonical_smiles` (standardized isomeric parent SMILES)
   - `nonisomeric_smiles` (stereochemistry-stripped parent SMILES, preventing cross-split enantiomer/diastereomer leakage)
   - `scaffold_group` (canonical Bemis-Murcko ring-linker scaffold)
   - `source_family` (patent family or peer-reviewed document cohort)
3. **Historical Reference Control Quarantine:** Named literature positive controls (`EPZ015666`, `GSK3326595`, `JNJ-64619178`, `LLY-283`, `MRTX1719`, `TNG908`, `AM-9747`) re-tested in 2022–2025 medicinal chemistry tables are quarantined from validation and test cohorts (`data/processed/exclusions_v2.csv`).
4. **Structural Crystal Temporal Cutoff (`< 2023-01-01`):** Only experimental RCSB PDB co-crystal structures released before `2023-01-01` (`4X61`, `5EML`, `7L1G`, `7S1S`, `8CSG`, `6CKC`, `7MX7`, `6RLL`) are permitted in Stage 2 featurization and docking. Post-2023 clinical co-crystals (`8VEO`, `8VEY` [`TNG908`], `9C10` [`AMG 193`], `9N3R` [`TNG462`], `9T43` [`AZD3470`]) are quarantined in `data/structures/aligned/diagnostic_post2023/` for Stage 3 biophysical validation figures only.
5. **Physical Isolation of Sealed Test Labels:** Test ground-truth targets reside exclusively in `data/processed/sealed_test_labels/` and are never read by `featurize_3d.py` or `train_and_seal.py`. All Stage 2 test predictions are cryptographically sealed in `receipts/prediction_seal.json` (`Commit 2`) prior to Stage 3 scoring (`Commit 3`).

---

## 3. Two Confirmatory Benchmark Tracks

### Track A — Cofactor-Cooperativity & `MTAP`-Selective Mechanism Prediction (`PRMT5•MTA` vs `PRMT5•SAM`)
- ** Cohort:** `430` train (`≤2021`, `159` scaffolds, `5` families) / `88` validation (`2019–2022`, `55` scaffolds, `4` families) / `98` test (`2024–2025`, `49` scaffolds, `4` untouched *J. Med. Chem.* series: `MECHSWITCH_JMC2025` [`n=55`], `TNG908_JMC2024` [`n=22`], `AZ_JMC2024` [`n=14`], `AM9747_JMC2025` [`n=7`]). Every single test molecule has direct head-to-head experimental `PRMT5•MTA` and `PRMT5•SAM` biochemical measurements.
- **Endpoints:**
  - **Primary binary (`y_coop_3x`):** $\Delta\text{pActivity} = \log_{10}(\text{Activity}_{\text{SAM}} / \text{Activity}_{\text{MTA}}) \ge 0.5$ ($\ge 3.16\times$ MTA cooperativity) with $\text{Activity}_{\text{MTA}} < 10\ \mu\text{M}$ (`44/98 = 44.9%` test prevalence); evaluated by **Average Precision (AP)**, **ROC-AUC**, **Top-20 Screening Hits**, and **Balanced Accuracy** at validation-locked threshold.
  - **Co-primary stringent clinical selectivity (`y_coop_10x`):** $\Delta\text{pActivity} \ge 1.0$ ($\ge 10\times$ MTA selectivity) with $\text{Activity}_{\text{MTA}} \le 1\ \mu\text{M}$ (`21/98 = 21.4%` test prevalence).
  - **Continuous selectivity ordering:** Pooled Spearman $\rho$ and censor-aware Concordance Index ($C$-index) against empirical $\Delta\text{pActivity}$.
- **Methods Compared:** `B0` (Constant reference), `P0` (**Zero-Shot 3D Cofactor-State Pocket Physics** — zero training labels used), `B1` (v1-style achiral 2D Morgan classifier), `M1` (Chiral 2D + physicochemical ensemble), `M2` (Chiral 2D + 3D multi-state cofactor-pocket physics hybrid).

### Track B — Assay-Conditioned Cross-Document Potency Ranking & Calibration (`2024–2025` Series)
- **Cohort:** `818` train (`≤2021`, `361` scaffolds, `5` families, `29.3%` `≤100 nM`) / `261` validation (`2019–2022`, `119` scaffolds, `4` families, `49.4%` `≤100 nM`) / `182` test (`2024–2025`, `75` scaffolds, `5` untouched peer-reviewed journal biochemical cohorts: `CHEMBL5544428`, `CHEMBL5509028`, `CHEMBL6088129`, `CHEMBL6090539`, `CHEMBL6157149`, `31.3%` `≤100 nM`).
- **Endpoints:**
  - **Primary:** Within-assay pair-weighted censor-aware **Concordance Index ($C$-index)**, **Pooled Average Precision (AP)**, and **Pooled ROC-AUC**.
  - **Secondary:** **Strict Virtual-Screening Enrichment AP** (`≤100 nM` actives [`n=57`] vs confirmed `>1 µM` weak/inactive compounds [`n=87`]), **Top-25 Hit Recovery**, and **Threshold Sensitivity / Balanced Accuracy** across shifted `[SAM]` assays.
- **Methods Compared:** `B0` (Constant reference), `B1` (v1-style uncalibrated non-chiral 2D classifier), `M1` (Chiral 2D + physicochemical + assay covariates within-source pairwise ranker + OOF base-rate/covariate calibration), `M2` (Chiral 2D + 3D multi-state pocket physics + within-source pairwise ranker).

---

## 4. Statistical Uncertainty

All test metrics and pairwise method contrasts (`M2 − B1`, `M2 − M1`, `M1 − B1`, `P0 − B1`) are accompanied by **2,000 paired Murcko-scaffold cluster bootstrap resamples** (`seed = 20261120`), resampling test Murcko scaffolds with replacement and retaining all members of each sampled scaffold simultaneously across all compared models.
