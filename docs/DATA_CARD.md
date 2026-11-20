# Data Card — PRMT5 v2 (`ChEMBL_37` + Original Patent Tables + RCSB PDB Crystal Ensemble)

## 1. Provenance & Licensing

- **Bioactivity Sources:** Derived from public `ChEMBL_37` (`2026-05-01` release; human PRMT5 `CHEMBL1795116` and human PRMT5/MEP50 complex `CHEMBL3137261`), cross-audited against original patent biochemical tables for Mirati family `72659905` (`US20210078994A1` / `US11479551B2`, earliest verified public table disclosure `2021-03-18`), Epizyme families (`50029211`, `50001259`, `50001260`), and Pharmablock (`WO2019173804A1` / `US11591326B2`). ChEMBL-derived data tables are distributed under [CC BY-SA 3.0](https://creativecommons.org/licenses/by-sa/3.0/).
- **3D Macromolecular Structures:** Experimental co-crystal structures retrieved from the [RCSB Protein Data Bank](https://www.rcsb.org/pages/usage-policy) under CC0 1.0 Public Domain Dedication (`.licenses/pdb_database_LICENSE.txt`; exact SHA-256 hashes in [`manifests/pdb_receptors.json`](../manifests/pdb_receptors.json)).

---

## 2. Dataset Summary & Split Isolation

All molecules are standardized via `scripts/chemistry_v2.py` (largest organic fragment selection, charge neutralization, RDKit canonical tautomerization, preserved stereochemistry in `canonical_smiles`, stereochemistry-stripped `nonisomeric_smiles`, and Bemis-Murcko `scaffold_group`). Across all split pairs in both tracks, overlap is strictly **`0`** on `canonical_smiles`, `nonisomeric_smiles`, `scaffold_group`, and `source_family`, and **`0`** against v1's quarantined 246 SK test compounds (`manifests/v1_sk_quarantine.json`).

### Track A: Cofactor-Cooperativity Cohort (`PRMT5•MTA` vs `PRMT5•SAM`)

| Split | Years | Source Families | Unique Parents | Murcko Scaffolds | Chiral Parents | Direct Paired | `y_coop_3x` Positives | `y_coop_10x` Positives |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| **Train** | `≤2021` | 5 (`MIRATI_72659905`, 3 Epizyme, `NUCPAT_FAM`) | 430 | 159 | 131 (30.5%) | 190 | 172 (40.0%) | 147 (34.2%) |
| **Validation** | `2019–2022` | 4 (`TNG_BMC2022`, `CALU6_JMC2022`, `PHARMABLOCK_2019`, `JMC2022_5154928`) | 88 | 55 | 29 (33.0%) | 8 | 18 (20.5%) | 16 (18.2%) |
| **Held-Out Test** | `2024–2025` | 4 (`MECHSWITCH_JMC2025`, `TNG908_JMC2024`, `AZ_JMC2024`, `AM9747_JMC2025`) | **98** | **49** | **69 (70.4%)** | **98 (100%)** | **44 (44.9%)** | **21 (21.4%)** |

**Held-Out Track A Test Breakdown by Peer-Reviewed *J. Med. Chem.* Study:**
- `MECHSWITCH_JMC2025` (`CHEMBL6153116`, PMID `39919252`, 2025): `55` paired parents (`29` scaffolds), fluorescence anisotropy peptide displacement `Ki` in presence of `MTA` (`CHEMBL6157148`) vs `SAM` (`CHEMBL6157149`); `20` `y_coop_3x` positives, `14` `y_coop_10x` positives, spanning SAM-preferring THIQ leads ($\Delta\text{p}K_i = -2.68$) to rational MTA-switched leads ($\Delta\text{p}K_i = +1.88$).
- `TNG908_JMC2024` (`CHEMBL5543141`, PMID `38595098`, 2024): `22` paired parents (`10` scaffolds), `PRMT5/MEP50` peptide displacement `IC50` with `MTA` (`CHEMBL5547645`) vs `SAM` (`CHEMBL5547646`); `8` `y_coop_3x` positives, `2` `y_coop_10x` positives.
- `AZ_JMC2024` (`CHEMBL5576120`, PMID `39080842`, 2024): `14` paired parents (`6` scaffolds), `PRMT5/MEP50` MTase-Glo `IC50` plus `MTA` (`CHEMBL5577992`) vs minus `MTA` (`CHEMBL5577993`); `11` `y_coop_3x` positives, `4` `y_coop_10x` positives.
- `AM9747_JMC2025` (`CHEMBL6115040`, PMID `40102181`, 2025): `7` paired parents (`4` scaffolds), `PRMT5/MEP50` MTase-Glo `IC50` with `SAM/MTA` (`CHEMBL6119736`) vs `SAM` (`CHEMBL6119737`); `5` `y_coop_3x` positives, `1` `y_coop_10x` positive.

---

### Track B: Cross-Document Biochemical Potency Cohort (`2024–2025` Untouched Series)

| Split | Years | Source Families | Unique Parents | Murcko Scaffolds | Chiral Parents | Active `≤100 nM` | Weak/Inactive `>1 µM` | Exact Numeric |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| **Train** | `≤2021` | 5 (3 Epizyme, `US9856252_2018`, `NUCPAT_FAM`) | 818 | 361 | 257 (31.4%) | 240 (29.3%) | 276 (33.7%) | 219 |
| **Validation** | `2019–2022` | 4 (`PHARMABLOCK_2019`, `JMC2022_5154928`, `US11524962_2022`, `US11365205_2022`) | 261 | 119 | 126 (48.3%) | 129 (49.4%) | 20 (7.7%) | 259 |
| **Held-Out Test** | `2024–2025` | 5 (`EJMC2024_5543090`, `EJMC2024_5506008`, `EJMC2025_6087401`, `JMC2025_6087464`, `JMC2025_6153116`) | **182** | **75** | **140 (76.9%)** | **57 (31.3%)** | **87 (47.8%)** | **136** |

**Held-Out Track B Test Breakdown by Peer-Reviewed Assay:**
- `CHEMBL5544428` (`EJMC2024_5543090`, PMID `38301329`, AlphaLISA): `35` novel ring-opened nucleoside parents (`12` scaffolds, `11 ≤100 nM`, `2 >1 µM`).
- `CHEMBL5509028` (`EJMC2024_5506008`, PMID `38518523`, SPA radiometric): `42` novel tetrahydro-$\beta$-carboline parents (`4` scaffolds, `0 ≤100 nM`, `42 >1 µM` confirmed weak/inactive series).
- `CHEMBL6088129` (`EJMC2025_6087401`, PMID `39488968`, AlphaLISA): `29` novel piperazine-nucleoside parents (`10` scaffolds, `21 ≤100 nM`, `0 >1 µM`).
- `CHEMBL6090539` (`JMC2025_6087464`, PMID `39722476`, `[3H]-SAM` enzymatic): `21` novel dihydroisoquinolinone parents (`20` scaffolds, `14 ≤100 nM`, `4 >1 µM`).
- `CHEMBL6157149` (`JMC2025_6153116`, PMID `39919252`, `PRMT5/MEP50` + `SAM` `Ki`): `55` novel parents (`29` scaffolds, `11 ≤100 nM`, `39 >1 µM`).

---

## 3. Experimental Crystal Structures (`manifests/pdb_receptors.json`)

All structures are aligned onto `PDB 4X61` chain A catalytic domain $C_\alpha$ atoms (`residues 300–630`, `331` $C_\alpha$ atoms, RMSD `0.41–0.77 Å`):
- **Pre-2023 Featurization & Docking Ensemble (`receptors_pre2023/`):**
  - *State 1 (`PRMT5•SAM` + THIQ substrate channel):* `4X61` (`EPZ015666`, `2.85 Å`), `5EML` (`Cmpd 10`, `2.39 Å`), `7L1G` (binary `SAM`, `2.47 Å`).
  - *State 2 (`PRMT5•MTA` cooperative induced pocket):* `7S1S` (`MRTX1719`, `2.62 Å`), `8CSG` (`Fragment 1`, `2.48 Å`).
  - *State 3 (`PRMT5` SAM-pocket nucleoside mimetics):* `6CKC` (`LLY-283`, `2.80 Å`), `7MX7` (`PF-06939999`, `2.49 Å`), `6RLL` (`JNJ-44064146`, `2.22 Å`).
- **Post-2023 Quarantined Diagnostic Co-Crystals (`diagnostic_post2023/`, never used in featurization/modeling):** `8VEO` (`MTA`, `2.03 Å`), `8VEY` (`TNG908`, `2.44 Å`), `9C10` (`AMG 193`, `2.85 Å`), `9N3R` (`TNG462`, `2.47 Å`), `9T43` (`AZD3470`, `2.18 Å`).

## 4. Known Limitations

1. Track A training negatives combine `Mirati` paired low-selectivity compounds with active (`≤100 nM`) pre-2022 State 1 (`PRMT5•SAM`) and State 3 (`PRMT5` nucleoside) inhibitors whose mechanism precludes MTA cooperativity, whereas **100% (`98/98`) of Track A test labels** are direct head-to-head experimental `PRMT5•MTA` vs `PRMT5•SAM` biochemical measurements.
2. Assay `[SAM]` concentrations for literature ChEMBL entries where the ChEMBL description omits exact micromolar concentration default to standard category protocol metadata (`1 µM` enzymatic / `10 µM` competitive peptide displacement); exact row-level patent concentrations are retained wherever disclosed.
3. Biochemical potency and cofactor cooperativity do not establish *in vivo* pharmacokinetics, safety, or clinical efficacy.
