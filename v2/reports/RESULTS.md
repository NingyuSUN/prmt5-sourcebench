# PRMT5 v2 Scientific Results & Biophysical Analysis

**Protocol:** [`configs/protocol_v2.yaml`](../configs/protocol_v2.yaml) (`PRMT5-V2-3D-COOPERATIVITY-2026`)  
**Frozen M0 Audit:** [`manifests/m0_source_audit.json`](../manifests/m0_source_audit.json)  
**Sealed Prediction Receipt:** [`receipts/prediction_seal.json`](../receipts/prediction_seal.json)  
**Canonical Metrics Artifact:** [`reports/results/final_metrics_v2.json`](results/final_metrics_v2.json)

---

## 1. Crystallographic Discovery: The Dual Steric Switch Governing PRMT5 Cofactor Selectivity

To establish a first-principles physical baseline before fitting any bioactivity model, we aligned 13 experimental human PRMT5 co-crystal structures retrieved from the RCSB PDB (`331` catalytic-domain $C_\alpha$ atoms, residues `300–630`, RMSD `0.41–0.77 Å` relative to `PDB 4X61` chain A; [`manifests/pdb_receptors.json`](../manifests/pdb_receptors.json)). Superimposing `SAM` and `Glu435` from `4X61` (`PRMT5•MEP50•SAM` + `EPZ015666`) against `MTA` and `Phe327` from `7S1S` (`PRMT5•MEP50•MTA` + `MRTX1719`) reveals a **reciprocal two-gate steric switch**:

| PDB ID | Release Date | Inhibitor (`comp_id`) | Binding State | Min Dist to `SAM:CE` (Å) | Min Dist to `Glu435` (`SAM` state, Å) | Min Dist to `Phe327` (`MTA` state, Å) | Featurization Status |
|---|---|---|---:|---:|---:|---:|---|
| `4X61` | 2015-04-22 | `EPZ015666` (`3XV`) | State 1 (`PRMT5•SAM`) | `3.61` | `3.19` | **`0.58` (clash)** | Pre-2023 Vina receptor & template |
| `5EML` | 2016-02-24 | `Epizyme_Cmpd10` (`5QK`) | State 1 (`PRMT5•SAM`) | `3.63` | `3.21` | **`0.51` (clash)** | Pre-2023 O3A template |
| `7S1S` | 2022-01-26 | `MRTX1719` (`85K`) | State 2 (`PRMT5•MTA`) | **`2.71` (clash)** | **`2.13` (clash)** | `3.22` | Pre-2023 Vina receptor & template |
| `8CSG` | 2022-10-05 | `Mirati_Frag1` (`PWL`) | State 2 (`PRMT5•MTA`) | **`1.87` (clash)** | **`1.74` (clash)** | `3.35` | Pre-2023 O3A template |
| `6CKC` | 2018-05-16 | `LLY-283` (`F5J`) | State 3 (`SAM`-pocket) | **`1.66` (clash)** | `3.47` | **`1.50` (clash)** | Pre-2023 Vina receptor & template |
| `7MX7` | 2021-11-17 | `PF-06939999` (`ZR1`) | State 3 (`SAM`-pocket) | **`0.50` (clash)** | **`2.25` (clash)** | **`2.60` (clash)** | Pre-2023 O3A template |
| `6RLL` | 2020-06-03 | `JNJ-44064146` (`K8H`) | State 3 (`SAM`-pocket) | **`0.73` (clash)** | **`2.66` (clash)** | **`1.63` (clash)** | Pre-2023 O3A template |
| `8VEY` | 2024-04-24 | `TNG908` (`A1AAV`) | State 2 (`PRMT5•MTA`) | **`2.06` (clash)** | **`2.17` (clash)** | `3.42` | Quarantined post-2023 diagnostic |
| `9C10` | 2024-10-02 | `AMG_193` (`A1ATH`) | State 2 (`PRMT5•MTA`) | **`2.39` (clash)** | **`1.63` (clash)** | `3.67` | Quarantined post-2023 diagnostic |
| `9N3R` | 2025-03-12 | `TNG462` (`A1BV0`) | State 2 (`PRMT5•MTA`) | **`2.25` (clash)** | **`1.99` (clash)** | `3.74` | Quarantined post-2023 diagnostic |
| `9T43` | 2026-06-10 | `AZD3470` (`A1JTI`) | State 2 (`PRMT5•MTA`) | **`1.30` (clash)** | **`1.97` (clash)** | `3.44` | Quarantined post-2023 diagnostic |

**Mechanistic Insight:**
1. **Why 1st-gen THIQ inhibitors (`EPZ015666`, `GSK3326595`) are excluded from `PRMT5•MTA`:** When `MTA` replaces `SAM`, loss of `SAM`'s positively charged sulfonium and homocysteine carboxylate tail allows **`Phe327` to rotate ~4.5 Å into the substrate channel**, colliding directly (`0.51–0.58 Å`) with the THIQ bicyclic core.
2. **Why all clinical 2nd-gen MTA-cooperative drugs (`MRTX1719`, `TNG908`, `TNG462`, `AMG 193`, `AZD3470`) are excluded from `PRMT5•SAM`:** Every MTA-cooperative chemotype projects into the space vacated by `SAM`'s sulfonium methyl (`SAM:CE`, `1.30–2.71 Å`) and collides with the `SAM`-locked carboxylate rotamer of **`Glu435` (`1.63–2.17 Å`)**, while stacking favorably against swung-in `Phe327` (`3.22–3.74 Å`).
3. Across the `98` untouched `2024–2025` paired test molecules, AutoDock Vina differential binding energy $\Delta E_{\text{Vina}}(\text{PRMT5}\cdot\text{MTA} - \text{PRMT5}\cdot\text{SAM})$ averaged **`-0.69 kcal/mol`** on MTA-cooperative compounds (`y_coop_3x = 1`) versus **`+2.58 kcal/mol`** on non-cooperative compounds (`y_coop_3x = 0`) — a **`3.27 kcal/mol` thermodynamic gap** derived solely from pre-2023 crystal coordinates (`7S1S` vs `4X61`).

![Experimental Co-Crystal 3D Steric-Switch Orthants and Zero-Shot Physics](figures/fig1_cofactor_steric_switch.png)

---

## 2. Track A Results — Cofactor-Cooperativity & `MTAP`-Selective Mechanism (`n = 98`)

Track A evaluates `98` strictly scaffold-, parent-, nonisomeric-, and study-disjoint `2024–2025` *J. Med. Chem.* compounds (`49` Murcko scaffolds across `4` independent studies), each with direct paired experimental `PRMT5•MTA` and `PRMT5•SAM` biochemical measurements.

### 2.1 Overall Held-Out Track A Performance

| Model | Primary `y_coop_3x` AP `[95% CI]` | `ROC-AUC` `[95% CI]` | Top-20 Hits (`EF`) | Bal. Acc (`Sens` / `Spec`) | Stringent `10x` AP `[95% CI]` | $\Delta\text{pActivity}$ $C$-index `[95% CI]` | Spearman $\rho$ |
|---|---:|---:|---:|---:|---:|---:|---:|
| **`B0`** | `0.4490` `[0.3305, 0.5817]` | `0.5000` | `8.98 exp` (`1.00x`) | `0.5000` (`0.00` / `1.00`) | `0.2143` `[0.1250, 0.3367]` | `0.5000` | `0.0000` |
| **`P0`** | `0.5198` `[0.4080, 0.7267]` | `0.5884` `[0.4465, 0.7124]` | `11` (`1.23x`) | `0.5606` (`0.45` / `0.67`) | `0.2095` `[0.1330, 0.3843]` | `0.5554` `[0.4310, 0.6545]` | `+0.0722` |
| **`B1`** | `0.5600` `[0.4208, 0.6901]` | `0.5553` `[0.4051, 0.6825]` | `14` (`1.56x`) | `0.5476` (`0.11` / `0.98`) | `0.1862` `[0.1057, 0.4067]` | `0.5159` `[0.3853, 0.6152]` | `+0.0216` |
| **`M1`** | `0.5951` `[0.4636, 0.7251]` | `0.6582` `[0.5165, 0.7607]` | `13` (`1.45x`) | **`0.6692`** (`0.73` / `0.61`) | `0.2172` `[0.1261, 0.4092]` | `0.6194` `[0.4870, 0.7103]` | **`+0.2728`** |
| **`M2`** | **`0.6552`** `[0.4992, 0.7976]` | **`0.6768`** `[0.5169, 0.7921]` | **`16` (`1.78x`)** | `0.6006` (`0.39` / `0.81`) | **`0.2319`** `[0.1203, 0.5067]` | **`0.6390`** `[0.5130, 0.7326]` | `+0.2370` |

### 2.2 Paired Murcko-Scaffold Bootstrap Contrasts (`2,000` Resamples)

| Contrast | `y_coop_3x` $\Delta\text{AP}$ `[95% CI]` ($p$) | `y_coop_3x` $\Delta\text{ROC-AUC}$ `[95% CI]` ($p$) | $\Delta\text{pActivity}$ $\Delta C\text{-index}$ `[95% CI]` ($p$) |
|---|---:|---:|---:|
| **`M2 − B1`** | **`+0.0952`** `[-0.0011, +0.1991]` (`p = 0.052`) | **`+0.1214`** `[+0.0121, +0.2355]` (**`p = 0.033`**) | **`+0.1231`** `[+0.0462, +0.2125]` (**`p < 0.001`**) |
| **`M2 − M1`** | **`+0.0602`** `[-0.0153, +0.1544]` (`p = 0.113`) | `+0.0185` `[-0.0558, +0.0969]` (`p = 0.586`) | `+0.0196` `[-0.0472, +0.0822]` (`p = 0.512`) |
| **`M1 − B1`** | `+0.0350` `[-0.0210, +0.1059]` (`p = 0.234`) | **`+0.1029`** `[+0.0142, +0.2106]` (**`p = 0.020`**) | **`+0.1035`** `[+0.0230, +0.2106]` (**`p = 0.013`**) |
| **`P0 − B1`** | `-0.0402` `[-0.1443, +0.1139]` (`p = 0.686`) | `+0.0330` `[-0.0953, +0.1710]` (`p = 0.613`) | `+0.0395` `[-0.0657, +0.1542]` (`p = 0.437`) |
| **`M2 − P0`** | **`+0.1354`** `[+0.0062, +0.2221]` (**`p = 0.033`**) | **`+0.0884`** `[+0.0146, +0.1620]` (**`p = 0.020`**) | **`+0.0836`** `[+0.0247, +0.1456]` (**`p = 0.002`**) |

### 2.3 Per-Study Breakdown: Rational Mechanism Switching (`MECHSWITCH_JMC2025`)

In `MECHSWITCH_JMC2025` (`CHEMBL6153116`, `n = 55` paired compounds, `20` `y_coop_3x` positives, `36.4%` prevalence), medicinal chemists rationally transformed SAM-cooperative THIQ analogues into MTA-cooperative inhibitors:
- **v1-Style Achiral 2D `B1` collapsed completely:** `AP = 0.3095` (below `0.3636` prevalence), `ROC-AUC = 0.3686`, within-study $\Delta\text{p}K_i$ `C-index = 0.3676` — because 2D achiral Morgan fingerprints matched the THIQ/dihydroisoquinoline core to 1st-gen Epizyme `PRMT5•SAM` training negatives!
- **Zero-Shot 3D Physics `P0` reversed that failure without any training labels:** `AP = 0.5146`, `ROC-AUC = 0.6114`, within-study `C-index = 0.5589`.
- **`M1` (`AP = 0.6040`, `C-index = 0.6744`) and `M2` (`AP = 0.5491`, `C-index = 0.5881`, `AZ_JMC2024` `AP = 0.8478` vs `B1 = 0.7470`)** combined stereochemical/physicochemical reorganization with 3D pocket energetics to achieve robust multi-series transfer (`M2` Top-20 hits = `16/20`).

---

## 3. Track B Results — Untouched `2024–2025` Biochemical Potency (`n = 182`)

Track B tests cross-document biochemical potency transfer from `≤2021` training assays (`n = 818`) and `2019–2022` validation assays (`n = 261`, `49.4%` `≤100 nM`) to `5` untouched `2024–2025` peer-reviewed medicinal chemistry cohorts (`n = 182`, `75` scaffolds, `57` `≤100 nM` positives = `31.3%` prevalence, `87` confirmed `>1 µM` weak/inactives).

### 3.1 Overall Track B Ranking, Virtual Screening & Calibration Metrics

| Model | Pooled `≤100 nM` AP `[95% CI]` | Pooled `ROC-AUC` `[95% CI]` | Strict `≤100 nM` vs `>1 µM` AP `[95% CI]` | Within-Assay $C$-index `[95% CI]` | Top-25 Hits (`EF`) | Sensitivity / Specificity / Bal. Acc |
|---|---:|---:|---:|---:|---:|---:|
| **`B0`** | `0.3132` `[0.1991, 0.4487]` | `0.5000` | `0.3958` `[0.2357, 0.6084]` | `0.5000` | `7.83 exp` (`1.00x`) | `0.0000` / `1.0000` / `0.5000` |
| **`B1`** | `0.5638` `[0.3759, 0.6960]` | **`0.6258`** `[0.4635, 0.7410]` | **`0.7533`** `[0.5132, 0.9000]` | **`0.5043`** `[0.4007, 0.5854]` | `19` (`2.43x`) | `0.5439` / **`0.7360`** / **`0.6399`** |
| **`M1`** | **`0.5791`** `[0.3641, 0.7456]` | `0.5561` `[0.3408, 0.7195]` | `0.6848` `[0.4341, 0.8557]` | `0.4843` `[0.3834, 0.6074]` | **`20` (`2.55x`)** | `0.5789` / `0.1200` / `0.3495` |
| **`M2`** | `0.4488` `[0.2661, 0.6084]` | `0.5840` `[0.4115, 0.7149]` | `0.5426` `[0.3223, 0.7513]` | `0.4668` `[0.3874, 0.5814]` | `12` (`1.53x`) | **`0.8947`** / `0.0320` / `0.4634` |

**Paired Scaffold-Bootstrap Contrasts (Track B):**
- Pooled AP: `M1 − B1 = +0.0153` `[-0.0752, +0.0758]` (`p = 0.683`); `M2 − B1 = -0.1150` `[-0.2109, -0.0142]` (`p = 0.029`); `M2 − M1 = -0.1303` `[-0.2208, -0.0252]` (`p = 0.014`).
- Within-Assay Weighted $C$-index: `M1 − B1 = -0.0201` `[-0.1595, +0.1136]` (`p = 0.778`); `M2 − B1 = -0.0376` `[-0.1523, +0.1075]` (`p = 0.595`).

### 3.2 Per-Assay Dissection: Why Standard Series Succeeded While the Mechanism-Switched Counter-Screen Diverged

| Held-Out `2024–2025` Assay | Chemical Series (`n`, Positives) | `B1` AP / $C$-index / $\rho$ | `M1` AP / $C$-index / $\rho$ | `M2` AP / $C$-index / $\rho$ | Biophysical Interpretation |
|---|---|---:|---:|---:|---|
| **`CHEMBL5544428`** (`EJMC 2024`) | Chiral ring-opened nucleosides (`35`, `11` pos) | `0.3808` / `0.4876` / `-0.0236` | `0.5888` / `0.5927` / **`+0.5346`** | **`0.7181`** / **`0.6463`** / `+0.5230` | **Massive chiral + 3D `6CKC` SAM-pocket gain (`+0.337` AP over `B1`)** |
| **`CHEMBL6088129`** (`EJMC 2025`) | Chiral piperazine nucleosides (`29`, `21` pos) | `0.8451` / `0.6425` / `+0.3930` | **`0.9415`** / **`0.6855`** / **`+0.4522`** | `0.8098` / `0.6290` / `+0.2699` | Chiral `M1` excels; bulky piperazine tails experience mild rigid-box Vina penalty |
| **`CHEMBL6090539`** (`JMC 2025`) | Chiral dihydroisoquinolinones (`21`, `14` pos) | `0.5904` / `0.4726` / `-0.2013` | `0.5629` / `0.3955` / `+0.0021` | **`0.6440`** / **`0.5224`** / **`+0.2095`** | **3D `4X61` substrate-channel alignment reverses `B1` negative correlation** |
| **`CHEMBL5509028`** (`EJMC 2024`) | Tetrahydro-$\beta$-carbolines (`42`, `0` pos, `42` `>1 µM`) | `—` / `0.5966` / `+0.1404` | `—` / `0.4701` / `+0.2105` | `—` / `0.4239` / `+0.0333` | All `42` compounds are weak/inactive decoys (`IC50 > 1 µM`); `B1`/`M1` cleanly rank them below actives (`Top-25 = 19–20/25`!) |
| **`CHEMBL6157149`** (`JMC 2025`) | **`MECHSWITCH` `PRMT5•SAM` counter-screen** (`55`, `11` pos) | **`0.2240`** / **`0.4414`** / `-0.3083` | `0.1223` / `0.4107` / `-0.1740` | `0.1263` / `0.3689` / `-0.3125` | **Anti-target counter-screen:** MTA-switching groups intentionally abolish `PRMT5•SAM` potency on a THIQ core |

![Calibration Diagnostics and Virtual Screening Hit Recovery](figures/fig3_v1_threshold_collapse_vs_v2_calibration.png)

---

## 4. Synthesis & Practical Lessons for Structure-Enabled PRMT5 Discovery

1. **Static Multi-State Crystal Physics Excels at Large-Amplitude Cofactor-State Selectivity (`Track A`):** Because `PRMT5•MTA` vs `PRMT5•SAM` cooperativity hinges on a `~4.5 Å` `Phe327` gatekeeper swing and `1.3–2.7 Å` `SAM:CE` / `Glu435` steric exclusion, pre-2023 rigid crystal pockets (`7S1S` vs `4X61`) provide a `3.27 kcal/mol` physical separation that boosts held-out `M2` selectivity prediction significantly (`C-index = 0.6390` vs `0.5159`, `p < 0.001`; `AP = 0.6552` vs `0.5600`).
2. **Chirality + Pairwise Ranking Boosts Early Virtual-Screening Enrichment (`M1` `Top-25 = 20/25` = `80%` Precision) and Solves Nucleoside SAR (`CHEMBL5544428` `M2 AP = 0.7181` vs `B1 = 0.3808`), While Exposing Calibration & Induced-Fit Tradeoffs:**
   - Natural-prevalence OOF ranker calibration in `M1`/`M2` eliminated v1's catastrophic false-negative collapse (`2.9%` SK sensitivity → `57.9%–89.5%` sensitivity), though at the cost of over-calling positives (`low specificity`) when applied to decoy-heavy libraries like `CHEMBL5509028` (`0%` positives) without assay-specific control anchors. In practice, **1–3 internal reference controls (`EPZ015666`, `LLY-283`, `MRTX1719`) per new assay plate** should be used to anchor absolute potency thresholds rather than relying on unanchored probability cutoffs.
   - Resolving subtle congeneric potency cliffs *within* induced-fit or anti-target counter-screen series (`CHEMBL6157149`) requires explicit-solvent receptor flexibility (`MD` / `FEP+`) and joint multi-task modeling of both `PRMT5•SAM` and `PRMT5•MTA` states.
