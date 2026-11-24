# Model Card — PRMT5 v2 Multi-State Crystal Pocket Physics & Pairwise Assay Rankers

## 1. Model Overview

PRMT5 v2 evaluates deterministic, CPU-reproducible structural chemistry and statistical learning tiers across two orthogonal drug-discovery tasks:
1. **Track A (`PRMT5•MTA` vs `PRMT5•SAM` Cofactor Cooperativity & Selectivity):** Predicting whether an unseen small molecule selectively inhibits the `MTAP`-deleted tumor target state (`PRMT5•MTA`) over the healthy-tissue state (`PRMT5•SAM`).
2. **Track B (Cross-Document PRMT5 Biochemical Potency Ranking & Calibration):** Ranking potent (`≤100 nM`) inhibitors and maintaining calibrated operating thresholds across shifted `[SAM]` concentrations (`75 nM` to `10 µM`) and assay formats on untouched `2024–2025` medicinal chemistry cohorts.

All models run deterministically on CPU (`≤8` threads, `seed = 20261120`) without external network access or proprietary checkpoints.

---

## 2. Representation & Feature Architecture

1. **Achiral 2D Baseline (`B1`):** `2,048`-bit Morgan fingerprint (`radius = 2`, `includeChirality = False`) computed on canonical stereochemistry-stripped SMILES (`nonisomeric_smiles`), replicating v1's non-chiral input bottleneck.
2. **Chiral 2D + Physicochemical Representation (`M1`):** `2,048`-bit chiral Morgan fingerprint (`radius = 2`, `includeChirality = True`) + `16` physicochemical and stereochemical descriptors (`MolWt`, `MolLogP`, `TPSA`, `NumHDonors`, `NumHAcceptors`, `NumRotatableBonds`, `FractionCSP3`, `NumAromaticRings`, `NumAliphaticRings`, `RingCount`, `HeavyAtomCount`, `FormalCharge`, `QED`, `NumChiralCenters`, `HasSConfig`, `HasRConfig`) + training-only chiral Tanimoto similarity differential (`max_sim_pos − max_sim_neg`).
3. **3D Multi-State Crystal Pocket Physics (`P0` and `M2`):**
   - **Conformer Ensemble:** `RDKit ETKDGv3` (`numConfs = 8`, `useSmallRingTorsions = True`, `randomSeed = 20261120`) minimized with `MMFF94s` (`maxIters = 200`), recording internal strain energy (`conformer_strain_s1/s2_kcal`) and 3D shape moments (`radius_of_gyration_3d`, `npr1_3d`, `npr2_3d`).
   - **Multi-Template Pharmacophore & Shape Superposition:** `RDKit CrippenO3A` (`o3a.Align()` + `o3a.Score()`) and `ShapeTanimotoDist` across pre-2023 experimental crystal ligands in the unified `4X61` catalytic frame:
     - *State 1 (`PRMT5•SAM` substrate channel):* `4X61:3XV` (`EPZ015666`), `5EML:5QK`
     - *State 2 (`PRMT5•MTA` cooperative pocket):* `7S1S:85K` (`MRTX1719`), `8CSG:PWL` (`Fragment 1`)
     - *State 3 (`PRMT5` SAM-pocket nucleoside):* `6CKC:F5J` (`LLY-283`), `7MX7:ZR1` (`PF-06939999`), `6RLL:K8H` (`JNJ-44064146`)
   - **AutoDock Vina 1.2.7 BFGS Pocket Optimization:** Local torsion + rigid-body energy minimization (`v.optimize`) in `4X61` (`PRMT5•SAM`), `7S1S` (`PRMT5•MTA`), and `6CKC` (`PRMT5` SAM-pocket), yielding state-specific binding energies, $\Delta E_{\text{Vina}}(\text{MTA} - \text{SAM})$, and cross-pose incompatibility energies ($E_{\text{SAM}|\text{MTA-pose}} - E_{\text{MTA}|\text{SAM-pose}}$).
   - **Explicit Cofactor & Gatekeeper Steric Repulsion:** Soft-core Lennard-Jones repulsion and interatomic contact counts against `4X61`'s `SAM` sulfonium/homocysteine tail (`N, CA, C, O, OXT, CB, CG, SD, CE`), SAM-state `Glu435`, and `7S1S`'s swung-in MTA-state `Phe327` aromatic ring (`steric_switch_differential`).

---

## 3. Model Tiers & Validation-Only Selection

All parameters, scalers, and out-of-fold meta-weights (`5-fold` Murcko-scaffold `GroupKFold` strictly inside `train`) are fitted exclusively on `train` rows. Candidate selection and decision thresholds are locked on `validation` in [`manifests/model_selection_lock.json`](../manifests/model_selection_lock.json) prior to test evaluation:

- **`B0` (Constant Reference):** Predicts training positive prevalence.
- **`P0` (Track A Only — Zero-Shot 3D Cofactor-State Pocket Physics):** Unsupervised physical chemistry composite combining standardized $-\Delta E_{\text{Vina}}(\text{MTA}-\text{SAM})$, cross-pose incompatibility, $\Delta S_{\text{O3A}}(\text{MTA}-\text{SAM})$, $\Delta S_{\text{shape}}(\text{MTA}-\text{SAM})$, and `steric_switch_differential`. **Uses zero bioactivity labels for ranking.**
- **`B1` (v1-Style Achiral 2D Classifier):** Scaled L2 Logistic Regression (`C ∈ {0.1, 1.0, 10.0}`, balanced class weights) on non-chiral Morgan fingerprints; single probability threshold locked by validation balanced accuracy.
- **`M1` (Chiral 2D + Assay-Conditioned Ranker/Ensemble):**
  - *Track A:* Chiral 2D + physicochemical ensemble (L2 Logistic, regularized `XGBClassifier`, and chiral Tanimoto differential; candidate selected by validation AP).
  - *Track B:* Within-assay pairwise potency ranker (`XGBRanker(objective='rank:pairwise')` grouped strictly by `assay_id` + within-assay demeaned `Ridge` + assay-covariate logistic) calibrated via `train` scaffold-OOF logistic calibration on `(rank_score, log10_sam_uM, is_displacement)` anchored to base rate, preventing validation-prevalence threshold collapse.
- **`M2` (Two-Branch Chiral 2D + 3D Multi-State Pocket Physics Hybrid):** Fuses the `M1` chiral 2D branch with a dedicated 3D crystal pocket physics branch (regularized linear/tree models on `phys3d` + a priori 3D physics composite) using `train` Murcko-scaffold OOF stacking, preventing 2,048-bit fingerprint dimensionality from diluting the 34 physical energy terms.

---

## 4. Intended Use & Out-of-Scope Limitations

- **Intended Use:** Methodological benchmarking of 3D cofactor-state steric physics, stereochemistry preservation, and within-assay pairwise ranking / `[SAM]`-conditioned calibration on retrospective PRMT5 medicinal chemistry series.
- **Out-of-Scope:** Clinical dosing, *in vivo* toxicology prediction, covalent PRMT5 Cys449 warhead reactivity modeling, or allosteric non-catalytic adaptor PPI (`pICln`/`RioK1` PBM-groove) inhibitor design without dedicated PBM-site receptor grids.
