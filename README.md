# PRMT5 SourceBench

**A source-audited benchmark for molecular model transfer across PRMT5 inhibitor generations and cofactor states.**

[中文说明](README.zh-CN.md) · [Catalytic transfer report](reports/RESULTS.md) · [催化口袋中文报告](reports/RESULTS.zh-CN.md) · [Cofactor & 3D report](v2/reports/RESULTS.md) · [三维物理与协同选择性中文报告](v2/reports/RESULTS.zh-CN.md) · [Data cards](docs/DATA_CARD.md) · [3D data card](v2/docs/DATA_CARD.md) · [Model cards](docs/MODEL_CARD.md) · [3D model card](v2/docs/MODEL_CARD.md) · [Reproduce](docs/REPRODUCIBILITY.md)

[![Validate benchmark](https://github.com/NingyuSUN/prmt5-sourcebench/actions/workflows/validate.yml/badge.svg)](https://github.com/NingyuSUN/prmt5-sourcebench/actions/workflows/validate.yml)

PRMT5 SourceBench evaluates when molecular representations transfer across chemically distinct patent and journal series targeting human PRMT5/MEP50. Across **1,514 source-audited unique parent structures** (curated from **17,697** `ChEMBL_37` records and original patent/journal tables) and **13 human PRMT5 crystal structures**, every evaluation split enforces zero overlap in isomeric parent identity, nonisomeric connectivity, Bemis–Murcko scaffold, and source family.

**Core findings:**
1. **First-generation catalytic patent transfer (`2017–2021`, `920` parents) exposes two coupled failure modes of achiral 2D adaptation.** Partially fine-tuned `MoLFormer` reached held-out average precision **`0.2553`** versus **`0.2805`** for constant prevalence (`ROC-AUC = 0.4518`, **`2/25`** top hits vs `7.01` expected at random), recalling only **`2/69`** positives at its validation threshold when SAM concentration shifted `133×` (`75 nM → 10 µM → 1 µM`) into a congeneric tetrahydroisoquinoline (`THIQ`) optimization series (`69` `≤100 nM` vs `174` `100–1000 nM`).
2. **Crystallographic cofactor-state physics (`PRMT5•MTA` vs `PRMT5•SAM`) resolves second-generation `MTAP`-deleted selectivity.** Superimposing `13` human PRMT5 co-crystals onto a common `331`-atom catalytic-domain $C_\alpha$ frame (`PDB 4X61`, RMSD `0.41–0.77 Å`) isolates a reciprocal two-gate steric switch: without SAM's homocysteine carboxylate tail in `PRMT5•MTA`, `Phe327` swings `~4.5 Å` into the substrate channel (`0.51–0.58 Å` clash with first-generation `THIQ` ligands), whereas MTA-cooperative inhibitors (`MRTX1719`, `TNG908`, `TNG462`, `AMG 193`, `AZD3470`) stack against `Phe327` (`3.22–3.74 Å`) and collide (`1.30–2.71 Å`) with `SAM:CE` and SAM-locked `Glu435` (`1.63–2.17 Å`). Pre-2023 `7S1S` versus `4X61` `AutoDock Vina` scoring separates `98` untouched `2024–2025` paired compounds by **`3.27 kcal/mol`** (`−0.69` vs `+2.58 kcal/mol`) zero-shot without bioactivity labels.
3. **Chiral 2D + 3D multi-state pocket physics transfers to `2024–2025` medicinal chemistry while delineating static-docking boundaries.** On **`98` paired `2024–2025` *J. Med. Chem.* `PRMT5•MTA` vs `PRMT5•SAM` compounds**, chiral 2D + 3D pocket physics reaches **`AP = 0.6552`** (`16/20` Top-20 hits vs `0.4490` reference), outperforming achiral 2D (`AP = 0.5600`) by **`ΔROC-AUC = +0.1214 [+0.0121, +0.2355]` (`p = 0.033`)** and continuous selectivity concordance **`ΔC-index = +0.1231 [+0.0462, +0.2125]` (`p < 0.001`)**. Across **`182` untouched `2024–2025` biochemical compounds**, chiral within-assay pairwise ranking reaches **`AP = 0.5791`** (`20/25` Top-25 hits = `2.55×` enrichment vs `0.3132` reference); adding 3D pocket features lifts chiral ring-opened nucleoside (`CHEMBL5544428`) `AP` from `0.3808` to **`0.7181`** (`Spearman ρ = +0.5230`), while transparently documenting anti-target counter-screen behavior (`CHEMBL6157149`).

![Crystallographic cofactor steric switch and zero-shot pocket physics](v2/reports/figures/fig1_cofactor_steric_switch.png)

## Study at a glance

| Component | Evidence |
|---|---|
| Public pharmacology & structural corpus | **1,514 unique parent structures** (`17,697` `ChEMBL_37` records + patent/journal tables) and **13 human PRMT5 PDB co-crystals** (`4X61` catalytic-domain $C_\alpha$ alignment; post-2023 crystals strictly quarantined from featurization) |
| Zero-leakage split contracts | Zero cross-split overlap in **isomeric parent, nonisomeric connectivity, Murcko scaffold, and source family** across all cohorts; the `246`-compound catalytic patent test set remains completely untouched by later cohorts ([quarantine manifest](v2/manifests/v1_sk_quarantine.json)) |
| Evaluation tracks | **1.** First-generation catalytic patent transfer (`920` parents) · **2.** Cofactor cooperativity `PRMT5•MTA` vs `PRMT5•SAM` (`616` parents, `98` paired test) · **3.** Post-2021 cross-series biochemical potency (`1,261` parents, `182` journal test) |
| Controlled model spectrum | Classical 2D fingerprints, frozen and partially fine-tuned `MoLFormer-XL-both-10pct`, zero-shot 3D crystal pocket physics, chiral 2D + within-assay pairwise `XGBRanker`, and chiral 2D + 3D multi-state pocket hybrids |
| Uncertainty & error audit | **2,000 paired Murcko-scaffold bootstrap resamples** per track, censor-aware concordance ($C$-index), `[SAM]` shift calibration diagnostics, stereochemical collision audits, and per-assay error decomposition |
| Automated verification | **90** catalytic replay/contract tests + **13** seal-chain, zero-leakage and tamper-rejection `pytest` checks; deterministic CPU-only scoring replay from tracked repository files |

## Scientific design and cohorts

PRMT5 inhibition spans three structurally distinct mechanisms: SAM-uncompetitive substrate-channel binders (`PRMT5•SAM` + `THIQ`), MTA-cooperative inhibitors that exploit elevated `MTA` in `MTAP`-deleted tumors (`PRMT5•MTA`), and SAM-pocket nucleoside mimetics. Because archival `ChEMBL` entries mix interval bins (such as `50 / 550 / 5500 nM` placeholders for `+++ / ++ / +` patent tables), `133×` differences in `[SAM]`, and anti-target counter-screens, each evaluation setting defines endpoints matched to source measurement fidelity:

| Evaluation setting | Train | Validation | Held-out test | Verified test window | Primary endpoint |
|---|---|---|---|---|---|
| **First-generation catalytic patent transfer** | `604` (`285` scaffolds) | `70` (`44` scaffolds) | **`246`** (`89` scaffolds, `1` patent family) | `2021-04-08` | Source-defined strong inhibition (`≤100 nM` / `<100 nM`, `28.0%` positive) |
| **Cofactor cooperativity (`PRMT5•MTA` vs `PRMT5•SAM`)** | `430` (`159` scaffolds) | `88` (`55` scaffolds) | **`98` paired** (`49` scaffolds, `4` *J. Med. Chem.* studies) | `2024–2025` | `≥3.16×` MTA selectivity ($\Delta\text{pActivity} \ge 0.5$, `44.9%` positive) & continuous $\Delta\text{pActivity}$ $C$-index |
| **Post-2021 cross-series biochemical potency** | `818` (`361` scaffolds) | `261` (`119` scaffolds) | **`182`** (`75` scaffolds, `5` journal assays) | `2024–2025` | Biochemical `IC50` / `Ki ≤ 100 nM` (`31.3%` positive) & within-assay $C$-index |

All configuration choices, stopping epochs, ensemble weights, and classification thresholds were locked on validation splits before scoring held-out test predictions ([catalytic lock](manifests/final_selection_lock.json), [3D protocol freeze](v2/receipts/protocol_freeze.json), [prediction seal](v2/receipts/prediction_seal.json)).

## Held-out performance across inhibitor generations

### 1. First-generation catalytic patent transfer (`n = 246`, `2021` patent series)

The molecular sequence encoder is [MoLFormer-XL-both-10pct](https://huggingface.co/ibm-research/MoLFormer-XL-both-10pct) (`361063d0ad524ef77cf39b08469f6be770dc550f`), comparing frozen linear probing and partial adaptation of the final two transformer blocks (`3` seeds averaged at the validation-selected learning rate) against classical 2D fingerprints.

![First-generation catalytic patent-transfer ranking and screening results](reports/figures/benchmark.png)

| Method | Test average precision | ROC-AUC | Positives among top 25 | Recalled positives at validation threshold |
|---|---:|---:|---:|---:|
| Constant prevalence reference | `0.2805` | `0.5000` | `7.01 expected` | `0 / 69` |
| Five-neighbor fingerprint | `0.2812` | `0.4955` | `7.00 expected` | `10 / 69` |
| Fingerprint logistic regression | `0.2704` | `0.4666` | `5` | `18 / 69` |
| Random forest (validation-selected comparator) | `0.2433` | `0.4337` | `3` | `10 / 69` |
| Frozen MoLFormer + linear head | `0.2551` | `0.4507` | `2` | `2 / 69` |
| Partially fine-tuned MoLFormer | `0.2553` | `0.4518` | `2` | `2 / 69` |
| Chiral fingerprint logistic sensitivity arm | `0.2749` | `0.4809` | `5` | `23 / 69` |

Partial fine-tuning minus random forest AP was **`+0.0120`** (paired 95% scaffold-bootstrap interval **`[−0.0337, +0.0824]`**), failing to exceed constant prevalence (`0.2805`). Structural audit explains why: `243/246` held-out compounds already retain the active `THIQ` substrate-pocket warhead (`69` `≤100 nM` vs `174` `100–1000 nM` on a benzoxazinone core), while training (`75 nM` SAM), validation (`10 µM` SAM, `55.7%` positive), and test (`1 µM` SAM, `28.0%` positive) span a `133×` cofactor shift under SAM-uncompetitive kinetics. [Full catalytic analysis →](reports/RESULTS.md)

### 2. Cofactor-state 3D pocket physics (`PRMT5•MTA` vs `PRMT5•SAM`, `n = 98` paired `2024–2025` compounds)

To test whether explicit pre-2023 active-site geometry resolves cofactor-dependent mechanism switching, multi-conformer `MMFF94s` ensembles and `AutoDock Vina` poses were scored strictly against pre-2023 crystal receptors (`4X61` `PRMT5•SAM`, `7S1S` `PRMT5•MTA`, `6CKC` nucleoside state) with zero post-2023 structural leakage.

![Held-out cofactor-cooperativity and post-2021 potency benchmarks](v2/reports/figures/fig2_track_a_b_benchmark_ci.png)

| Method | Primary `y_coop_3x` AP (95% CI) | ROC-AUC (95% CI) | Positives among top 20 (`EF`) | Stringent `10×` AP | $\Delta\text{pActivity}$ $C$-index (95% CI) |
|---|---:|---:|---:|---:|---:|
| Constant prevalence reference (`B0`) | `0.4490` `[0.3305, 0.5817]` | `0.5000` | `8.98 expected` (`1.00×`) | `0.2143` | `0.5000` |
| Zero-shot 3D crystal pocket physics (`P0`, no labels) | `0.5198` `[0.4080, 0.7267]` | `0.5884` `[0.4465, 0.7124]` | `11` (`1.23×`) | `0.2095` | `0.5554` `[0.4310, 0.6545]` |
| Non-chiral 2D Morgan logistic (`B1`) | `0.5600` `[0.4208, 0.6901]` | `0.5553` `[0.4051, 0.6825]` | `14` (`1.56×`) | `0.1862` | `0.5159` `[0.3853, 0.6152]` |
| Chiral 2D + physicochemical ensemble (`M1`) | `0.5951` `[0.4636, 0.7251]` | `0.6582` `[0.5165, 0.7607]` | `13` (`1.45×`) | `0.2172` | `0.6194` `[0.4870, 0.7103]` |
| **Chiral 2D + 3D multi-state pocket physics (`M2`)** | **`0.6552`** `[0.4992, 0.7976]` | **`0.6768`** `[0.5169, 0.7921]` | **`16` (`1.78×`)** | **`0.2319`** | **`0.6390`** `[0.5130, 0.7326]` |

Adding 3D multi-state crystal pocket features (`M2`) improves monotonically across every metric tier, outperforming non-chiral 2D (`B1`) by **`ΔAP = +0.0952 [−0.0011, +0.1991]` (`p = 0.052`)**, **`ΔROC-AUC = +0.1214 [+0.0121, +0.2355]` (`p = 0.033`)**, and continuous $\Delta\text{pActivity}$ concordance across `3,421` comparable pairs by **`ΔC-index = +0.1231 [+0.0462, +0.2125]` (`p < 0.001`)**. Pure zero-shot 3D pocket physics (`P0`, using no bioactivity labels) already exceeds supervised achiral 2D on continuous selectivity concordance (`0.5554` vs `0.5159`), and `M2` rescues the rational SAM-to-MTA mechanism-switching series (`CHEMBL6153116`, where `B1` falls to `AP = 0.3095` below `0.3636` prevalence, whereas `M2` reaches **`AP = 0.6094`**).

### 3. Post-2021 cross-series biochemical potency (`n = 182`, `2024–2025` peer-reviewed series)

| Method | Pooled `≤100 nM` AP (95% CI) | Strict `≤100 nM` vs `>1 µM` AP | Within-assay $C$-index | Positives among top 25 (`EF`) | Sensitivity at validation threshold |
|---|---:|---:|---:|---:|---:|
| Constant prevalence reference (`B0`) | `0.3132` `[0.1991, 0.4487]` | `0.3958` | `0.5000` | `7.83 expected` (`1.00×`) | `0.0000` |
| Non-chiral 2D Morgan logistic (`B1`) | `0.5638` `[0.3759, 0.6960]` | `0.7533` | `0.5043` | `19` (`2.43×`) | `0.5439` |
| **Chiral 2D + within-assay `XGBRanker` (`M1`)** | **`0.5791`** `[0.3641, 0.7456]` | `0.6848` | `0.4843` | **`20` (`2.55×`)** | `0.5789` |
| Chiral 2D + 3D multi-state physics + ranker (`M2`) | `0.4488` `[0.2661, 0.6084]` | `0.5426` | `0.4668` | `12` (`1.53×`) | **`0.8947`** |

Across heterogeneous `2024–2025` journal assays, chiral 2D pairwise ranking (`M1`) achieves the strongest pooled screening enrichment (**`20/25`** top hits, `2.55×` over random), while explicit 3D pocket features (`M2`) exhibit sharp assay-mechanism specificity:
- **Large gains on stereochemically gated nucleoside and lactam pockets:** On chiral ring-opened adenosine analogues (`CHEMBL5544428`, `n = 35`), `M2` raises `AP` from `0.3808` (`B1`, `Spearman ρ = −0.0236`) to **`0.7181`** (`Spearman ρ = +0.5230`). On chiral dihydroisoquinolinones (`CHEMBL6090539`, `n = 21`), `M2` lifts `AP` from `0.5904` (`ρ = −0.2013`) to **`0.6440`** (`ρ = +0.2095`).
- **Boundary on anti-target counter-screens:** Pooled `M2` enrichment drops because `CHEMBL6157149` (`n = 33`) is the `PRMT5•SAM` counter-screen of a rational SAM-to-MTA switching series in which `THIQ` analogues were modified specifically to collide with `SAM` and abolish `PRMT5•SAM` potency — an effect captured directly by the paired cofactor-cooperativity model above rather than single-assay potency pooling. [Full cofactor & 3D report →](v2/reports/RESULTS.md)

## Quick start: reproduce all published metrics on CPU

Use Python 3.12 in a clean virtual environment. Both saved-prediction replay pipelines run entirely on CPU using only tracked files in this repository; no GPU, pretrained checkpoints, or network credentials are required.

```bash
git clone https://github.com/NingyuSUN/prmt5-sourcebench.git
cd prmt5-sourcebench
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements/baseline.txt
python -m pip install pytest pyyaml pandas

# 1. Replay first-generation catalytic patent transfer (metrics + 2,000 scaffold resamples)
python scripts/reproduce_metrics.py --root . --output-dir outputs/replay_catalytic

# 2. Replay cofactor-state 3D pocket physics & 2024-2025 cross-series benchmarks
python v2/scripts/evaluate_locked.py --root v2 --output-dir outputs/replay_cofactor_3d

# 3. Run seal-chain, zero-leakage and tamper-rejection tests
pytest v2/tests -q
python -m unittest discover -s scripts -p 'test_*.py' -v
```

*(Running `python -m unittest discover -s tests -p 'test_*.py' -v` additionally checks the `MoLFormer` encoder training contracts and requires CPU PyTorch: `python -m pip install torch==2.6.0 --index-url https://download.pytorch.org/whl/cpu`.)*

Both replay commands verify SHA-256 input locks and refuse to overwrite non-empty output directories. They validate deterministic metric and bootstrap reproduction from sealed predictions rather than full model retraining from raw external archives. [Detailed reproduction scope →](docs/REPRODUCIBILITY.md)

## Repository guide

| Path | Contents |
|---|---|
| [`data/`](data/README.md) | Frozen `920`-parent catalytic patent cohort, exclusions, nearest neighbors, `1,722` model–compound predictions, bootstrap draws and error table |
| [`scripts/`](scripts/) | Catalytic cohort curation, `MoLFormer` partial adaptation, validation locks, scoring replay and biological error analysis |
| [`reports/RESULTS.md`](reports/RESULTS.md) · [`中文`](reports/RESULTS.zh-CN.md) | Detailed analysis of first-generation catalytic patent transfer and `[SAM]` threshold miscalibration |
| [`docs/`](docs/) · [`manifests/`](manifests/) · [`PROJECT_INDEX.md`](PROJECT_INDEX.md) | Catalytic data card, model card, protocol, provenance hashes and release evidence index |
| [`v2/configs/protocol_v2.yaml`](v2/configs/protocol_v2.yaml) | Preregistered cofactor-cooperativity and post-2021 cross-series protocol, split rules and pre-2023 PDB cutoffs |
| [`v2/data/`](v2/data/) · [`v2/manifests/`](v2/manifests/) · [`v2/receipts/`](v2/receipts/) | Standardized structures, `13`-crystal alignment metadata, sealed test labels/predictions, and SHA-256 protocol/prediction seals |
| [`v2/scripts/`](v2/scripts/) · [`v2/tests/`](v2/tests/) | 3D multi-state `MMFF94s`/`AutoDock Vina` pocket featurization, `XGBRanker`/hybrid training, locked evaluator and contract tests |
| [`v2/reports/RESULTS.md`](v2/reports/RESULTS.md) · [`中文`](v2/reports/RESULTS.zh-CN.md) | Complete crystallographic steric-switch audit, `2024–2025` held-out results, bootstrap contrasts and per-assay diagnostics |
| [`v2/docs/`](v2/docs/) | Cofactor & 3D pocket benchmark [protocol](v2/docs/PROTOCOL.md), [data card](v2/docs/DATA_CARD.md) and [model card](v2/docs/MODEL_CARD.md) |

## Interpretation and scope

This benchmark characterizes when 2D chemical representations suffice and when explicit cofactor-state 3D pocket geometry is required across historical and modern PRMT5 drug discovery:
- **Observational cross-source confounding:** Public patent and journal series change chemical scaffolds, `[SAM]` / `[MTA]` concentrations, substrate format, and publication dates together; retrospective holdouts measure transfer across real disclosure boundaries rather than isolated single-variable interventions.
- **Static crystal receptor approximation:** Pre-2023 `AutoDock Vina` and steric-clash features use rigid receptor coordinates (`4X61`, `7S1S`, `6CKC`) without molecular dynamics or induced-fit side-chain relaxation. While sufficient to capture the `~4.5 Å` `Phe327` swing and `SAM:CE`/`Glu435` steric exclusion, static scoring over-penalizes flexible pockets and anti-target counter-screens (`CHEMBL6157149`).
- **Pretraining and translational boundaries:** Exact pretraining-corpus membership of `MoLFormer` is unknown, and all results remain strictly computational; no prospective wet-lab hit discovery, cellular selectivity, ADMET, or clinical validation is claimed.

## Attribution and licensing

Code and original documentation: [MIT](LICENSE). Released `ChEMBL`-derived bioactivity tables and curated crosswalks: [CC BY-SA 3.0](data/LICENSE.md), with source document and patent identifiers preserved. Experimental macromolecular crystal structures are obtained from the [RCSB Protein Data Bank](https://www.rcsb.org/pages/usage-policy) under CC0 1.0 Public Domain ([PDB license notice](v2/.licenses/pdb_database_LICENSE.txt)). Third-party patent pages and pretrained model weights are not redistributed. Chemistry standardization builds on the pinned [DrugPerturb](https://github.com/NingyuSUN/drugperturb/tree/cfa256c201ef45ad54efdc78eb2cef529f211d51) pipeline. See [third-party notices](THIRD_PARTY_NOTICES.md) and [citation metadata](CITATION.cff).
