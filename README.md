# PRMT5 SourceBench

**A source-audited benchmark for molecular model transfer across PRMT5 patent series.**

[Results](reports/RESULTS.md) · [中文报告](reports/RESULTS.zh-CN.md) · [Data card](docs/DATA_CARD.md) · [Model card](docs/MODEL_CARD.md) · [Reproduce](docs/REPRODUCIBILITY.md)

PRMT5 SourceBench tests whether a pretrained molecular representation can rank strong biochemical inhibitors in a chemically distinct patent series. It connects original activity tables to standardized structures, separates patent families and scaffolds, compares classical models with frozen and partially fine-tuned MoLFormer, and releases the complete evaluation cohort and predictions.

**Main result:** partial fine-tuning did not establish useful screening enrichment in the held-out family. Its average precision was **0.255**, compared with **0.280** for the constant-prevalence reference. The top 25 predictions contained **2 positives**; random selection would yield **7.01 positives on average**. All comparison methods and uncertainty estimates are retained.

![Held-out ranking and screening results](reports/figures/benchmark.png)

[![Validate benchmark](https://github.com/NingyuSUN/prmt5-sourcebench/actions/workflows/validate.yml/badge.svg)](https://github.com/NingyuSUN/prmt5-sourcebench/actions/workflows/validate.yml)

## At a glance

| Component | Evidence |
|---|---|
| Public pharmacology data | **920 unique parent structures**, five patent families, original activity intervals and source-date evidence |
| Evaluation design | **604 train / 70 validation / 246 test**; no shared parent, nonisomeric identity, Murcko scaffold or patent family between splits |
| Model development | **10** reference/fingerprint configurations, **3** frozen-head configurations, **6** partial fine-tuning trials across two learning rates and three seeds |
| Controlled comparison | Same test compounds for seven methods; validation-only selection of configurations, epochs and operating thresholds |
| Uncertainty and interpretation | **2,000 paired scaffold resamples**, similarity strata, stereochemical collisions, source assay conditions and molecular descriptors |
| Verification | **90 contract tests passed**, including a relocated-checkout replay and altered-input rejection; [validation record](reports/public_release_validation.json) |
| Public reproduction | Frozen cohort, all **1,722 model–compound predictions**, original bootstrap draws and a CPU-only scoring command; no weights or credentials needed |

## Scientific question and dataset

The label is **source-defined strong biochemical inhibition at approximately 100 nM**: early letter-coded A bins are `<100 nM`, while later numeric/interval sources use `≤100 nM`. Original bin boundaries remain in the data. Values such as 50/550/5500 nM in ChEMBL can represent intervals; they are not precise IC50 measurements suitable for pooled regression.

| Split | Compounds | Positive / negative | Scaffolds | Verified publication bounds |
|---|---:|---:|---:|---|
| Training | 604 | 164 / 440 | 285 | 2017–2018 |
| Validation | 70 | 39 / 31 | 44 | 2019-09-12 |
| Final test | 246 | 69 / 177 | 89 | 2021-04-08 |

All included assays concern PRMT5/MEP50-mediated histone H4 methylation, but SAM concentration, assay readout and substrate conditions differ. Source, chemistry, date and experimental protocol change together. The benchmark measures retrospective **transfer across sources**, not an isolated time effect or a shared absolute-potency scale. Publication dates are verified public-by bounds, not exhaustive first-disclosure dates for every structure–label pair.

## Models and held-out performance

The molecular encoder is [MoLFormer-XL-both-10pct](https://huggingface.co/ibm-research/MoLFormer-XL-both-10pct), pinned to revision `361063d0ad524ef77cf39b08469f6be770dc550f`. Partial adaptation updates the final two transformer blocks and a linear classifier. Earlier layers, pretrained random-feature buffers and training-derived scaling remain fixed. Final adapted predictions average all three seeds at the validation-selected learning rate; no best-seed selection or validation refit is used.

| Method | Test average precision | ROC-AUC | Positives among top 25 |
|---|---:|---:|---:|
| Constant prevalence reference | 0.2805 | 0.5000 | 7.01 expected |
| Five-neighbor fingerprint | 0.2812 | 0.4955 | 7.00 expected |
| Fingerprint logistic regression | 0.2704 | 0.4666 | 5 |
| Random forest | 0.2433 | 0.4337 | 3 |
| Frozen MoLFormer + linear head | 0.2551 | 0.4507 | 2 |
| Partially fine-tuned MoLFormer | 0.2553 | 0.4518 | 2 |
| Chiral fingerprint logistic regression | 0.2749 | 0.4809 | 5 |

Average precision measures concentration of positives near the top of a ranking; it is not accuracy or a calibrated hit probability. The first two rows use expected hit counts because their ranking boundaries contain ties. The encoder counts are observed integers. Random forest was the strongest primary classical method on **validation**, so it remained the prespecified comparator even though other methods scored higher on the final test.

Adaptation minus random forest AP was **+0.0120**, with a paired 95% scaffold-bootstrap interval of **[−0.0337, +0.0824]**. This does not establish reliable improvement. Intervals describe uncertainty within one test family, not variation across future laboratories. [Full results, thresholds and error analysis →](reports/RESULTS.md)

## Quick start: reproduce the published metrics

Use Python 3.12 in a new virtual environment. This path runs entirely on CPU and uses only files included in the repository.

```bash
git clone https://github.com/NingyuSUN/prmt5-sourcebench.git
cd prmt5-sourcebench
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements/baseline.txt
python scripts/reproduce_metrics.py --root . --output-dir outputs/replay
```

The command verifies frozen hashes and compound coverage, recomputes ranking, screening and threshold metrics plus 2,000 paired scaffold resamples, and compares them with the published results. It refuses to overwrite prior output. **This is saved-prediction reproduction, not model retraining.** Full training code and recorded environments are included; source-page snapshots and trained checkpoints are not distributed. [Reproduction scope and commands →](docs/REPRODUCIBILITY.md)

## Repository guide

| Path | Contents |
|---|---|
| [`data/`](data/README.md) | Frozen cohort, exclusions, nearest-training neighbors, predictions, bootstrap draws and per-compound error table |
| [`scripts/`](scripts/) | Source curation, model training, selection locks, prediction, scoring and biological analysis |
| [`reports/RESULTS.md`](reports/RESULTS.md) | Biology-facing interpretation, comparison table and limitations |
| [`docs/`](docs/) | Data card, model card, protocol, implementation and reproduction details |
| [`manifests/`](manifests/) | Source identifiers, acquisition hashes, model revision and frozen experiment bindings |
| [`PROJECT_INDEX.md`](PROJECT_INDEX.md) | Canonical evidence and release index |

## Interpretation and scope

This study provides a reproducible account of where a fixed small-data adaptation strategy failed to transfer. Assay and chemical differences are plausible contributors, but their causal effects are not identified. Exact pretraining-corpus membership is unknown. There is one held-out patent family and no prospective screen, new wet experiment, cellular efficacy, selectivity, ADMET or clinical validation. No validated drug candidate or antibody discovery result is claimed.

## Attribution and licensing

Code and original documentation: [MIT](LICENSE). Released ChEMBL-derived data and derived data tables: [CC BY-SA 3.0](data/LICENSE.md), with source identifiers and patent references retained. Third-party patent pages and pretrained weights are not redistributed. Chemistry processing reuses the pinned [DrugPerturb](https://github.com/NingyuSUN/drugperturb/tree/cfa256c201ef45ad54efdc78eb2cef529f211d51) implementation. See [third-party notices](THIRD_PARTY_NOTICES.md) and [citation metadata](CITATION.cff).
