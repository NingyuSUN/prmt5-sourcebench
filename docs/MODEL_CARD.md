# Model card — PRMT5 SourceBench

## Intended use and outcome

Research comparison of structural models for ranking source-defined strong PRMT5 biochemical inhibitors in an independent patent series. The label means source-defined strong biochemical inhibition around 100 nM: early A bins are <100 nM and later sources use ≤100 nM under their respective protocols. It does not imply a common potency scale, calibrated screening probability, cellular efficacy, selectivity, binding mechanism or clinical benefit. See the [final report](../reports/RESULTS.md) for measured performance and limitations.

## Data and representation

604 training parents (164 positive), 70 validation parents (39 positive), 246 test parents (69 positive). Parent, scaffold, patent family and nonisomeric identity do not overlap across splits. Training-only fitting and scaling; validation-only configuration/epoch/threshold selection; no validation refit. Full original intervals, exclusions and disclosure-date bounds remain in the [protocol](PROTOCOL.md).

Primary inputs discard stereochemical markers to match the documented pretrained representation. All 920 inputs tokenize without unknown symbols or truncation (29–71 tokens). Within-split nonisomeric collisions include conflicting labels: the empirical minimum classification errors from this representation alone are 20/604 training, 6/70 validation and 4/246 test. These are representation bounds, not observed model errors. A prespecified chiral-fingerprint logistic sensitivity model retains stereochemical information.

## Models actually fitted

Ten fingerprint/reference configurations were fitted, then repeated once to check deterministic CPU prediction reproducibility. Repetitions are not independent scientific experiments. The frozen encoder head used three training-only scaling/logistic configurations. Six partial-adaptation trials completed, totaling 39 epochs across trials. All alternative validation runs and checkpoints are preserved.

Base encoder: `ibm-research/MoLFormer-XL-both-10pct`, revision `361063d0ad524ef77cf39b08469f6be770dc550f`; 44,375,040 loaded encoder parameters. Fine-tuning updated **7,093,248 encoder parameters plus 769 classifier parameters**; 37,281,792 encoder parameters stayed frozen. The training scaler and 12 pretrained random-feature buffers stayed fixed. The encoder remained in deterministic evaluation mode while gradients updated the last two blocks.

The frozen head selected C = 1. Partial adaptation selected encoder learning rate 1e-05 by mean validation AP across seeds 17, 29 and 43. Final adapted predictions average those three probabilities. The primary classical comparator is the validation-selected random forest (500 trees, leaf size 5). Scores are uncalibrated.

| Encoder learning rate | Seed | Epochs run | Selected epoch | Validation AP |
|---:|---:|---:|---:|---:|
| 1e-05 | 17 | 7 | 3 | 0.740025 |
| 1e-05 | 29 | 7 | 3 | 0.734477 |
| 1e-05 | 43 | 7 | 3 | 0.737400 |
| 3 e-05 | 17 | 5 | 1 | 0.738503 |
| 3 e-05 | 29 | 7 | 3 | 0.698651 |
| 3 e-05 | 43 | 6 | 2 | 0.740437 |

Frozen-head validation AP: 0.733859; adapted-ensemble validation AP: 0.740232; selected random-forest validation AP: 0.683772. These are model-selection results on reused validation data, not independent-test evidence.

## Artifacts and reproducibility

- [Complete encoder development record](../reports/models/encoder_development.json), [baseline record](../reports/models/validation_baselines.json).
- [Input/runtime preflight](../reports/models/encoder_preflight.json) and [pinned assets](../manifests/pretrained_model.json).
- [Model selection lock](../manifests/model_selection_lock.json) and [final scoring lock](../manifests/final_selection_lock.json).
- Full baseline models, all six best/last checkpoints and all validation predictions are stored under ignored `data/models/`. Adapted checkpoints contain parameter deltas and require the original pinned base weights. They safely load with `weights_only=True`; baseline joblib files require their recorded sklearn runtime.
- [Reproduction commands and delivery scope](REPRODUCIBILITY.md). Exact saved-artifact replay is distinct from fresh-environment retraining.

## Known limitations

Source, assay conditions, chemistry and publication date change together. Only one patent family is held out. Scaffold-bootstrap intervals quantify uncertainty inside this family, not variation across future sources. First-ever disclosure for every structure–label pair and original-figure reconstruction are incomplete. The modern encoder cannot establish what a historically available model would predict. Unlabelled pretraining structure overlap is unknown; no PRMT5-label training was identified in the inspected publisher documentation. No wet experiments, prospective screening, ADMET evaluation or external cellular validation were performed.

Project code and documentation use the [MIT license](../LICENSE); external data and model licenses remain separate. ChEMBL-derived data, original patent sources and MoLFormer assets must retain their own attribution and applicable terms.
