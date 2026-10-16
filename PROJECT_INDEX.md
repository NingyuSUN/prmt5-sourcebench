# PRMT5 SourceBench — evidence index

Updated: 2026-10-16, Asia/Tokyo. This is the canonical entry for the public benchmark.

- **Objective:** Document a completed, source-audited comparison of fingerprint models and pretrained molecular adaptation for PRMT5 biochemical activity classification.
- **Status:** Fixed experiment complete. Public data, code and saved-prediction replay passed 90 contract tests; the original metrics and all 2,000 bootstrap draws reproduced byte-for-byte.
- **Measured result:** 246 held-out compounds / 69 positives / 89 scaffolds / one family. Adapted AP 0.255312, frozen AP 0.255069, forest AP 0.243269, prevalence reference 0.280488. Adapted top-25 hits: 2. No useful screening enrichment or reliable improvement established.
- **Limits:** Source, chemistry, time and assay conditions are confounded; only one held-out family; exact first-disclosure and pretraining-corpus membership are incomplete. No prospective, cellular or clinical validation.
- **Further research:** Requires a new protocol and an untouched evaluation source; this final test is not a tuning set.

## Evidence map

| Topic | Artifact |
|---|---|
| Results | [English report](reports/RESULTS.md), [Chinese report](reports/RESULTS.zh-CN.md), [figure](reports/figures/benchmark.png) |
| Data and source provenance | [Data card](docs/DATA_CARD.md), [released tables](data/README.md), [family/date ledger](reports/curation/source_family_dates.csv), [cohort checks](reports/curation/cohort_gate.json) |
| Model methods and selection | [Model card](docs/MODEL_CARD.md), [scientific protocol](docs/PROTOCOL.md), [baseline runs](reports/models/validation_baselines.json), [encoder runs](reports/models/encoder_development.json) |
| Fixed test results | [Metrics](reports/results/final_metrics.json), [biological analysis](reports/results/biological_analysis.json), [prediction lock](manifests/final_selection_lock.json) |
| Reproduction | [Guide](docs/REPRODUCIBILITY.md), [public replay contract](docs/PUBLIC_SCORING_REPLAY.md), [frozen inputs](manifests/public_scoring_replay.json) |
| Verification | [Public replay validation](reports/public_release_validation.json), [original experiment checks](reports/FINAL_VALIDATION.md), [original scoring replay](reports/results/scoring_replay.json) |
| Attribution | [Code license](LICENSE), [data license](data/LICENSE.md), [third-party notices](THIRD_PARTY_NOTICES.md), [citation](CITATION.cff) |

## 2026-10-16 — Public benchmark release

| Field | Record |
|---|---|
| Date / timezone | 2026-10-16, Asia/Tokyo. Exact release commit and CI verification times are available from GitHub. |
| Purpose | Publish the completed study with unchanged scientific results and a CPU-only scoring replay. No new model fitting or final-test selection. |
| Inputs / code | Frozen 920-parent cohort, 1,722 prediction rows, original model/selection records, 2,000 bootstrap draws. Exact original code and artifact hashes: manifests/public_scoring_replay.json. Run scripts/reproduce_metrics.py from the release root. |
| Results / limitations | Same seven-model negative-transfer result and single-family uncertainty as the original experiment. Public replay verifies predictions/metrics, not model generation or biological efficacy. Source snapshots and weights are not distributed. |
| Outputs | Relative links above. Public protocol is an editorial copy; original and published hashes are recorded in manifests/protocol_publication.json. |
| Git / GitHub | Public release target: https://github.com/NingyuSUN/prmt5-sourcebench, branch main. Pre-publication state: locally prepared; remote synchronization must be verified after push. This row does not assert that ignored weights, raw snapshots or environment files are uploaded. |

The initial public release contains completed-study evidence. Earlier feasibility drafts and execution records are not part of this repository. The index's own version is available with `git log -1 -- PROJECT_INDEX.md`.
