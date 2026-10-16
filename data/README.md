# Released data and predictions

These are unchanged copies of frozen experiment artifacts. Scientific values and original CSV bytes are preserved. See [data licensing](LICENSE.md), [data card](../docs/DATA_CARD.md) and [source-family ledger](../reports/curation/source_family_dates.csv).

| File | Unit and scope |
|---|---|
| `processed/cohort.csv` | 920 standardized parents; identity, split, label, scaffold, patent family, publication bound and source evidence |
| `processed/exclusions.csv` | Excluded source records/parents; counts are not independent experimental samples |
| `processed/nearest_neighbors.csv` | Validation/test nearest-training comparisons using the prespecified chiral fingerprint |
| `predictions/final_predictions.csv` | 1,722 rows: seven methods × 246 test parents; scores are not calibrated hit probabilities |
| `predictions/bootstrap_draws.json` | Original 2,000 paired scaffold resamples for all methods |
| `analysis/rows.csv` | 1,722 model–parent error rows with source bins, similarity and descriptors |

The [public scoring manifest](../manifests/public_scoring_replay.json) binds replay inputs and expected outputs. Source HTML/images, checkpoints and caches are not included. Their historical hashes remain in manifests; a hash is provenance, not a claim that its file is distributed.
