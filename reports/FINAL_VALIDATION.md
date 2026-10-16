# Final validation record

Verified 2026-10-16T16:22:57+09:00 (Asia/Tokyo). The fixed experiment is complete. Engineering correctness and scientific interpretation are assessed separately.

## Engineering and reproduction

- **84 tests passed:** 61 data/CPU/locking/metric contracts plus 23 encoder-training and post-score-analysis tests. These use small synthetic examples; they do not establish drug activity.
- All ten CPU candidate fits were repeated with identical validation predictions; saved-model reload reproduced every validation CSV hash.
- Real encoder preflight covered all 920 inputs, pinned source/weights and fixed random-feature buffers. All six adaptation trials completed; frozen-parameter digests and initialization equivalence passed. Safe checkpoint serialization was verified through actual `weights_only=True` loading.
- Model choices and both inference arms were sealed before final scoring. All seven models cover the same 246 parents; no final-test tuning or threshold change occurred.
- The two saved-prediction scoring outputs replay exactly. The source snapshot rebuilds the exact frozen cohort: 11 of 15 comparison outputs are byte-identical; the other four have explicitly verified ordering-only differences. See [source replay](curation/source_snapshot_replay.json) and [scoring replay](results/scoring_replay.json).
- The final result figure was checked against the frozen metrics.

## Independent source and scientific checks

A separate read-only review linked the original SK A1 biochemical activity table to all 248 retained example links representing 246 test parents. Earlier parsed tables matched all 611 retained example links to the biochemical activity column; multiple aliases do not add independent compounds. Strong `+++` and `A` labels were correctly positive. The review kept cellular readouts separate.

Independent ID-based pairing reproduced all seven AP, AUC, top-25 and confusion-count results. An independently implemented tie-grouped AP algorithm reproduced all 2,000 paired scaffold resamples and interval endpoints within 1e-12. All three seed vectors reconstructed the sealed ensemble within 1e-15. No conclusion-changing label, pairing or scoring defect was identified.

The supported conclusion is absence of useful screening enrichment and no reliable improvement over the preselected random forest **in this held-out family**. It does not establish statistical equivalence, universal model failure, or a causal effect of SAM concentration or stereochemistry. Source, assay, chemistry and date remain confounded.

## Remaining limits

Fresh-environment installation and complete live acquisition are unverified. Every chemical figure and earliest worldwide structure–label disclosure were not independently reconstructed. Exact pretraining membership is unknown. No new wet assay, prospective screen, cellular efficacy, selectivity, ADMET or clinical validation was performed. Negative predictive results are preserved in [English](RESULTS.md) and [Chinese](RESULTS.zh-CN.md).

Machine-readable record: [validation_final.json](validation_final.json). Original weights and source snapshots are not distributed. Released cohort and predictions support the separate public metric replay described in the reproduction guide.
