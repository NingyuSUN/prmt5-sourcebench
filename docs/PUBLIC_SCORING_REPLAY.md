# Reproduce saved-prediction metrics

From a clean checkout, install the CPU baseline dependencies and run:

```sh
python -m pip install -r requirements/baseline.txt
python scripts/reproduce_metrics.py --output-dir outputs/replay
python -m unittest discover -s tests -p test_public_replay.py -v
```

Choose a new output directory for every replay. Existing output directories are
rejected. The commands require no model downloads, model weights, GPU, or network
access after dependency installation. RDKit is required by the original shared
module; the scoring calculations use NumPy and scikit-learn.

The entrypoint checks the byte hashes of the released cohort, saved predictions,
original final selection lock, six original frozen scripts, expected final
metrics, and expected bootstrap draws. It uses the original ranking,
classification, prediction-coverage and paired scaffold-bootstrap functions.
All seven methods are evaluated on the same 246 held-out compounds (69 positive,
177 negative), with the original validation-selected thresholds, 89 observed
scaffolds, 2,000 bootstrap draws, seed 20261016, and the original random-forest
comparator. Thresholds, model choices and test labels are not changed.

Every recomputed metric and bootstrap value is compared with the frozen result.
Only floating-point differences within relative and absolute tolerances of
1e-12 are allowed. The output contains `final_metrics.json`,
`bootstrap_draws.json`, and `verification.json`; the last file states whether
values also matched exactly. On the recorded baseline versions, the first two
files reproduce the original bytes.

This is a saved-prediction scoring replay. It does not verify model fitting,
regenerate predictions, reconstruct source tables, check the private training
receipt chain, or establish prospective biochemical or cellular validity. The
original protocol hash in the immutable selection lock identifies the experiment
freeze; any editorial publication of that protocol has a separate recorded
binding. Uncertainty is conditional on the single observed test patent family
and does not estimate variation across future sources.
