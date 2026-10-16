# Reproducing PRMT5 SourceBench

The public release supports a complete **saved-prediction metric replay** on CPU. It includes the original cohort, all seven methods' final scores, frozen validation thresholds, bootstrap draws, and the scoring implementation. Model training and raw-source reconstruction are separate operations with additional inputs.

## Run the public replay

Use Python 3.12 and a fresh virtual environment:

```bash
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements/baseline.txt
python scripts/reproduce_metrics.py --root . --output-dir outputs/replay
```

The command verifies the byte hashes of eleven frozen input/code artifacts, checks exactly seven methods × 246 test parents, and recomputes AP, ROC-AUC, tie-aware top-25 hits, validation-threshold confusion counts and 2,000 paired scaffold resamples. Every recomputed value is compared with the published artifact; float tolerance is 1e-12. Existing output directories are refused. No fitting, model selection, download or GPU operation occurs.

Outputs are `final_metrics.json`, `bootstrap_draws.json` and `verification.json`. The first two match the original bytes under the recorded baseline versions. [Detailed replay contract](PUBLIC_SCORING_REPLAY.md); [input manifest](../manifests/public_scoring_replay.json).

Run the CPU contract suites:

```bash
python -m unittest discover -s scripts -p 'test_*.py' -v
python -m unittest discover -s tests -p test_biological_analysis.py -v
python -m unittest discover -s tests -p test_public_replay.py -v
```

The separate synthetic encoder-training suite requires PyTorch; it uses a tiny encoder and does not download or train MoLFormer:

```bash
python -m pip install torch==2.6.0 --index-url https://download.pytorch.org/whl/cpu
python -m unittest discover -s tests -p test_encoder_training.py -v
```

Engineering tests do not establish biochemical validity or predictive usefulness. All test performance remains the frozen negative result.

## Artifact availability

| Layer | Included here | Additional requirements |
|---|---|---|
| Saved-prediction scoring | Complete cohort, seven prediction vectors, original thresholds, metric code, expected metrics/draws | CPU dependencies only |
| Biological interpretation | Similarity table, exclusions, per-compound errors, descriptor summaries and code | Existing frozen metrics; see the data/model cards for scientific limits |
| Source reconstruction | Acquisition manifests, original-table curation code, reviewed counts and hashes | Original ChEMBL snapshots, patent pages/images and reviewed manual transcription; these are not distributed |
| Model training/inference | Original fitting, prediction and lock code; all validation trial summaries | Pinned encoder assets, preprocessing artifacts, original protocol-bound experiment workspace and the separate training runtime |

The original 226-member source archive was verified to reconstruct the frozen cohort byte-for-byte. Eleven of fifteen output comparisons matched bytes; four auxiliary tables/receipts matched only under the documented ordering rules in [source replay evidence](../reports/curation/source_snapshot_replay.json). The source archive itself is not distributed. Full original source download-to-training reproduction is **not** provided by the public quick start.

All included cohort, predictions, bootstrap draws, original metrics, model-development records and six hash-bound scoring/training scripts preserve their original bytes. The public scientific protocol is an editorial copy with administrative statements removed; scientific definitions, data splits, model budgets and evaluation rules are unchanged. Both protocol hashes are in [protocol publication metadata](../manifests/protocol_publication.json). Historical locks retain the original protocol hash; they have not been rewritten to manufacture a new experiment identity. Accordingly, the archival `evaluate_locked.py` entry also requires the original protocol bytes and complete model/receipt chain. Use `reproduce_metrics.py` for this public release.

## Recorded training environments

| Operation | Recorded versions |
|---|---|
| Data / CPU models / scoring | Python 3.12.14; NumPy 2.5.4; SciPy 1.18.1; scikit-learn 1.9.1; RDKit 2026.09.1; joblib 1.6.0 |
| Encoder training / inference | Python 3.10.16; PyTorch 2.6.0+cu124; Transformers 5.12.1; tokenizers 0.22.2; NumPy 2.2.6; SciPy 1.15.3; scikit-learn 1.7.2; RDKit 2026.09.1 |

[`baseline.txt`](../requirements/baseline.txt) and [`encoder.txt`](../requirements/encoder.txt) describe separate recorded direct dependencies, not a complete platform/wheel lock. Do not deserialize baseline estimators in the encoder sklearn runtime. GPU work requires an actually passing CUDA tensor probe; a visible GPU name alone is insufficient. [`visualization.txt`](../requirements/visualization.txt) supplies the optional figure-rendering dependency.

The unchanged original training commands, applicable only in a complete protocol-bound reconstruction workspace, are:

```bash
python scripts/prepare_encoder.py --project-root . --model-dir data/model/molformer --stage assets --with-weights
python scripts/prepare_encoder.py --project-root . --model-dir data/model/molformer --stage inputs
python scripts/prepare_encoder.py --project-root . --model-dir data/model/molformer --stage runtime
python scripts/train_baselines.py --root . --preflight reports/models/encoder_preflight.json
python scripts/train_encoder.py --root . --model-dir data/model/molformer --device cuda --micro-batch-size 8
```

Use the appropriate interpreter for each command. These are an implementation reference, not a promise that a public checkout contains every prerequisite. Existing training outputs are protected; new experiments belong in a separate workspace. Final-test labels must not drive further model selection. The encoder loader verifies the pinned revision, source files, weights, tokenizer and all twelve pretrained random-feature buffers. [Implementation contract](ENCODER_IMPLEMENTATION.md).

## Source acquisition and limits

Original data manifests record URLs, request parameters and SHA256 values. ChEMBL's live API reflects its current release, so fetching today is not an exact replay of the recorded ChEMBL 37 snapshot. [Official release/API documentation](https://chembl.gitbook.io/chembl-interface-documentation/frequently-asked-questions/general-questions).

The archival `fetch_chembl.py`, `fetch_documents.py` and `fetch_compound_records.py` wrappers require an external client exposing `execute(request)`; that client is not included. There is no integrated downloader for all included patent pages or an automatic substitute for the reviewed 75-cell WO2019173804A1 transcription. Original-source reconstruction therefore requires restored snapshots or a new explicitly documented curation exercise. Missing inputs must not be replaced by silently changing expected hashes.

The pinned chemistry source is [DrugPerturb chemistry.py](https://github.com/NingyuSUN/drugperturb/blob/cfa256c201ef45ad54efdc78eb2cef529f211d51/src/drugperturb/chemistry.py), SHA256 `0a059a1ec5528ff140ad5634a277f6fb96022acbfec4ab9aa8f3bb450698627b`; retain its upstream MIT notice. Public metric replay does not require this module.

No complete fresh-environment live acquisition and GPU retraining run is claimed. Saved-prediction replay cannot verify how weights were learned, and it does not independently replicate an assay. One-family uncertainty, experimental-condition differences, source-date bounds and unknown pretraining membership remain scientific limitations regardless of engineering reproducibility.
