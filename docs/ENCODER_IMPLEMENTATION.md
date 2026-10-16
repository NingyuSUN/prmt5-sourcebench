# Encoder implementation and artifact contract

This implements protocol PRMT5-SSB-1, without changing its cohort or selection budget. Unit tests use a tiny synthetic encoder; passing them does not establish PRMT5 performance or pretrained-runtime readiness.

## Gates and fitting scope

`scripts/train_encoder.py` requires explicit `--root`, `--model-dir`, and `--device` (`cpu` or `cuda[:N]`). It uses `benchmark_common.load_frozen_development`, which returns 604 training and 70 validation parents and withholds final-test labels. The pinned input/runtime preflight and the reviewed local safetensors loader must pass before any classifier fitting or optimizer step. Model downloading and remote-code review belong to `prepare_encoder.py`; training uses its hash-verified `local_files_only` loader.

Only canonical nonisomeric SMILES are encoded. A hash-verified cache may contain vectors for all source structures; the fitting code selects development IDs and never reads final-test labels. Token truncation and unknown tokens fail the affected arm. Pooling is the reviewed encoder's `pooler_output`, an attention-mask mean including BOS/EOS.

## Frozen classifier and fair adaptation initialization

Three training-only StandardScaler/logistic pipelines use C = {0.1, 1, 10}. Validation AP selects C, with smaller C winning exact ties. The chosen scaler mean/scale are copied into fixed PyTorch buffers. Its coefficient and intercept initialize the same linear head for every adaptation run. Scaling/head arithmetic uses float64 to reproduce sklearn's frozen classifier; encoder weights and activations retain their original dtype. The scaler is never refit during adaptation.

Arithmetic equivalence on cached vectors is checked. Independently, every freshly loaded trial encoder must reproduce the chosen frozen classifier's 70 validation probabilities within atol/rtol 1e-5 before its optimizer exists. This checks fresh encoder buffers as well as the head transformation. Mismatch stops the arm; it does not silently initialize a different model.

## Partial adaptation and selection

Only `encoder.layer.10.*`, `encoder.layer.11.*`, and the linear head may learn. Actual names and parameter counts are saved. Shared/tied parameters reaching earlier blocks are rejected. The encoder stays in `eval()` throughout optimization, which disables dropout and keeps deterministic random features fixed while retaining autograd for the final two blocks. Earlier parameter hashes are compared before and after each trial.

Weighted binary cross-entropy uses weights N/(2Nclass) from training labels only. Effective batch size is 32; smaller micro batches accumulate losses divided by the current effective batch size. The final smaller batch keeps all training rows. No mixed precision, clipping, augmentation or architecture search is added.

AdamW uses encoder LR 1e-5 or 3e-5, head LR 1e-4 and weight decay 0.01. Each LR runs seeds {17, 29, 43}, at most 20 epochs with patience 4. Strict AP improvements select a copied checkpoint; tied AP retains the earlier epoch. LR selection averages the three seeds' selected-epoch validation AP and ties prefer lower LR. No best seed is selected. The primary adapted score averages the three selected-LR seed probabilities; its operating threshold uses validation labels only.

## Saved state and scoring interface

Default output is ignored `data/models/encoder/`. Existing nonempty runs are preserved; use a new output directory for another run.

- `frozen_head.npz`: training scaler and selected logistic coefficient/intercept; `frozen_head.json`: all C candidates, validation predictions/metrics and selected threshold.
- `lr_<rate>/seed_<seed>/best.pt` and `last.pt`: partial encoder parameter deltas, all encoder buffers, fixed scaler/head state, matching optimizer state, epoch, metrics, Python/NumPy/PyTorch RNG states and frozen source/model/code bindings. The original pinned base encoder is required; these are not standalone full-model weights.
- Each trial also saves its complete loss/AP trajectory, selected epoch, trainable-name/count evidence, validation probabilities and artifact hashes.
- `encoder_development.json`: all six trials, the mean-seed LR decision, ensemble validation metrics, sealed threshold and artifact checksums. It contains no final-test scores.

For later scoring, load the reviewed base encoder, create `ScaledLinearHead` from `frozen_head.npz`, configure the final two blocks to recover the expected parameter-name contract, load each selected seed's `best.pt` with `torch.load(...,weights_only=True)`, and call `restore_training_weights`. `predict_probabilities` evaluates structures; average the three seed probabilities and apply the already sealed threshold. Do not fit or select anything on final-test labels.

Source hashes are checked before/after development fitting; metadata JSON rejects NaN/Infinity. Partial files/checkpoints are retained on failure, and completion is reported only after all six trials and selection succeed. The report and checkpoint context record actual Python, NumPy, sklearn, RDKit, tokenizers, PyTorch and Transformers versions. RDKit/tokenizers must match the input gate; PyTorch/Transformers must match the runtime preflight. The runtime environment and model acquisition receipts remain separate evidence.

Example, from the canonical project root in the reviewed environment:

```bash
python scripts/train_encoder.py --root . --model-dir data/model/molformer-361063d0ad524ef77cf39b08469f6be770dc550f --device cuda --micro-batch-size 8
python -m unittest discover -s tests -p test_encoder_training.py -v
```

Model quality is reported in the completed results; archived inference uses the sealed final scoring stage and source-aware interpretation. This code does not turn retrospective biochemical labels into cellular, experimental or clinical validation.
