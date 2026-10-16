#!/usr/bin/env python3
"""Generate locked test scores without fitting, label access, or evaluation metrics.

Baselines must run in their original sklearn environment. Encoder inference uses
reviewed local weights and the fixed training scaler, without sklearn unpickling.
"""
from __future__ import annotations
import argparse
from collections import Counter
import csv
import hashlib
import io
import json
import os
from pathlib import Path
import shutil
import tempfile

import numpy as np

PROTOCOL_SHA256 = '283d38f275ef4c91d17dc87f1049213936543353d873198ba219b9a0b5890d1f'
COHORT_SHA256 = '17772cc5f360464abb2d5b913f1055ca60921d71315dfc03aa93baa1d8d89c96'
BASELINE_NAMES = {'constant_prevalence', 'five_neighbor', 'fingerprint_logistic',
                  'fingerprint_logistic_chiral', 'random_forest'}
ENCODER_NAMES = {'frozen_encoder', 'adapted_encoder'}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def inside(root, relative):
    relative = Path(relative)
    require(not relative.is_absolute(), 'Artifact path must be project-relative')
    path = (Path(root) / relative).resolve()
    require(path.is_relative_to(Path(root).resolve()), 'Artifact path escapes project')
    return path


def test_inputs(rows, expected_counts=None):
    """Project explicit structure fields only; never inspect target columns."""
    ids, counts, result = set(), Counter(), []
    for row in rows:
        split, identifier = row['split'], row['compound_id']
        require(split in ('train', 'validation', 'test'), 'Unknown cohort split')
        require(identifier and identifier not in ids, 'Missing or duplicate compound ID')
        ids.add(identifier)
        counts[split] += 1
        if split == 'test':
            require(bool(row['canonical_smiles']), 'Missing test structure')
            result.append({'compound_id': identifier, 'split': 'test',
                           'canonical_smiles': row['canonical_smiles']})
    if expected_counts is not None:
        require(dict(counts) == expected_counts, 'Frozen split counts changed')
    require(bool(result), 'No test inputs')
    return result


def load_test_inputs(root):
    root = Path(root)
    require(sha256(root/'docs/PROTOCOL.md') == PROTOCOL_SHA256, 'Frozen protocol hash mismatch')
    cohort = (root/'data/processed/cohort.csv').read_bytes()
    require(hashlib.sha256(cohort).hexdigest() == COHORT_SHA256, 'Frozen cohort hash mismatch')
    # Whole-file hashing binds the curated freeze. CSV values for the target
    # column are never indexed, parsed, copied into model inputs, or returned.
    rows = csv.DictReader(io.StringIO(cohort.decode('utf-8'), newline=''))
    return test_inputs(rows, {'train':604, 'validation':70, 'test':246})


def checked_scores(scores, count):
    scores = np.asarray(scores, dtype=float)
    require(scores.shape == (count,), 'Wrong score count or shape')
    require(np.isfinite(scores).all() and ((scores >= 0) & (scores <= 1)).all(),
            'Scores must be finite probabilities in [0,1]')
    return scores


def mean_seed_scores(seed_scores, count):
    require(set(seed_scores) == {17, 29, 43}, 'All three fixed seeds are required')
    return np.mean([checked_scores(seed_scores[seed], count) for seed in (17,29,43)], axis=0)


def prediction_rows(inputs, predictions, expected_names):
    require(set(predictions) == set(expected_names), 'Wrong selected model coverage')
    ids = [row['compound_id'] for row in inputs]
    require(len(ids) == len(set(ids)), 'Duplicate prediction IDs')
    output = []
    for name in sorted(predictions):
        scores = checked_scores(predictions[name], len(inputs))
        output.extend({'model':name, 'compound_id':identifier, 'split':'test',
                       'score':format(float(score), '.17g')} for identifier, score in zip(ids,scores))
    return output


def require_runtime(observed, expected, keys):
    for key in keys:
        require(key in observed and key in expected and str(observed[key]) == str(expected[key]),
                'Runtime version differs from development: '+key)


def verified_artifact(root, relative, artifact_hashes):
    require(relative in artifact_hashes, 'Artifact absent from locked development hashes')
    path = inside(root, relative)
    require(sha256(path) == artifact_hashes[relative], 'Locked development artifact changed: '+relative)
    return path


def baseline_predictions(root, inputs, models, artifact_hashes):
    import joblib
    import sklearn
    from rdkit import rdBase
    from threadpoolctl import threadpool_limits
    from benchmark_common import features
    require(set(models) == BASELINE_NAMES, 'Five prespecified baseline families required')
    report_path = 'data/models/baselines/validation_baselines.json'
    report = json.loads(verified_artifact(root, report_path, artifact_hashes).read_text())
    versions = {'numpy':np.__version__, 'scikit_learn':sklearn.__version__,
                'rdkit':rdBase.rdkitVersion, 'joblib':joblib.__version__}
    require_runtime(versions, report['runtime'], versions)
    matrices, predictions, consumed = {}, {}, {report_path:artifact_hashes[report_path]}
    with threadpool_limits(limits=1):
        for name, spec in sorted(models.items()):
            chirality = spec['chirality']
            require(isinstance(chirality, bool) and chirality == (name == 'fingerprint_logistic_chiral'),
                    'Fingerprint chirality differs from fixed model arm')
            relative = spec['model_relative_path']
            path = verified_artifact(root, relative, artifact_hashes)
            consumed[relative] = artifact_hashes[relative]
            # Loading is permitted only after same-runtime and sealed hash checks.
            model = joblib.load(path)
            require(np.array_equal(model.classes_, [0,1]), 'Positive class convention changed')
            if chirality not in matrices:
                matrices[chirality] = features([row['canonical_smiles'] for row in inputs], chirality=chirality)
            result = np.asarray(model.predict_proba(matrices[chirality]))
            require(result.shape == (len(inputs), 2), 'Unexpected baseline probability matrix')
            predictions[name] = checked_scores(result[:,1], len(inputs))
    return predictions, {'runtime':versions, 'consumed_model_artifacts':consumed}


def encoder_predictions(root, model_dir, inputs, models, artifact_hashes, device):
    import torch
    import tokenizers
    from rdkit import rdBase
    from prepare_encoder import load_encoder_bundle, REVISION
    from train_encoder import (ScaledLinearHead, configure_last_two_blocks, restore_training_weights,
        nonisomeric_smiles, predict_probabilities, fixed_feature_map_digest, frozen_parameter_digest, seed_everything)
    require(set(models) == ENCODER_NAMES, 'Both prespecified encoder arms required')
    require(device in ('cpu','cuda') or device.startswith('cuda:'), 'Unsupported explicit inference device')
    # Match fitting's deterministic numeric policy before any CUDA work.
    os.environ['CUBLAS_WORKSPACE_CONFIG'] = ':4096:8'
    seed_everything(17, device)
    if torch.device(device).type == 'cuda':
        require(torch.cuda.is_available(), 'Requested CUDA unavailable; no fallback')
    preflight_path = 'reports/models/encoder_preflight.json'
    preflight = json.loads(verified_artifact(root, preflight_path, artifact_hashes).read_text())
    require(preflight.get('status') == 'PASS' and preflight['runtime_gate']['status'] == 'PASS',
            'Completed encoder runtime preflight required')
    require(preflight['protocol_sha256'] == PROTOCOL_SHA256 and preflight['cohort_sha256'] == COHORT_SHA256,
            'Encoder preflight input binding mismatch')
    require(preflight['runtime_gate']['repeat_pass'] and preflight['runtime_gate']['different_batch_pass'],
            'Encoder deterministic evaluation preflight required')
    report_relative = 'data/models/encoder/encoder_development.json'
    report = json.loads(verified_artifact(root, report_relative, artifact_hashes).read_text())
    versions = {'numpy':np.__version__, 'torch':str(torch.__version__),
                'rdkit':rdBase.rdkitVersion, 'tokenizers':tokenizers.__version__}
    require_runtime(versions, report['runtime_versions'], versions)
    require(versions['rdkit'] == preflight['input_gate']['rdkit_version'] and
            versions['tokenizers'] == preflight['input_gate']['tokenizers_version'],
            'Representation runtime differs from input audit')
    model, tokenizer, loader = load_encoder_bundle(model_dir, device)
    require(loader['revision'] == REVISION and loader['deterministic_eval'], 'Reviewed deterministic encoder required')
    require(loader['asset_manifest_sha256'] == preflight['model_asset_manifest_sha256'] and
            loader['weights_sha256'] == preflight['weights_sha256'], 'Encoder weights differ from preflight')
    require_runtime(loader, preflight['runtime_gate']['loader'], ('torch','transformers'))
    require_runtime(loader, report['model_receipt'], ('torch','transformers'))
    consumed = {preflight_path:artifact_hashes[preflight_path], report_relative:artifact_hashes[report_relative]}
    for spec in models.values():
        require(spec['model_asset_sha'] == loader['asset_manifest_sha256'], 'Selected encoder asset hash mismatch')
    head_relative = models['frozen_encoder']['head_path']
    require(models['adapted_encoder']['head_path'] == head_relative, 'Adaptation must share fixed frozen head source')
    head_path = verified_artifact(root, head_relative, artifact_hashes)
    consumed[head_relative] = artifact_hashes[head_relative]
    with np.load(head_path, allow_pickle=False) as arrays:
        require(set(arrays.files) == {'mean','scale','coefficient','intercept'}, 'Unexpected frozen head fields')
        head_arrays = {key:arrays[key].copy() for key in arrays.files}
    def make_head():
        return ScaledLinearHead(**head_arrays, device=device)
    smiles = nonisomeric_smiles(inputs)
    frozen_head = make_head()
    predictions = {'frozen_encoder':checked_scores(
        predict_probabilities(model, frozen_head, tokenizer, smiles, device), len(inputs))}
    # Repeat a structure-only subset before restoring deltas; no labels/metrics.
    probe = predict_probabilities(model, frozen_head, tokenizer, smiles[:2], device, batch_size=1)
    require(np.allclose(probe, predictions['frozen_encoder'][:2], atol=1e-5, rtol=1e-5),
            'Initial encoder prediction is batch-dependent')
    spec = models['adapted_encoder']
    require(spec['seeds'] == [17,29,43], 'Selected adapted ensemble must retain all three seeds')
    require(spec['learning_rate'] in (1e-5,3e-5), 'Unprescribed learning rate')
    checkpoints = spec['checkpoint_paths']
    require(len(checkpoints) == 3 and len(set(checkpoints)) == 3, 'Three distinct selected checkpoints required')
    seed_predictions, details = {}, []
    for relative in checkpoints:
        path = verified_artifact(root, relative, artifact_hashes)
        consumed[relative] = artifact_hashes[relative]
        checkpoint = torch.load(path, map_location='cpu', weights_only=True)
        require(checkpoint['format'] == 'partial_encoder_delta_v1', 'Unexpected checkpoint format')
        context = checkpoint['context']
        seed = context['seed']
        require(seed in spec['seeds'] and seed not in seed_predictions, 'Unexpected or duplicate checkpoint seed')
        require(context['encoder_learning_rate'] == spec['learning_rate'], 'Checkpoint learning rate differs from selection')
        require(context['model_revision'] == REVISION and
                context['model_receipt']['asset_manifest_sha256'] == loader['asset_manifest_sha256'],
                'Checkpoint base model differs from reviewed weights')
        require(context['bindings']['protocol_sha256'] == PROTOCOL_SHA256 and
                context['bindings']['cohort_sha256'] == COHORT_SHA256, 'Checkpoint cohort/protocol binding mismatch')
        head = make_head()
        trainable = configure_last_two_blocks(model, head, expected_blocks=12)
        require(trainable == context['trainable'], 'Checkpoint trainable block contract differs')
        before_features, before_frozen = fixed_feature_map_digest(model), frozen_parameter_digest(model)
        restore_training_weights(model, head, checkpoint)
        require(fixed_feature_map_digest(model) == before_features, 'Checkpoint changes fixed random features')
        require(frozen_parameter_digest(model) == before_frozen, 'Checkpoint changes frozen earlier layers')
        require(torch.equal(head.mean, torch.as_tensor(head_arrays['mean'], device=device, dtype=torch.float64)) and
                torch.equal(head.scale, torch.as_tensor(head_arrays['scale'], device=device, dtype=torch.float64)),
                'Checkpoint changes fixed training scaler')
        seed_predictions[seed] = checked_scores(predict_probabilities(model, head, tokenizer, smiles, device), len(inputs))
        details.append({'seed':seed, 'epoch':checkpoint['epoch'], 'checkpoint':relative,
                        'sha256':consumed[relative], 'fixed_feature_map_sha256':before_features})
    predictions['adapted_encoder'] = mean_seed_scores(seed_predictions, len(inputs))
    return predictions, {'runtime':versions, 'model_loader':loader,
                         'consumed_model_artifacts':consumed, 'selected_checkpoints':details,
                         'initial_batch_difference':float(np.max(np.abs(probe-predictions['frozen_encoder'][:2]))),
                         'seed_probabilities':{str(seed):scores.tolist() for seed,scores in seed_predictions.items()}}


def publish_output(destination, rows, receipt):
    destination = Path(destination)
    require(not destination.exists(), 'Prediction output exists; refusing overwrite')
    destination.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix='.predictions_pending_', dir=destination.parent))
    try:
        path = staging/'predictions.csv'
        with path.open('w', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=['model','compound_id','split','score'])
            writer.writeheader()
            writer.writerows(rows)
        receipt = dict(receipt, predictions_sha256=sha256(path))
        (staging/'receipt.json').write_text(json.dumps(receipt, indent=2, sort_keys=True, allow_nan=False)+'\n')
        require(not destination.exists(), 'Prediction output appeared during computation')
        os.rename(staging, destination)
    except BaseException:
        shutil.rmtree(staging, ignore_errors=True)
        raise
    return receipt


def run(root, kind, model_dir=None, device=None):
    from selection_lock import load_locked_model_selection
    root = Path(root).resolve()
    require(kind in ('baselines','encoder'), 'Unknown prediction arm')
    destination = root/'data/predictions'/kind
    require(not destination.exists(), 'Prediction output exists; refusing overwrite')
    lock_relative = 'manifests/model_selection_lock.json'
    lock = load_locked_model_selection(root)
    lock_hash = sha256(root/lock_relative)
    inputs = load_test_inputs(root)
    names = BASELINE_NAMES if kind == 'baselines' else ENCODER_NAMES
    models = {name:lock['models'][name] for name in sorted(names)}
    source_hashes = {**lock['code_sha256'], **lock['development_artifacts'],
                     lock_relative:lock_hash, 'docs/PROTOCOL.md':PROTOCOL_SHA256,
                     'data/processed/cohort.csv':COHORT_SHA256}
    # Recheck before execution in case an artifact changed after lock loading.
    for relative, digest in source_hashes.items():
        require(sha256(inside(root, relative)) == digest, 'Locked input changed before prediction: '+relative)
    if kind == 'baselines':
        require(device is None and model_dir is None, 'Baseline arm uses its local fitted model paths')
        predictions, details = baseline_predictions(root, inputs, models, lock['development_artifacts'])
    else:
        require(model_dir is not None and device is not None, 'Encoder requires explicit --model-dir and --device')
        model_dir = Path(model_dir).resolve()
        asset_files = {str(path.relative_to(model_dir)):sha256(path)
                       for path in model_dir.iterdir() if path.is_file()}
        predictions, details = encoder_predictions(root, model_dir, inputs, models,
                                                   lock['development_artifacts'], device)
        require(asset_files == {str(path.relative_to(model_dir)):sha256(path)
                                for path in model_dir.iterdir() if path.is_file()},
                'Encoder assets changed during prediction')
        details['model_asset_files_sha256'] = asset_files
    rows = prediction_rows(inputs, predictions, names)
    for relative, digest in source_hashes.items():
        require(sha256(inside(root, relative)) == digest, 'Locked input changed during prediction: '+relative)
    # Re-derive all validation-only choices before publishing any output.
    require(load_locked_model_selection(root) == lock, 'Model selection changed during prediction')
    receipt = {'status':'LOCKED_TEST_PREDICTIONS_NO_METRICS', 'kind':kind,
               'protocol_sha256':PROTOCOL_SHA256, 'cohort_sha256':COHORT_SHA256,
               'model_selection_lock_sha256':lock_hash, 'code_sha256':lock['code_sha256'],
               'development_artifacts':lock['development_artifacts'],
               'models':models, 'test_input_count':len(inputs), 'prediction_row_count':len(rows),
               'compound_ids':[row['compound_id'] for row in inputs],
               'target_columns_accessed':False, 'fitting_performed':False, 'metrics_computed':False,
               **details}
    return publish_output(destination, rows, receipt)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--kind', choices=('baselines','encoder'), required=True)
    parser.add_argument('--model-dir', type=Path)
    parser.add_argument('--device', help='Encoder only: explicit cpu, cuda, or cuda:N; no fallback')
    args = parser.parse_args()
    receipt = run(args.root, args.kind, args.model_dir, args.device)
    print(json.dumps({key:receipt[key] for key in
        ('status','kind','test_input_count','prediction_row_count','predictions_sha256')}))


if __name__ == '__main__':
    main()
