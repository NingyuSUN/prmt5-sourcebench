#!/usr/bin/env python3
"""Fit development-only frozen heads and last-two-block MoLFormer adaptation.

Real fitting requires the pinned cohort/protocol and a completed input/runtime
preflight. Unit tests inject a tiny encoder and never download model weights.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import random
import tempfile

import numpy as np
import torch
from torch import nn
from torch.nn import functional as F

C_GRID = (.1, 1., 10.)
LEARNING_RATES = (1e-5, 3e-5)
SEEDS = (17, 29, 43)
MAX_EPOCHS = 20
PATIENCE = 4
EFFECTIVE_BATCH = 32
HEAD_LR = 1e-4
WEIGHT_DECAY = .01
MAX_TOKENS = 202


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def snapshot_hashes(paths):
    return {name: sha256(path) for name, path in sorted(paths.items())}


def verify_snapshot(paths, expected):
    observed = snapshot_hashes(paths)
    require(observed == expected, 'Frozen source/model/code bytes changed during fitting')
    return observed


def labels_array(labels):
    labels = np.asarray(labels)
    require(labels.ndim == 1 and len(labels) > 0, 'Expected nonempty binary labels')
    require(set(np.unique(labels)) <= {0, 1}, 'Nonbinary labels')
    return labels.astype(np.int64)


def balanced_class_weights(labels):
    labels = labels_array(labels)
    counts = np.bincount(labels, minlength=2)
    require(np.all(counts > 0), 'Training needs both classes')
    return len(labels) / (2. * counts)


class ScaledLinearHead(nn.Module):
    """The selected logistic head, with a fixed training-only scaler.

    Float64 scaling/head arithmetic preserves the fitted sklearn reference.
    Encoder activations and parameters retain their original dtype.
    """
    def __init__(self, mean, scale, coefficient, intercept, device='cpu'):
        super().__init__()
        mean, scale = np.asarray(mean), np.asarray(scale)
        require(mean.ndim == 1 and mean.shape == scale.shape, 'Wrong scaler shape')
        require(np.isfinite(mean).all() and np.isfinite(scale).all() and (scale > 0).all(),
                'Invalid fixed scaler')
        self.register_buffer('mean', torch.as_tensor(mean.copy(), dtype=torch.float64))
        self.register_buffer('scale', torch.as_tensor(scale.copy(), dtype=torch.float64))
        self.linear = nn.Linear(len(mean), 1, dtype=torch.float64)
        with torch.no_grad():
            self.linear.weight.copy_(torch.as_tensor(np.asarray(coefficient).reshape(1, -1), dtype=torch.float64))
            self.linear.bias.copy_(torch.as_tensor(np.asarray(intercept).reshape(1), dtype=torch.float64))
        require(all(torch.isfinite(p).all() for p in self.parameters()), 'Invalid initial classifier')
        self.to(device)

    @classmethod
    def from_pipeline(cls, pipeline, device='cpu'):
        scaler, classifier = pipeline.named_steps['scaler'], pipeline.named_steps['classifier']
        require(np.array_equal(classifier.classes_, [0, 1]), 'Classifier positive-class convention changed')
        return cls(scaler.mean_, scaler.scale_, classifier.coef_, classifier.intercept_, device)

    def forward(self, features):
        require(features.ndim == 2 and features.shape[1] == len(self.mean), 'Wrong encoder feature shape')
        scaled = (features.to(dtype=self.mean.dtype) - self.mean) / self.scale
        return self.linear(scaled).squeeze(-1)


def select_frozen_head(train_features, train_labels, validation_features, validation_labels,
                       validation_ids, fit_fn, metric_fn):
    train_features = np.asarray(train_features, dtype=np.float64)
    validation_features = np.asarray(validation_features, dtype=np.float64)
    train_labels, validation_labels = labels_array(train_labels), labels_array(validation_labels)
    require(train_features.ndim == 2 and validation_features.ndim == 2 and
            train_features.shape[1] == validation_features.shape[1], 'Wrong feature matrices')
    require(len(train_features) == len(train_labels) and len(validation_features) == len(validation_labels),
            'Feature/label count mismatch')
    require(np.isfinite(train_features).all() and np.isfinite(validation_features).all(), 'Nonfinite features')
    candidates = []
    pipelines = {}
    for C in C_GRID:
        pipeline = fit_fn(train_features, train_labels, C=C, split='train')
        scaler = pipeline.named_steps['scaler']
        require(np.all(np.asarray(scaler.n_samples_seen_) == len(train_labels)), 'Scaler fitted outside training')
        require(np.allclose(scaler.mean_, train_features.mean(0), rtol=1e-12, atol=1e-12),
                'Scaler mean differs from training features')
        scores = np.asarray(pipeline.predict_proba(validation_features)[:, 1], dtype=float)
        ap = float(metric_fn(validation_labels, scores, validation_ids, fraction=.1)['average_precision'])
        require(math.isfinite(ap), 'Nonfinite frozen validation AP')
        pipelines[C] = pipeline
        candidates.append({'C': C, 'validation_average_precision': ap,
                           'validation_probabilities': scores.tolist()})
    best = max(candidates, key=lambda row: (row['validation_average_precision'], -row['C']))
    return {'pipeline': pipelines[best['C']], 'selected_C': best['C'], 'candidates': candidates}


def configure_last_two_blocks(model, head, expected_blocks=None):
    require(hasattr(model, 'encoder') and hasattr(model.encoder, 'layer'),
            'Reviewed encoder.layer structure missing')
    require(getattr(model.config, 'deterministic_eval', False), 'Deterministic encoder evaluation is required')
    blocks = model.encoder.layer
    require(len(blocks) >= 2, 'Encoder has fewer than two blocks')
    if expected_blocks is not None:
        require(len(blocks) == expected_blocks, 'Unexpected encoder block count')
    block_ids = [{id(parameter) for parameter in block.parameters()} for block in blocks]
    selected_ids = block_ids[-2] | block_ids[-1]
    require(bool(block_ids[-2]) and bool(block_ids[-1]), 'Empty selected transformer block')
    earlier_ids = set().union(*block_ids[:-2]) if len(blocks) > 2 else set()
    require(not (selected_ids & earlier_ids) and not (block_ids[-2] & block_ids[-1]),
            'Selected block has shared parameters with another block')
    for parameter in model.parameters():
        parameter.requires_grad_(id(parameter) in selected_ids)
    for parameter in head.parameters():
        parameter.requires_grad_(True)
    prefixes = tuple(f'encoder.layer.{i}.' for i in (len(blocks)-2, len(blocks)-1))
    names = [name for name, parameter in model.named_parameters() if parameter.requires_grad]
    require(names and all(name.startswith(prefixes) for name in names), 'Unexpected trainable encoder parameter')
    require({id(p) for p in model.parameters() if p.requires_grad} == selected_ids,
            'Selected parameter registration differs')
    model.eval()
    return {'block_path': 'encoder.layer', 'encoder_blocks': len(blocks),
            'trainable_block_indices': [len(blocks)-2, len(blocks)-1],
            'encoder_parameter_names': names,
            'head_parameter_names': [name for name, p in head.named_parameters() if p.requires_grad],
            'encoder_parameter_count': sum(p.numel() for p in model.parameters() if p.requires_grad),
            'head_parameter_count': sum(p.numel() for p in head.parameters() if p.requires_grad),
            'frozen_encoder_parameter_count': sum(p.numel() for p in model.parameters() if not p.requires_grad),
            'encoder_mode': 'eval_with_gradients', 'fixed_scaler_buffers': ['mean', 'scale'],
            'head_dtype': str(head.linear.weight.dtype)}


def make_optimizer(model, head, learning_rate):
    encoder_parameters = [p for p in model.parameters() if p.requires_grad]
    head_parameters = [p for p in head.parameters() if p.requires_grad]
    require(bool(encoder_parameters) and bool(head_parameters), 'Missing trainable parameter group')
    return torch.optim.AdamW([{'params': encoder_parameters, 'lr': learning_rate},
                              {'params': head_parameters, 'lr': HEAD_LR}], weight_decay=WEIGHT_DECAY)


def token_batch(tokenizer, smiles, device):
    values = tokenizer(list(smiles), padding=True, truncation=False, return_tensors='pt')
    require('input_ids' in values and 'attention_mask' in values, 'Tokenizer masks/IDs missing')
    require(values['input_ids'].shape[1] <= MAX_TOKENS, 'Overlength molecule; truncation is forbidden')
    require(not bool((values['input_ids'] == tokenizer.unk_token_id).any()), 'Unknown tokenizer input')
    return {key: value.to(device) for key, value in values.items()}


def pooled_features(model, inputs):
    pooled = model(**inputs).pooler_output
    require(pooled.ndim == 2 and torch.isfinite(pooled).all(), 'Invalid pooled encoder output')
    return pooled


def train_epoch(model, head, tokenizer, smiles, labels, optimizer, class_weights, device,
                micro_batch_size=8):
    labels = labels_array(labels)
    require(len(smiles) == len(labels), 'Training count mismatch')
    require(1 <= micro_batch_size <= EFFECTIVE_BATCH, 'Invalid micro batch')
    model.eval()  # Disables dropout/random-feature train behavior without disabling autograd.
    head.train()
    order = torch.randperm(len(labels)).tolist()
    total_loss = 0.
    for start in range(0, len(order), EFFECTIVE_BATCH):
        effective = order[start:start+EFFECTIVE_BATCH]
        optimizer.zero_grad(set_to_none=True)
        for offset in range(0, len(effective), micro_batch_size):
            indices = effective[offset:offset+micro_batch_size]
            inputs = token_batch(tokenizer, [smiles[i] for i in indices], device)
            logits = head(pooled_features(model, inputs))
            target = torch.as_tensor(labels[indices], dtype=logits.dtype, device=device)
            weights = torch.as_tensor(class_weights, dtype=logits.dtype, device=device)[target.long()]
            loss = F.binary_cross_entropy_with_logits(logits, target, weight=weights, reduction='sum')
            require(bool(torch.isfinite(loss)), 'Nonfinite training loss')
            (loss / len(effective)).backward()
            total_loss += float(loss.detach().cpu())
        gradients = [p.grad for group in optimizer.param_groups for p in group['params'] if p.grad is not None]
        require(gradients and bool(torch.stack([torch.isfinite(g).all() for g in gradients]).all()),
                'Nonfinite or missing training gradients')
        optimizer.step()
    return total_loss / len(labels)


def predict_probabilities(model, head, tokenizer, smiles, device, batch_size=32):
    model.eval()
    head.eval()
    output = []
    with torch.no_grad():
        for start in range(0, len(smiles), batch_size):
            logits = head(pooled_features(model, token_batch(tokenizer, smiles[start:start+batch_size], device)))
            require(bool(torch.isfinite(logits).all()), 'Nonfinite classifier logits')
            output.append(torch.sigmoid(logits).cpu().numpy())
    probabilities = np.concatenate(output)
    require(np.isfinite(probabilities).all(), 'Nonfinite predicted probabilities')
    return probabilities


def verify_initialization(model, head, tokenizer, smiles, reference_probabilities, device):
    require(reference_probabilities is not None, 'Frozen initialization reference is required')
    reference = np.asarray(reference_probabilities, dtype=float)
    observed = predict_probabilities(model, head, tokenizer, smiles, device)
    require(reference.shape == observed.shape and np.isfinite(reference).all(),
            'Invalid frozen initialization reference')
    require(np.allclose(observed, reference, atol=1e-5, rtol=1e-5),
            'Fresh encoder initialization differs from the selected frozen classifier')
    return float(np.max(np.abs(observed-reference)))


def fixed_feature_map_digest(model):
    digest = hashlib.sha256()
    for name, buffer in model.named_buffers():
        if name.endswith('.feature_map.weight'):
            digest.update(name.encode())
            digest.update(buffer.detach().cpu().contiguous().view(torch.uint8).numpy().tobytes())
    return digest.hexdigest()


def checkpoint_primitives(value):
    """Rebuild metadata with builtin types; never pickle version/scalar classes."""
    if isinstance(value, np.generic):
        return checkpoint_primitives(value.item())
    if isinstance(value, np.ndarray):
        return checkpoint_primitives(value.tolist())
    if value is None:
        return None
    if isinstance(value, str):
        return str(value)  # TorchVersion is a str subclass with an unsafe pickle global.
    if isinstance(value, bool):
        return bool(value)
    if isinstance(value, int):
        return int(value)
    if isinstance(value, float):
        require(math.isfinite(value), 'Nonfinite checkpoint metadata')
        return float(value)
    if isinstance(value, dict):
        require(all(isinstance(key, str) for key in value), 'Checkpoint metadata keys must be strings')
        return {str(key): checkpoint_primitives(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [checkpoint_primitives(item) for item in value]
    raise TypeError('Unsupported checkpoint metadata type: ' + type(value).__name__)


def cpu_copy(value):
    if torch.is_tensor(value):
        return value.detach().cpu().clone()
    if isinstance(value, dict):
        return {key: cpu_copy(item) for key, item in value.items()}
    if isinstance(value, list):
        return [cpu_copy(item) for item in value]
    if isinstance(value, tuple):
        return tuple(cpu_copy(item) for item in value)
    return value


def capture_rng(device):
    numpy_state = np.random.get_state()
    return {'python': random.getstate(),
            'numpy': {'algorithm': str(numpy_state[0]), 'keys': [int(value) for value in numpy_state[1]],
                      'position': int(numpy_state[2]), 'has_gaussian': int(numpy_state[3]),
                      'cached_gaussian': float(numpy_state[4])},
            'torch_cpu': torch.get_rng_state().clone(),
            'torch_cuda': [state.cpu().clone() for state in torch.cuda.get_rng_state_all()]
                          if torch.device(device).type == 'cuda' else []}


def capture_training_state(model, head, optimizer, epoch, metrics, device):
    return {'format': 'partial_encoder_delta_v1', 'epoch': int(epoch), 'metrics': checkpoint_primitives(metrics),
            'encoder_state': {name: cpu_copy(p) for name, p in model.named_parameters() if p.requires_grad},
            'encoder_buffers': {name: cpu_copy(buffer) for name, buffer in model.named_buffers()},
            'head_state': cpu_copy(head.state_dict()), 'optimizer_state': cpu_copy(optimizer.state_dict()),
            'rng_state': capture_rng(device)}


def restore_training_weights(model, head, checkpoint):
    parameters = {name: p for name, p in model.named_parameters() if p.requires_grad}
    require(set(parameters) == set(checkpoint['encoder_state']), 'Checkpoint trainable encoder names differ')
    buffers = dict(model.named_buffers())
    require(set(buffers) == set(checkpoint['encoder_buffers']), 'Checkpoint encoder buffer names differ')
    with torch.no_grad():
        for name, parameter in parameters.items():
            parameter.copy_(checkpoint['encoder_state'][name].to(parameter.device))
        for name, buffer in buffers.items():
            buffer.copy_(checkpoint['encoder_buffers'][name].to(buffer.device))
    head.load_state_dict(checkpoint['head_state'], strict=True)
    model.eval()


def frozen_parameter_digest(model):
    digest = hashlib.sha256()
    for name, parameter in model.named_parameters():
        if not parameter.requires_grad:
            digest.update(name.encode())
            digest.update(str((tuple(parameter.shape), parameter.dtype)).encode())
            digest.update(parameter.detach().cpu().contiguous().view(torch.uint8).numpy().tobytes())
    return digest.hexdigest()


def seed_everything(seed, device):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.device(device).type == 'cuda':
        os.environ.setdefault('CUBLAS_WORKSPACE_CONFIG', ':4096:8')
        torch.cuda.manual_seed_all(seed)
        torch.backends.cuda.matmul.allow_tf32 = False
        torch.backends.cudnn.allow_tf32 = False
    torch.use_deterministic_algorithms(True)


def atomic_save(path, writer):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode='wb', dir=path.parent, prefix='.'+path.name+'.', delete=False) as handle:
            temporary = Path(handle.name)
            writer(handle)
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def save_json(path, value):
    atomic_save(path, lambda handle: handle.write((json.dumps(value, indent=2, sort_keys=True, allow_nan=False)+'\n').encode()))


def run_trial(model, tokenizer, head, train_smiles, train_labels, validation_smiles, validation_labels,
              validation_ids, learning_rate, seed, device, output_dir, metric_fn,
              max_epochs=MAX_EPOCHS, patience=PATIENCE, micro_batch_size=8, expected_blocks=None,
              checkpoint_context=None, initial_validation_probabilities=None):
    require(1 <= max_epochs <= MAX_EPOCHS and patience >= 1, 'Trial epoch/patience limit invalid')
    safe_context = checkpoint_primitives(checkpoint_context or {})
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    seed_everything(seed, device)
    trainable = configure_last_two_blocks(model, head, expected_blocks)
    initial_difference = verify_initialization(model, head, tokenizer, validation_smiles,
                                               initial_validation_probabilities, device)
    before_frozen = frozen_parameter_digest(model)
    before_feature_maps = fixed_feature_map_digest(model)
    optimizer = make_optimizer(model, head, learning_rate)
    weights = balanced_class_weights(train_labels)
    best_ap, best_state, best_scores, waiting = -math.inf, None, None, 0
    trace = []
    for epoch in range(1, max_epochs+1):
        loss = train_epoch(model, head, tokenizer, train_smiles, train_labels, optimizer, weights,
                           device, micro_batch_size)
        scores = predict_probabilities(model, head, tokenizer, validation_smiles, device)
        ap = float(metric_fn(validation_labels, scores, validation_ids, fraction=.1)['average_precision'])
        require(math.isfinite(ap), 'Nonfinite adapted validation AP')
        improved = ap > best_ap  # Tied epochs retain the earlier checkpoint.
        metrics = {'average_precision': ap, 'weighted_training_loss': float(loss)}
        state = capture_training_state(model, head, optimizer, epoch, metrics, device)
        state['context'] = dict(safe_context)
        state['context'].update({'initial_validation_max_abs_probability_difference': initial_difference, 'seed': seed, 'encoder_learning_rate': learning_rate,
                                 'head_learning_rate': HEAD_LR, 'weight_decay': WEIGHT_DECAY,
                                 'trainable': trainable})
        state['context'] = checkpoint_primitives(state['context'])
        atomic_save(output_dir/'last.pt', lambda handle: torch.save(state, handle))
        if improved:
            best_ap, best_state, best_scores, waiting = ap, state, scores.copy(), 0
            atomic_save(output_dir/'best.pt', lambda handle: torch.save(best_state, handle))
        else:
            waiting += 1
        trace.append({'epoch': epoch, 'weighted_training_loss': float(loss),
                      'validation_average_precision': ap, 'is_new_best': improved,
                      'consecutive_nonimproving_epochs': waiting})
        save_json(output_dir/'trajectory.json', trace)
        print(json.dumps({'stage':'partial_fine_tuning', 'learning_rate':learning_rate,
                          'seed':seed, **trace[-1]}, allow_nan=False), flush=True)
        if waiting >= patience:
            break
    require(best_state is not None, 'No selected epoch')
    restore_training_weights(model, head, best_state)
    require(frozen_parameter_digest(model) == before_frozen, 'Frozen earlier parameters changed')
    require(fixed_feature_map_digest(model) == before_feature_maps, 'Deterministic random feature-map buffers changed')
    result = {'learning_rate': learning_rate, 'seed': seed, 'best_epoch': best_state['epoch'],
              'best_validation_average_precision': best_ap, 'training_trace': trace,
              'initial_validation_max_abs_probability_difference': initial_difference,
              'initialization_atol':1e-5, 'initialization_rtol':1e-5,
              'fixed_feature_map_sha256_before_and_after':before_feature_maps,
              'validation_probabilities': best_scores.tolist(), 'trainable': trainable,
              'frozen_parameters_sha256_before_and_after': before_frozen,
              'balanced_training_class_weights': weights.tolist(),
              'artifacts': {'best.pt':sha256(output_dir/'best.pt'),
                            'last.pt':sha256(output_dir/'last.pt'),
                            'trajectory.json':sha256(output_dir/'trajectory.json')}}
    save_json(output_dir/'trial.json', result)
    return result


def select_learning_rate(trials):
    require(len(trials) == len(LEARNING_RATES)*len(SEEDS), 'Missing learning rates or seeds')
    means = {}
    for learning_rate in LEARNING_RATES:
        group = [row for row in trials if row['learning_rate'] == learning_rate]
        require(sorted(row['seed'] for row in group) == sorted(SEEDS), 'Missing or duplicate seeds')
        scores = [float(row['best_validation_average_precision']) for row in group]
        require(np.isfinite(scores).all(), 'Nonfinite seed validation AP')
        means[learning_rate] = float(np.mean(scores))
    selected = max(LEARNING_RATES, key=lambda lr: (means[lr], -lr))
    return {'selected_learning_rate': selected,
            'mean_validation_average_precision_by_learning_rate': {str(lr):means[lr] for lr in LEARNING_RATES},
            'selection_rule': 'Mean best-epoch validation AP across all three seeds; ties lower LR; no best-seed selection',
            'selected_seeds': list(SEEDS)}


def nonisomeric_smiles(rows):
    from rdkit import Chem
    result = []
    for row in rows:
        molecule = Chem.MolFromSmiles(row['canonical_smiles'])
        require(molecule is not None, 'Invalid frozen development parent')
        result.append(Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=False))
    return result


def development_features(root, preflight, model, tokenizer, rows, smiles, device):
    """Read only unlabelled cached vectors and select development IDs, or encode development inputs."""
    runtime = preflight['runtime_gate']
    if runtime.get('feature_cache'):
        path = root / runtime['feature_cache']
        require(sha256(path) == runtime['features_sha256'], 'Frozen feature-cache hash mismatch')
        with np.load(path, allow_pickle=False) as cached:
            ids = cached['compound_ids'].astype(str).tolist()
            require(len(ids) == len(set(ids)), 'Duplicate cached feature IDs')
            lookup = {compound_id:i for i,compound_id in enumerate(ids)}
            indices = [lookup[row['compound_id']] for row in rows]
            features = np.asarray(cached['features'][indices], dtype=np.float64)
            require(cached['splits'][indices].astype(str).tolist() == [row['split'] for row in rows],
                    'Cached development split mismatch')
    else:
        from prepare_encoder import extract_features
        features = extract_features(model, tokenizer, smiles, device, 32).astype(np.float64)
    require(features.shape == (len(rows),768) and np.isfinite(features).all(), 'Invalid development features')
    return features


def save_predictions(path, ids, scores):
    require(len(ids) == len(scores), 'Prediction ID/score count mismatch')
    def write(handle):
        text = 'compound_id,probability\n' + ''.join(f'{compound_id},{float(score):.17g}\n'
                                                    for compound_id,score in zip(ids,scores))
        handle.write(text.encode())
    atomic_save(path, write)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--model-dir', type=Path, required=True)
    parser.add_argument('--device', required=True, help='Explicit cpu, cuda, or cuda:N; no automatic fallback')
    parser.add_argument('--micro-batch-size', type=int, default=8, choices=(1,2,4,8,16,32))
    parser.add_argument('--output-dir', type=Path)
    args = parser.parse_args()
    root, model_dir = args.root.resolve(), args.model_dir.resolve()
    device = torch.device(args.device)
    require(device.type in ('cpu','cuda'), 'Only the explicit reviewed CPU/CUDA execution paths are supported')
    if device.type == 'cuda':
        require(torch.cuda.is_available(), 'Requested CUDA is unavailable; no CPU fallback')
    from benchmark_common import load_frozen_development, require_preflight, fit_scaled_logistic
    from benchmark_common import ranking_metrics, choose_threshold
    from prepare_encoder import load_encoder_bundle, REVISION, EXPECTED_FILES
    # These gates run before fitting any classifier or taking any optimizer step.
    development, bindings = load_frozen_development(root)
    preflight_path = root/'reports/models/encoder_preflight.json'
    bindings.update(require_preflight(preflight_path))
    preflight = json.loads(preflight_path.read_text())
    require(preflight.get('status') == 'PASS' and
            preflight.get('runtime_gate',{}).get('status') == 'PASS', 'Encoder runtime gate has not passed')
    require(preflight['runtime_gate'].get('repeat_pass') and
            preflight['runtime_gate'].get('different_batch_pass'), 'Encoder feature determinism has not passed')
    seed_everything(SEEDS[0], device)
    model, tokenizer, loader = load_encoder_bundle(model_dir, str(device))
    require(loader['revision'] == REVISION and loader['deterministic_eval'], 'Unreviewed encoder contract')
    require(preflight['model_asset_manifest_sha256'] == loader['asset_manifest_sha256'] and
            preflight['weights_sha256'] == loader['weights_sha256'], 'Preflight/model bytes differ')
    runtime_loader = preflight['runtime_gate']['loader']
    require(all(loader[key] == runtime_loader[key] for key in ('torch','transformers')),
            'Encoder runtime versions differ from the reviewed preflight')
    source_paths = {name:root/name for name in
        ('docs/PROTOCOL.md','data/processed/cohort.csv','reports/curation/cohort_gate.json',
         'reports/models/encoder_preflight.json','scripts/train_encoder.py',
         'scripts/benchmark_common.py','scripts/prepare_encoder.py')}
    source_paths.update({'model_assets/'+name:model_dir/name for name in
                         list(EXPECTED_FILES)+['asset_manifest.json','model.safetensors']})
    if preflight['runtime_gate'].get('feature_cache'):
        name = preflight['runtime_gate']['feature_cache']
        source_paths[name] = root/name
    source_hashes_before = snapshot_hashes(source_paths)
    for source_name, binding_name in (
            ('docs/PROTOCOL.md','protocol_sha256'),
            ('data/processed/cohort.csv','cohort_sha256'),
            ('reports/curation/cohort_gate.json','cohort_gate_sha256'),
            ('reports/models/encoder_preflight.json','encoder_preflight_sha256')):
        require(source_hashes_before[source_name] == bindings[binding_name],
                'Frozen source changed between loading and snapshot: ' + source_name)
    require(source_hashes_before['model_assets/asset_manifest.json'] == loader['asset_manifest_sha256'] and
            source_hashes_before['model_assets/model.safetensors'] == loader['weights_sha256'],
            'Reviewed model assets changed between loading and snapshot')
    import platform
    import sklearn
    import tokenizers
    from rdkit import rdBase
    versions = {'python':platform.python_version(), 'numpy':np.__version__,
                'sklearn':sklearn.__version__, 'rdkit':rdBase.rdkitVersion,
                'tokenizers':tokenizers.__version__, 'torch':str(torch.__version__),
                'transformers':loader['transformers']}
    input_gate = preflight['input_gate']
    require(versions['rdkit'] == input_gate['rdkit_version'] and
            versions['tokenizers'] == input_gate['tokenizers_version'],
            'Representation runtime differs from the verified input gate')
    rows = development['train'] + development['validation']
    smiles = nonisomeric_smiles(rows)
    features = development_features(root, preflight, model, tokenizer, rows, smiles, str(device))
    n_train = len(development['train'])
    train_labels = labels_array([row['binary_y'] for row in development['train']])
    validation_labels = labels_array([row['binary_y'] for row in development['validation']])
    validation_ids = [row['compound_id'] for row in development['validation']]
    selected = select_frozen_head(features[:n_train], train_labels, features[n_train:], validation_labels,
                                  validation_ids, fit_scaled_logistic, ranking_metrics)
    pipeline = selected['pipeline']
    reference_head = ScaledLinearHead.from_pipeline(pipeline, str(device))
    with torch.inference_mode():
        reference_logits = reference_head(torch.as_tensor(features, device=device)).cpu().numpy()
    require(np.allclose(reference_logits, pipeline.decision_function(features), rtol=1e-5, atol=1e-5),
            'Adaptation initialization does not reproduce the selected frozen classifier')
    output = args.output_dir or root/'data/models/encoder'
    require(not output.exists() or not any(output.iterdir()), 'Output directory already contains a run; preserve it and use a new directory')
    output.mkdir(parents=True, exist_ok=True)
    frozen_score = pipeline.predict_proba(features[n_train:])[:,1]
    frozen_summary = {'selected_C': selected['selected_C'], 'candidates': selected['candidates'],
                      'runtime_versions':versions, 'bindings':bindings,
                      'validation_metrics':ranking_metrics(validation_labels, frozen_score, validation_ids),
                      'selected_threshold':choose_threshold(validation_labels, frozen_score, split='validation'),
                      'scaler_fit_split':'train', 'feature_dtype':'float64',
                      'initial_logit_equivalence_atol':1e-5, 'initial_logit_equivalence_rtol':1e-5}
    atomic_save(output/'frozen_head.npz', lambda handle: np.savez_compressed(handle,
        mean=pipeline.named_steps['scaler'].mean_, scale=pipeline.named_steps['scaler'].scale_,
        coefficient=pipeline.named_steps['classifier'].coef_, intercept=pipeline.named_steps['classifier'].intercept_))
    save_json(output/'frozen_head.json', frozen_summary)
    save_predictions(output/'frozen_validation_predictions.csv', validation_ids, frozen_score)
    trials = []
    context = {'bindings':bindings, 'model_revision':REVISION, 'model_receipt':loader,
               'runtime_versions':versions,
               'source_hashes_before':source_hashes_before,
               'code_sha256':sha256(Path(__file__)), 'fixed_scaler_fit_split':'train',
               'frozen_head_C':selected['selected_C']}
    del reference_head, model
    for learning_rate in LEARNING_RATES:
        for seed in SEEDS:
            verify_snapshot(source_paths, source_hashes_before)
            seed_everything(seed, device)
            model, tokenizer, current_loader = load_encoder_bundle(model_dir, str(device))
            require(current_loader['asset_manifest_sha256'] == loader['asset_manifest_sha256'], 'Model assets changed between trials')
            head = ScaledLinearHead.from_pipeline(pipeline, str(device))
            trial_dir = output/f'lr_{learning_rate:g}'/f'seed_{seed}'
            trial = run_trial(model, tokenizer, head, smiles[:n_train], train_labels, smiles[n_train:],
                validation_labels, validation_ids, learning_rate, seed, str(device), trial_dir, ranking_metrics,
                micro_batch_size=args.micro_batch_size, expected_blocks=12, checkpoint_context=context,
                initial_validation_probabilities=frozen_score)
            trial['artifact_directory'] = str(trial_dir.relative_to(output))
            save_predictions(trial_dir/'validation_predictions.csv', validation_ids, trial['validation_probabilities'])
            trials.append(trial)
            del head, model
            if device.type == 'cuda':
                torch.cuda.empty_cache()
    source_hashes_after = verify_snapshot(source_paths, source_hashes_before)
    selection = select_learning_rate(trials)
    chosen = [trial for trial in trials if trial['learning_rate'] == selection['selected_learning_rate']]
    ensemble = np.mean([trial['validation_probabilities'] for trial in chosen], axis=0)
    threshold = choose_threshold(validation_labels, ensemble, split='validation')
    save_predictions(output/'adapted_validation_predictions.csv', validation_ids, ensemble)
    report = {'status':'DEVELOPMENT_COMPLETE_NO_FINAL_TEST_SCORES', 'bindings':bindings,
              'protocol':'PRMT5-SSB-1', 'device':str(device), 'model_receipt':loader,
              'runtime_versions':versions,
              'source_hashes_before':source_hashes_before, 'source_hashes_after':source_hashes_after,
              'frozen_head':frozen_summary, 'partial_fine_tuning':selection, 'trials':trials,
              'adapted_validation_metrics':ranking_metrics(validation_labels, ensemble, validation_ids),
              'selected_threshold':threshold, 'primary_prediction_rule':'Mean probability across seeds17,29,43 at selected LR',
              'training_parents':len(train_labels), 'validation_parents':len(validation_labels),
              'fixed_budget':{'max_epochs':MAX_EPOCHS,'patience':PATIENCE,'effective_batch':EFFECTIVE_BATCH,
                              'micro_batch':args.micro_batch_size,'encoder_learning_rates':list(LEARNING_RATES),
                              'head_learning_rate':HEAD_LR,'weight_decay':WEIGHT_DECAY,'seeds':list(SEEDS)},
              'code_sha256':sha256(Path(__file__)),
              'artifacts':{str(path.relative_to(output)):sha256(path) for path in sorted(output.rglob('*')) if path.is_file()}}
    save_json(output/'encoder_development.json', report)
    print(json.dumps({'status':report['status'], 'selected_C':selected['selected_C'],
                      'selected_learning_rate':selection['selected_learning_rate'],
                      'selected_seeds':selection['selected_seeds']}, indent=2, allow_nan=False), flush=True)


if __name__ == '__main__':
    main()
