"""Frozen PRMT5 benchmark contracts: development-only fitting and tie-aware metrics."""
from __future__ import annotations

import csv
import hashlib
import io
import json
import math
from collections import Counter
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from sklearn.base import BaseEstimator, ClassifierMixin
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

PROTOCOL_SHA256 = '283d38f275ef4c91d17dc87f1049213936543353d873198ba219b9a0b5890d1f'
ENCODER_REVISION = '361063d0ad524ef77cf39b08469f6be770dc550f'
COHORT_SHA256 = '17772cc5f360464abb2d5b913f1055ca60921d71315dfc03aa93baa1d8d89c96'
EXPECTED_COUNTS = {'train': 604, 'validation': 70, 'test': 246}
LOGISTIC_C = (0.1, 1.0, 10.0)
FOREST_LEAVES = (2, 5)


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _vectors(y, score):
    y = np.asarray(y)
    score = np.asarray(score, dtype=float)
    if y.ndim != 1 or score.ndim != 1 or len(y) == 0 or len(y) != len(score):
        raise ValueError('nonempty equal-length one-dimensional labels and scores required')
    if not np.isin(y, [0, 1]).all() or not np.isfinite(score).all():
        raise ValueError('binary labels and finite scores required')
    return y.astype(np.int8), score


def _matrix(x):
    x = np.asarray(x)
    if x.ndim != 2 or len(x) == 0 or x.shape[1] == 0 or not np.isfinite(x).all():
        raise ValueError('nonempty finite feature matrix required')
    return x


def _train_arrays(x, y, split='train'):
    if split != 'train':
        raise ValueError('fitting is restricted to train rows')
    x = _matrix(x)
    y, _ = _vectors(y, np.zeros(len(x)))
    if len(np.unique(y)) != 2:
        raise ValueError('both training classes required')
    return x, y


def features(smiles, chirality=False):
    """Morgan radius2,2048bits; achiral policy explicitly removes stereochemistry."""
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048, includeChirality=chirality)
    fingerprints = []
    for text in smiles:
        molecule = Chem.MolFromSmiles(text)
        if molecule is None:
            raise ValueError(f'invalid structure: {text}')
        if not chirality:
            molecule = Chem.MolFromSmiles(Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=False))
        fingerprints.append(generator.GetFingerprintAsNumPy(molecule))
    if not fingerprints:
        raise ValueError('no structures supplied')
    return np.asarray(fingerprints, dtype=np.uint8)


def ranking_metrics(y, score, ids, fraction=.1):
    """AP/AUC and expected/min/max hits under uniform ordering of the boundary tie."""
    y, score = _vectors(y, score)
    ids = [str(identifier) for identifier in ids]
    if len(ids) != len(y) or len(set(ids)) != len(ids):
        raise ValueError('one unique compound ID per row required')
    if not np.isfinite(fraction) or not 0 < fraction <= 1:
        raise ValueError('screening fraction must be in (0,1]')
    n = len(y)
    positive = int(y.sum())
    ap = float(average_precision_score(y, score)) if positive else 0.0
    auc = float(roc_auc_score(y, score)) if 0 < positive < n else None
    k = min(n, max(1, math.ceil(fraction*n)))
    boundary = float(np.sort(score)[-k])
    above = score > boundary
    tied = score == boundary
    above_n = int(above.sum())
    above_hits = int(y[above].sum())
    tie_n = int(tied.sum())
    tie_hits = int(y[tied].sum())
    slots = k-above_n
    order = sorted(range(n), key=lambda i: (-score[i], ids[i]))[:k]
    return {
        'n': n, 'positive': positive, 'negative': n-positive,
        'prevalence': positive/n, 'average_precision': ap, 'roc_auc': auc,
        'single_class': positive in (0,n),
        'screening': {
            'fraction': float(fraction), 'k': k, 'boundary_score': boundary,
            'strictly_above_count': above_n, 'strictly_above_hits': above_hits,
            'boundary_tie_count': tie_n, 'boundary_tie_hits': tie_hits,
            'slots_in_boundary_tie': slots,
            'expected_hits': above_hits + slots*tie_hits/tie_n,
            'minimum_hits': above_hits + max(0, slots-(tie_n-tie_hits)),
            'maximum_hits': above_hits + min(slots, tie_hits),
            'random_selection_expected_hits': k*positive/n,
            'deterministic_ids': [ids[i] for i in order],
            'deterministic_hits': int(y[order].sum()),
            'arbitrary_order_within_boundary_tie': tie_n > slots,
            'ordering_rule': 'descending score, then lexicographic compound_id; ID order is arbitrary within ties',
            'expected_hits_interpretation': 'fractional expectation over boundary-tie order, not an observed integer count'
        }
    }


def threshold_metrics(y, score, threshold):
    y, score = _vectors(y, score)
    if not np.isfinite(threshold):
        raise ValueError('threshold must be finite')
    predicted = score >= threshold
    tp = int(np.sum(predicted & (y == 1)))
    fp = int(np.sum(predicted & (y == 0)))
    fn = int(np.sum(~predicted & (y == 1)))
    tn = int(np.sum(~predicted & (y == 0)))
    sensitivity = tp/(tp+fn) if tp+fn else None
    specificity = tn/(tn+fp) if tn+fp else None
    return {'threshold': float(threshold), 'n': len(y), 'tp': tp, 'fn': fn, 'fp': fp, 'tn': tn,
            'sensitivity': sensitivity, 'specificity': specificity,
            'balanced_accuracy': (sensitivity+specificity)/2 if sensitivity is not None and specificity is not None else None}


def choose_threshold(y, score, split='validation'):
    if split != 'validation':
        raise ValueError('threshold selection requires validation split')
    y, score = _vectors(y, score)
    if len(np.unique(y)) != 2:
        raise ValueError('threshold selection requires both validation classes')
    thresholds = np.r_[np.unique(score), np.nextafter(score.max(), np.inf)]
    if not np.isfinite(thresholds).all():
        raise ValueError('score magnitude prevents finite all-negative threshold')
    # Compare the equivalent integer objective TP*N_negative + TN*N_positive.
    # This avoids floating-point drift falsely resolving a balanced-accuracy tie.
    positives, negatives = int(y.sum()), int((1-y).sum())
    best = None
    for threshold in thresholds:
        predicted = score >= threshold
        tp = int(np.sum(predicted & (y == 1)))
        tn = int(np.sum(~predicted & (y == 0)))
        candidate = (tp*negatives + tn*positives, float(threshold))
        if best is None or candidate > best:
            best = candidate
    return best[1]


def logistic_estimator(C):
    if C not in LOGISTIC_C:
        raise ValueError('C outside frozen budget')
    # sklearn1.9 deprecates explicit penalty='l2'; l1_ratio=0 is its L2 equivalent.
    return LogisticRegression(C=C, l1_ratio=0.0, class_weight='balanced', max_iter=5000,
                              solver='lbfgs', random_state=17)


def fit_scaled_logistic(x, y, C, split='train'):
    x, y = _train_arrays(x, y, split)
    return Pipeline([('scaler', StandardScaler()), ('classifier', logistic_estimator(C))]).fit(x, y)


def select_validation_candidate(y, predictions, configurations, family, split='validation'):
    if split != 'validation':
        raise ValueError('model selection requires validation split')
    if len(predictions) != len(configurations) or not configurations:
        raise ValueError('matching nonempty predictions/configurations required')
    ranking = []
    for index, (score, config) in enumerate(zip(predictions, configurations)):
        yy, ss = _vectors(y, score)
        if len(np.unique(yy)) != 2:
            raise ValueError('model selection requires both validation classes')
        if family == 'logistic' and set(config) == {'C'} and config['C'] in LOGISTIC_C:
            preference = -config['C']
        elif family == 'forest' and set(config) == {'min_samples_leaf'} and config['min_samples_leaf'] in FOREST_LEAVES:
            preference = config['min_samples_leaf']
        else:
            raise ValueError('configuration outside frozen selection budget')
        ranking.append((float(average_precision_score(yy, ss)), preference, -index))
    return -max(ranking)[2]


class ConstantPrevalence(ClassifierMixin, BaseEstimator):
    def fit(self, x, y):
        x, y = _train_arrays(x, y)
        self.prevalence_ = float(y.mean())
        self.n_features_in_ = x.shape[1]
        self.classes_ = np.array([0,1])
        return self

    def predict_proba(self, x):
        x = _matrix(x)
        if x.shape[1] != self.n_features_in_:
            raise ValueError('feature dimension mismatch')
        return np.tile([1-self.prevalence_, self.prevalence_], (len(x),1))


class TieAwareKNN(ClassifierMixin, BaseEstimator):
    def __init__(self, k=5):
        self.k = k

    def fit(self, x, y):
        x, y = _train_arrays(x, y)
        if not isinstance(self.k, int) or self.k < 1 or self.k > len(y):
            raise ValueError('neighbor count must be within training size')
        if not np.isin(x, [0,1]).all():
            raise ValueError('Tanimoto requires binary fingerprints')
        self.x_ = x.astype(np.float64)
        self.y_ = y.copy()
        self.n_features_in_ = x.shape[1]
        self.classes_ = np.array([0,1])
        return self

    def predict_proba(self, x):
        x = _matrix(x)
        if x.shape[1] != self.n_features_in_ or not np.isin(x, [0,1]).all():
            raise ValueError('invalid query fingerprint shape or values')
        x = x.astype(np.float64)
        intersection = x @ self.x_.T
        union = x.sum(axis=1)[:,None] + self.x_.sum(axis=1)[None,:] - intersection
        similarity = np.divide(intersection, union, out=np.ones_like(union), where=union != 0)
        probabilities = []
        for row in similarity:
            boundary = np.partition(row, -self.k)[-self.k]
            probabilities.append(float(self.y_[row >= boundary].mean()))
        probability = np.asarray(probabilities)
        return np.c_[1-probability, probability]


def development_rows(rows):
    """Only train/validation labels are read; test label strings are never parsed."""
    ids = set()
    output = {'train': [], 'validation': []}
    for row in rows:
        if row.get('split') not in EXPECTED_COUNTS:
            raise ValueError(f"unknown split: {row.get('split')}")
        identifier = row.get('compound_id')
        if not identifier or identifier in ids:
            raise ValueError('missing/duplicate compound_id')
        ids.add(identifier)
        if row['split'] == 'test':
            continue
        if str(row.get('binary_y')) not in ('0','1'):
            raise ValueError('invalid development label')
        clean = dict(row)
        clean['binary_y'] = int(row['binary_y'])
        output[row['split']].append(clean)
    return output


def load_frozen_development(root):
    root = Path(root)
    protocol = root/'docs/PROTOCOL.md'
    cohort = root/'data/processed/cohort.csv'
    gate_path = root/'reports/curation/cohort_gate.json'
    protocol_bytes = protocol.read_bytes()
    if hashlib.sha256(protocol_bytes).hexdigest() != PROTOCOL_SHA256:
        raise ValueError('frozen protocol hash mismatch')
    cohort_bytes = cohort.read_bytes()
    if hashlib.sha256(cohort_bytes).hexdigest() != COHORT_SHA256:
        raise ValueError('frozen cohort hash mismatch')
    gate_bytes = gate_path.read_bytes()
    gate = json.loads(gate_bytes)
    if gate.get('status') != 'PASS_FOR_NARROW_SOURCE_SHIFT_BINARY_BENCHMARK':
        raise ValueError('cohort gate is not PASS')
    if gate.get('outputs', {}).get('data/processed/cohort.csv') != COHORT_SHA256:
        raise ValueError('cohort gate/output hash mismatch')
    rows = list(csv.DictReader(io.StringIO(cohort_bytes.decode('utf-8'), newline='')))
    dev = development_rows(rows)
    counts = Counter(row['split'] for row in rows)
    if dict(counts) != EXPECTED_COUNTS:
        raise ValueError('frozen split count mismatch')
    # Hashing raw bytes binds the entire freeze, but no test labels leave this loader.
    bindings = {'protocol_sha256': PROTOCOL_SHA256, 'cohort_sha256': COHORT_SHA256,
                'cohort_gate_sha256': hashlib.sha256(gate_bytes).hexdigest()}
    return dev, bindings


def require_preflight(path):
    """CPU gate is input PASS plus pinned bindings/revision, independent of GPU runtime."""
    raw = Path(path).read_bytes()
    receipt = json.loads(raw)
    if receipt.get('protocol_sha256') != PROTOCOL_SHA256:
        raise ValueError('preflight protocol hash mismatch')
    if receipt.get('cohort_sha256') != COHORT_SHA256:
        raise ValueError('preflight cohort hash mismatch')
    input_gate = receipt.get('input_gate', {})
    if input_gate.get('status') != 'PASS':
        raise ValueError('encoder input preflight must explicitly report input PASS')
    revision = receipt.get('model_revision', receipt.get('revision'))
    if revision != ENCODER_REVISION or receipt.get('revision', revision) != ENCODER_REVISION:
        raise ValueError('preflight encoder revision mismatch')
    return {'encoder_preflight_sha256': hashlib.sha256(raw).hexdigest(),
            'encoder_input_gate_sha256': hashlib.sha256(json.dumps(input_gate, sort_keys=True).encode()).hexdigest()}
