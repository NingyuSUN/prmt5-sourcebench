"""Fixed-budget PRMT5 CPU baselines: train only, select/score validation only."""
from __future__ import annotations

import argparse
import csv
import copy
import json
import os
import platform
import tempfile
import time
import traceback
import warnings
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import joblib
import numpy as np
import rdkit
import sklearn
from sklearn.base import clone
from sklearn.ensemble import RandomForestClassifier
from sklearn.exceptions import ConvergenceWarning
from threadpoolctl import threadpool_limits

from benchmark_common import (COHORT_SHA256, PROTOCOL_SHA256, ConstantPrevalence,
    TieAwareKNN, choose_threshold, features, load_frozen_development,
    logistic_estimator, ranking_metrics, require_preflight,
    select_validation_candidate, sha256, threshold_metrics)

REPEAT_ATOL = 1e-12


def candidate_specs():
    specs = [{'candidate_id':'constant_prevalence', 'family':'constant_prevalence', 'chirality':False, 'params':{}},
             {'candidate_id':'five_neighbor', 'family':'five_neighbor', 'chirality':False, 'params':{'k':5}}]
    for chiral in (False, True):
        family = 'fingerprint_logistic_chiral' if chiral else 'fingerprint_logistic'
        for C in (.1,1,10):
            specs.append({'candidate_id':f'{family}_C{str(C).replace(".","p")}', 'family':family,
                          'chirality':chiral, 'params':{'C':C}})
    for leaf in (2,5):
        specs.append({'candidate_id':f'random_forest_leaf{leaf}', 'family':'random_forest',
                      'chirality':False, 'params':{'min_samples_leaf':leaf}})
    return specs


def make_model(spec):
    if spec not in candidate_specs():
        raise ValueError('candidate specification outside fixed protocol budget')
    family = spec['family']
    if family == 'constant_prevalence':
        return ConstantPrevalence()
    if family == 'five_neighbor':
        return TieAwareKNN(k=5)
    if family.startswith('fingerprint_logistic'):
        return logistic_estimator(spec['params']['C'])
    return RandomForestClassifier(n_estimators=500, min_samples_leaf=spec['params']['min_samples_leaf'],
        max_features='sqrt', class_weight='balanced', random_state=17, n_jobs=1)


def write_json(path, value):
    """Atomic JSON replacement; disallow nonfinite numeric artefacts."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary = tempfile.mkstemp(dir=path.parent, prefix='.'+path.name)
    try:
        with os.fdopen(descriptor, 'w') as handle:
            json.dump(value, handle, indent=2, allow_nan=False)
            handle.write('\n')
        os.replace(temporary, path)
    except BaseException:
        if os.path.exists(temporary):
            os.unlink(temporary)
        raise


def fit_twice(model, train_x, train_y, validation_x):
    """Repeat the same CPU fit, not another tuning trial; retain first fitted model."""
    with threadpool_limits(limits=1), warnings.catch_warnings():
        warnings.simplefilter('error', ConvergenceWarning)
        first = clone(model).fit(train_x, train_y)
        first_score = first.predict_proba(validation_x)[:,1]
        second = clone(model).fit(train_x, train_y)
        second_score = second.predict_proba(validation_x)[:,1]
    if not np.isfinite(first_score).all() or not np.isfinite(second_score).all():
        raise ValueError('nonfinite validation predictions')
    difference = float(np.max(np.abs(first_score-second_score)))
    if not np.allclose(first_score, second_score, rtol=0, atol=REPEAT_ATOL):
        raise ValueError(f'CPU repetition failed: max difference {difference}')
    return first, first_score, difference


def fit_candidates(development, output_dir):
    if set(development) != {'train','validation'}:
        raise ValueError('only train and validation groups may enter fitting')
    for split, rows in development.items():
        if not rows or any(row['split'] != split for row in rows):
            raise ValueError('split rows do not match their development role')
    train, validation = development['train'], development['validation']
    train_y = np.array([row['binary_y'] for row in train], dtype=np.int8)
    validation_y = np.array([row['binary_y'] for row in validation], dtype=np.int8)
    validation_ids = [row['compound_id'] for row in validation]
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    matrices = {chirality: (
        features([row['canonical_smiles'] for row in train], chirality=chirality),
        features([row['canonical_smiles'] for row in validation], chirality=chirality)) for chirality in (False,True)}
    candidate_records, predictions = [], {}
    for spec in candidate_specs():
        started = time.perf_counter()
        train_x, validation_x = matrices[spec['chirality']]
        model, score, difference = fit_twice(make_model(spec), train_x, train_y, validation_x)
        identifier = spec['candidate_id']
        joblib.dump(model, output_dir/(identifier+'.joblib'))
        prediction_path = output_dir/(identifier+'_validation.csv')
        with prediction_path.open('w', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=['compound_id','split','binary_y','score'])
            writer.writeheader()
            writer.writerows({'compound_id':row['compound_id'], 'split':'validation',
                              'binary_y':int(label), 'score':format(float(value),'.17g')}
                             for row,label,value in zip(validation,validation_y,score))
        predictions[identifier] = score
        candidate_records.append({**spec, 'fitted_split':'train', 'fitted_n':len(train),
            'selected_on_split':'validation', 'scored_splits':['validation'],
            'validation':ranking_metrics(validation_y, score, validation_ids),
            'runtime_seconds':time.perf_counter()-started,
            'cpu_repeated_fit_max_abs_difference':difference,
            'cpu_repeated_fit_atol':REPEAT_ATOL, 'cpu_repeated_fit_rtol':0,
            'model_file':identifier+'.joblib', 'validation_prediction_file':prediction_path.name,
            'model_parameters':model.get_params(deep=False)})
    selected = {}
    for family in dict.fromkeys(record['family'] for record in candidate_records):
        group = [record for record in candidate_records if record['family'] == family]
        if len(group) == 1:
            chosen = group[0]
        else:
            index = select_validation_candidate(validation_y,
                [predictions[record['candidate_id']] for record in group],
                [record['params'] for record in group],
                family='forest' if family == 'random_forest' else 'logistic')
            chosen = group[index]
        identifier = chosen['candidate_id']
        threshold = choose_threshold(validation_y, predictions[identifier], split='validation')
        selected[family] = {'candidate_id':identifier, 'params':chosen['params'],
            'chirality':chosen['chirality'], 'model_file':chosen['model_file'],
            'validation_prediction_file':chosen['validation_prediction_file'],
            'validation_average_precision':chosen['validation']['average_precision'],
            'threshold':threshold,
            'validation_threshold_metrics':threshold_metrics(validation_y, predictions[identifier], threshold)}
    primary = {key:value for key,value in selected.items() if key != 'fingerprint_logistic_chiral'}
    strongest = sorted(primary, key=lambda key: (-primary[key]['validation_average_precision'],key))[0]
    return {'status':'VALIDATION_MODELS_READY', 'candidate_count':len(candidate_records),
        'candidate_fit_repetitions':2, 'train_n':len(train), 'validation_n':len(validation),
        'training_label_counts':dict(Counter(str(int(label)) for label in train_y)),
        'training_class_weights':{str(label):len(train_y)/(2*int(np.sum(train_y==label))) for label in (0,1)},
        'fingerprint_policy':{'radius':2,'bits':2048,'primary_chirality':False,'sensitivity_chirality':True},
        'logistic_penalty':'L2 (sklearn1.9 l1_ratio=0)', 'scaler_policy':'fingerprint baselines unscaled; encoder-head scaler is train-only in benchmark_common',
        'scored_splits':['validation'], 'test_scored':False, 'candidates':candidate_records, 'selected':selected,
        'strongest_primary_classical_baseline':strongest,
        'strongest_classical_tie_rule':'highest validation AP; exact ties use lexicographic family name; chiral sensitivity excluded',
        'selection_rules':{'logistic':'highest validation AP; smaller C on exact ties',
                           'forest':'highest validation AP; larger leaf on exact ties',
                           'threshold':'maximum validation balanced accuracy; higher threshold on exact ties'}}


def public_report(report):
    """Publish portable metadata while retaining the original sealed run receipt."""
    published = copy.deepcopy(report)
    published.get('runtime', {}).pop('executable', None)
    return published


def run(root, preflight=None):
    root = Path(root).resolve()
    preflight = Path(preflight or root/'reports/models/encoder_preflight.json').resolve()
    # Gate before feature creation or estimator fitting. Missing or stale PASS fails closed.
    preflight_bindings = require_preflight(preflight)
    development, bindings = load_frozen_development(root)
    bindings.update(preflight_bindings)
    destination = root/'data/models/baselines'
    if destination.exists() or (root/'reports/models/validation_baselines.json').exists():
        raise FileExistsError(f'refusing to overwrite an existing baseline run: {destination}')
    destination.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(tempfile.mkdtemp(prefix='.baselines_pending_', dir=destination.parent))
    preflight_snapshot = json.loads(preflight.read_text()) if preflight.exists() else None
    if preflight_snapshot is not None:
        write_json(staging/'input_preflight_at_start.json', preflight_snapshot)
    started = datetime.now(timezone.utc).isoformat()
    write_json(staging/'RUNNING.json', {'status':'RUNNING','started_utc':started,'bindings':bindings})
    try:
        report = fit_candidates(development, staging)
        _, current = load_frozen_development(root)
        current.update(require_preflight(preflight))
        if {key:value for key,value in current.items() if key != 'encoder_preflight_sha256'} != {key:value for key,value in bindings.items() if key != 'encoder_preflight_sha256'}:
            raise ValueError('bound inputs changed during fitting')
        report.update({'protocol_sha256':PROTOCOL_SHA256,'cohort_sha256':COHORT_SHA256,
            'bindings':bindings, 'completion_input_bindings':current,
            'input_preflight_scope':'input gate PASS; encoder runtime status does not gate CPU models',
            'input_preflight_overall_status_at_start':preflight_snapshot.get('status') if preflight_snapshot else None,
            'started_utc':started,'completed_utc':datetime.now(timezone.utc).isoformat(),
            'runtime':{'python':platform.python_version(),'executable':os.sys.executable,
                       'numpy':np.__version__,'scikit_learn':sklearn.__version__,'rdkit':rdkit.__version__,
                       'joblib':joblib.__version__,'threadpool_limit':1,'random_forest_n_jobs':1},
            'code_sha256':{name:sha256(Path(__file__).parent/name) for name in ['benchmark_common.py','train_baselines.py','test_benchmark.py']}})
        write_json(staging/'validation_baselines.json', report)
        files = {str(path.relative_to(staging)):sha256(path) for path in sorted(staging.rglob('*')) if path.is_file()}
        seal = {'status':'SEALED_VALIDATION_ONLY','protocol_sha256':PROTOCOL_SHA256,'cohort_sha256':COHORT_SHA256,
                'bindings':bindings,'test_scored':False,'scored_splits':['validation'], 'files':files,
                'completed_utc':report['completed_utc']}
        write_json(staging/'SEALED.json', seal)
        os.replace(staging, destination)
        published_report = public_report(report)
        published_report['sealed_artifacts'] = 'data/models/baselines/SEALED.json'
        published_report['sealed_artifacts_sha256'] = sha256(destination/'SEALED.json')
        write_json(root/'reports/models/validation_baselines.json', published_report)
        return published_report
    except BaseException as error:
        failure_dir = staging if staging.exists() else destination
        if (failure_dir/'SEALED.json').exists():
            os.replace(failure_dir/'SEALED.json', failure_dir/'SEAL_INVALIDATED.json')
        write_json(failure_dir/'FAILED.json', {'status':'FAILED','bindings':bindings,
            'error_type':type(error).__name__, 'error':str(error), 'traceback':traceback.format_exc(),
            'failed_utc':datetime.now(timezone.utc).isoformat()})
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument('--preflight', type=Path)
    args = parser.parse_args()
    result = run(args.root, args.preflight)
    print(json.dumps({'status':result['status'],'test_scored':result['test_scored'],
                      'selected':result['selected']}, indent=2, allow_nan=False))


if __name__ == '__main__':
    main()
