#!/usr/bin/env python3
"""Score sealed final-test predictions; this entry never fits or selects a model."""
import argparse
import collections
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from sklearn.metrics import average_precision_score
from benchmark_common import ranking_metrics, threshold_metrics
from selection_lock import validate_selection_lock


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_csv(path):
    with Path(path).open(newline='') as handle:
        return list(csv.DictReader(handle))


def paired_cluster_bootstrap(y, predictions, groups, comparator, draws=2000, seed=20261016):
    """Resample observed scaffold groups with replacement, paired across methods.

    Resampling compounds independently would understate correlation among members
    of a scaffold. This still cannot measure variation across future patent sources.
    """
    y = np.asarray(y, dtype=int)
    groups = np.asarray(groups)
    scores = {name: np.asarray(value, dtype=float) for name, value in predictions.items()}
    require(len(y) == len(groups) and len(y) > 0, 'Invalid bootstrap rows')
    require(set(np.unique(y)) <= {0, 1}, 'Nonbinary labels')
    require(comparator in scores, 'Missing frozen comparator')
    require(all(v.shape == y.shape and np.isfinite(v).all() for v in scores.values()), 'Invalid scores')
    require(draws > 0, 'Positive number of draws required')
    unique = np.unique(groups)
    members = [np.flatnonzero(groups == group) for group in unique]
    rng = np.random.default_rng(seed)
    sampled_ap = {name: [] for name in scores}
    valid_draw_ids = []
    invalid = 0
    for draw in range(draws):
        sampled_groups = rng.integers(0, len(unique), size=len(unique))
        indices = np.concatenate([members[i] for i in sampled_groups])
        sampled_y = y[indices]
        if len(np.unique(sampled_y)) < 2:
            invalid += 1
            continue
        valid_draw_ids.append(draw)
        for name, value in scores.items():
            sampled_ap[name].append(float(average_precision_score(sampled_y, value[indices])))
    require(bool(valid_draw_ids), 'No valid two-class bootstrap resample')
    reference = np.array(sampled_ap[comparator])
    results = {}
    for name, values in sampled_ap.items():
        array = np.asarray(values)
        results[name] = {
            'average_precision_percentile_95': np.quantile(array, [.025, .975]).tolist(),
            'paired_ap_difference_vs_comparator_percentile_95': np.quantile(array-reference, [.025, .975]).tolist(),
            'point_ap_difference_vs_comparator': float(average_precision_score(y, scores[name])-average_precision_score(y, scores[comparator])),
        }
    return {
        'resampling_unit': 'Observed final-test Murcko scaffold, all members retained with multiplicity',
        'scope': 'Conditional uncertainty within one patent family; not uncertainty across future sources',
        'seed': seed, 'requested_draws': draws, 'valid_draws': len(valid_draw_ids),
        'invalid_single_class_draws': invalid, 'scaffold_count': len(unique),
        'comparator': comparator, 'models': results,
    }, {'draw_ids': valid_draw_ids, **sampled_ap}


def validate_prediction_rows(rows, test_rows, model_names):
    expected = {r['compound_id'] for r in test_rows}
    require(len(expected) == len(test_rows), 'Duplicate test IDs')
    result = collections.defaultdict(dict)
    for row in rows:
        name, cid = row['model'], row['compound_id']
        require(name in model_names, 'Unselected model in final predictions')
        require(row['split'] == 'test', 'Final scorer accepts test predictions only')
        require(cid in expected, 'Prediction not in frozen test cohort')
        require(cid not in result[name], 'Duplicate prediction')
        score = float(row['score'])
        require(math.isfinite(score) and 0 <= score <= 1, 'Invalid classifier score')
        result[name][cid] = score
    require(set(result) == set(model_names), 'Missing selected model')
    require(all(set(values) == expected for values in result.values()), 'Unequal prediction coverage')
    return {name: np.array([result[name][r['compound_id']] for r in test_rows]) for name in model_names}


def score_frozen(root, lock_path):
    root = Path(root)
    lock = json.loads(Path(lock_path).read_text())
    validate_selection_lock(root, lock)
    prediction_path = root/lock['predictions_path']
    all_rows = load_csv(root/'data/processed/cohort.csv')
    test_rows = [r for r in all_rows if r['split'] == 'test']
    require(len(test_rows) == 246 and len({r['source_family'] for r in test_rows}) == 1, 'Unexpected test source')
    require(sum(int(r['binary_y']) for r in test_rows) == 69, 'Unexpected test labels')
    names = list(lock['models'])
    require(len(names) == len(set(names)), 'Duplicated model list')
    scores = validate_prediction_rows(load_csv(prediction_path), test_rows, names)
    y = np.array([int(r['binary_y']) for r in test_rows])
    ids = [r['compound_id'] for r in test_rows]
    metrics = {}
    for name in names:
        threshold = float(lock['models'][name]['validation_threshold'])
        require(math.isfinite(threshold), 'Invalid fixed threshold')
        metrics[name] = {
            'ranking': ranking_metrics(y, scores[name], ids, fraction=.1),
            'validation_selected_threshold': threshold,
            'classification': threshold_metrics(y, scores[name], threshold),
        }
    uncertainty, draws = paired_cluster_bootstrap(y, scores,
        [r['scaffold_group'] for r in test_rows], lock['classical_comparator'])
    report = {
        'status': 'SCORED_LOCKED_PREDICTIONS', 'cohort_sha256': lock['cohort_sha256'],
        'protocol_sha256': lock['protocol_sha256'], 'selection_lock_sha256': sha(lock_path),
        'predictions_sha256': lock['predictions_sha256'],
        'n': len(y), 'positive': int(y.sum()), 'negative': int((1-y).sum()),
        'source_families': 1, 'scaffolds': len({r['scaffold_group'] for r in test_rows}),
        'metrics': metrics, 'uncertainty': uncertainty,
        'limitations': ['One observed patent family, protocol/time/chemistry/source confounded.',
            'Source-defined biochemical label, not cellular efficacy or clinical outcome.',
            'Scores are not calibrated prospective hit probabilities.',
            'Modern structure-pretrained model; historical model availability is not claimed.'],
    }
    return report, draws


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--lock', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path,
                        help='Separate empty destination for deterministic scoring replay')
    args = parser.parse_args()
    if args.output_dir is None:
        report_path = args.root/'reports/results/final_metrics.json'
        draws_path = args.root/'data/predictions/bootstrap_draws.json'
    else:
        report_path = args.output_dir/'final_metrics.json'
        draws_path = args.output_dir/'bootstrap_draws.json'
    require(not report_path.exists() and not draws_path.exists(),
            'Scoring output exists; preserve it and use --output-dir for replay')
    report, draws = score_frozen(args.root, args.lock)
    for path, value in ((draws_path, draws), (report_path, report)):
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open('x') as handle:
            handle.write(json.dumps(value, indent=2, allow_nan=False)+'\n')
    print(json.dumps({'status':report['status'],'n':report['n'],'models':list(report['metrics'])}))


if __name__ == '__main__':
    main()
