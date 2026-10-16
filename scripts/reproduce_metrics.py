#!/usr/bin/env python3
"""Recompute the released seven-method metrics from saved final predictions.

This replay does not fit models, select thresholds, regenerate predictions, or
validate the private training/receipt chain. It calls the original frozen metric
functions without changing the original source files or their recorded hashes.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np

# Filled from the released byte-level replay manifest, independently of --root.
MANIFEST_SHA256 = 'ed3ed19d5d41b89643453e7bbac2b9a379439333c84c6ca3dc3e91f8effbac47'
MANIFEST_PATH = 'manifests/public_scoring_replay.json'


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def verify_artifacts(root):
    root = Path(root).resolve()
    path = root / MANIFEST_PATH
    require(sha256(path) == MANIFEST_SHA256, 'Replay manifest hash mismatch')
    manifest = json.loads(path.read_text())
    require(manifest['schema_version'] == 1, 'Unsupported replay manifest')
    for relative, record in manifest['files'].items():
        artifact = (root / relative).resolve()
        require(artifact.is_relative_to(root), 'Artifact escapes release root')
        require(artifact.is_file(), 'Missing public replay artifact: ' + relative)
        require(artifact.stat().st_size == record['bytes'] and sha256(artifact) == record['sha256'],
                'Public replay artifact hash mismatch: ' + relative)
    return manifest


def compare_values(actual, expected, path='result'):
    """Check every value, allowing only negligible floating-point runtime drift."""
    if isinstance(expected, dict):
        require(isinstance(actual, dict) and actual.keys() == expected.keys(), 'Keys differ: ' + path)
        for key in expected:
            compare_values(actual[key], expected[key], path + '.' + key)
    elif isinstance(expected, list):
        require(isinstance(actual, list) and len(actual) == len(expected), 'Length differs: ' + path)
        for index, value in enumerate(expected):
            compare_values(actual[index], value, path + '[' + str(index) + ']')
    elif isinstance(expected, float):
        require(isinstance(actual, (int, float)) and math.isclose(actual, expected, rel_tol=1e-12, abs_tol=1e-12),
                'Numeric value differs: ' + path)
    else:
        require(type(actual) is type(expected) and actual == expected, 'Value differs: ' + path)


def reproduce(root):
    root = Path(root).resolve()
    manifest = verify_artifacts(root)
    # Normal imports preserve the exact frozen implementations and their hashes.
    # No encoder runtime is imported: baseline.txt is sufficient, including RDKit.
    import benchmark_common
    import evaluate_locked
    import selection_lock
    for module in (benchmark_common, evaluate_locked, selection_lock):
        relative = 'scripts/' + Path(module.__file__).name
        require(sha256(module.__file__) == manifest['files'][relative]['sha256'],
                'Imported scoring module hash mismatch: ' + relative)
    from benchmark_common import ranking_metrics, threshold_metrics
    from evaluate_locked import load_csv, paired_cluster_bootstrap, validate_prediction_rows
    from selection_lock import MODEL_NAMES

    lock_path = root / 'manifests/final_selection_lock.json'
    lock = json.loads(lock_path.read_text())
    require(lock['status'] == 'LOCKED_BEFORE_FINAL_TEST_SCORING', 'Unexpected original lock status')
    require(set(lock['models']) == MODEL_NAMES, 'Unexpected method set')
    files = manifest['files']
    require(lock['cohort_sha256'] == files['data/processed/cohort.csv']['sha256'], 'Cohort binding mismatch')
    require(lock['predictions_path'] == 'data/predictions/final_predictions.csv' and
            lock['predictions_sha256'] == files[lock['predictions_path']]['sha256'], 'Prediction binding mismatch')
    for relative, digest in lock['code_sha256'].items():
        require(files[relative]['sha256'] == digest, 'Original code binding mismatch: ' + relative)
    rows = [row for row in load_csv(root / 'data/processed/cohort.csv') if row['split'] == 'test']
    require(len(rows) == 246 and len({row['source_family'] for row in rows}) == 1, 'Unexpected test source')
    require(all(row['binary_y'] in ('0', '1') for row in rows), 'Nonbinary test label')
    y = np.array([int(row['binary_y']) for row in rows])
    require(int(y.sum()) == 69, 'Unexpected test class count')
    ids = [row['compound_id'] for row in rows]
    scores = validate_prediction_rows(load_csv(root / lock['predictions_path']), rows, list(lock['models']))
    metrics = {}
    for name in lock['models']:
        threshold = float(lock['models'][name]['validation_threshold'])
        require(math.isfinite(threshold), 'Invalid frozen threshold')
        metrics[name] = {
            'ranking': ranking_metrics(y, scores[name], ids, fraction=.1),
            'validation_selected_threshold': threshold,
            'classification': threshold_metrics(y, scores[name], threshold),
        }
    uncertainty, draws = paired_cluster_bootstrap(
        y, scores, [row['scaffold_group'] for row in rows], lock['classical_comparator'])
    report = {
        'status': 'SCORED_LOCKED_PREDICTIONS', 'cohort_sha256': lock['cohort_sha256'],
        'protocol_sha256': lock['protocol_sha256'], 'selection_lock_sha256': sha256(lock_path),
        'predictions_sha256': lock['predictions_sha256'],
        'n': len(y), 'positive': int(y.sum()), 'negative': int((1-y).sum()),
        'source_families': 1, 'scaffolds': len({row['scaffold_group'] for row in rows}),
        'metrics': metrics, 'uncertainty': uncertainty,
        'limitations': ['One observed patent family, protocol/time/chemistry/source confounded.',
            'Source-defined biochemical label, not cellular efficacy or clinical outcome.',
            'Scores are not calibrated prospective hit probabilities.',
            'Modern structure-pretrained model; historical model availability is not claimed.'],
    }
    expected = json.loads((root / 'reports/results/final_metrics.json').read_text())
    expected_draws = json.loads((root / 'data/predictions/bootstrap_draws.json').read_text())
    compare_values(report, expected)
    compare_values(draws, expected_draws, 'bootstrap_draws')
    verification = {
        'status': 'PUBLIC_SAVED_PREDICTION_REPLAY_VERIFIED',
        'scope': 'Released predictions, original thresholds, seven-method scoring and paired scaffold bootstrap only',
        'manifest_sha256': MANIFEST_SHA256,
        'original_protocol_freeze_sha256': lock['protocol_sha256'],
        'protocol_note': 'Original freeze identity retained; editorial publication of protocol is a separate binding',
        'private_training_and_prediction_receipt_chain_verified': False,
        'raw_source_reconstruction_performed': False,
        'model_fitting_or_selection_performed': False,
        'all_metric_and_draw_values_matched': True,
        'float_relative_tolerance': 1e-12, 'float_absolute_tolerance': 1e-12,
        'metrics_exact_json_values': report == expected,
        'bootstrap_exact_json_values': draws == expected_draws,
        'n': len(y), 'models': list(metrics), 'bootstrap_draws': uncertainty['valid_draws'],
    }
    return report, draws, verification


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path(__file__).resolve().parents[1],
                        help='Public release root (defaults to this checkout)')
    parser.add_argument('--output-dir', type=Path, required=True,
                        help='New destination directory; existing directories are never overwritten')
    args = parser.parse_args()
    require(not args.output_dir.exists(), 'Output directory exists; choose a new replay destination')
    report, draws, verification = reproduce(args.root)
    args.output_dir.mkdir(parents=True, exist_ok=False)
    for filename, value in (('final_metrics.json', report), ('bootstrap_draws.json', draws),
                            ('verification.json', verification)):
        (args.output_dir / filename).write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    print(json.dumps(verification, allow_nan=False))


if __name__ == '__main__':
    main()
