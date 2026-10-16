"""Public replay integration and corruption checks; no model weights required."""
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from reproduce_metrics import compare_values, verify_artifacts
from evaluate_locked import paired_cluster_bootstrap, validate_prediction_rows


class PublicReplay(unittest.TestCase):
    def test_release_reproduces_every_metric_and_bootstrap_draw_outside_checkout(self):
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            # A moved checkout has neither private models nor execution receipts.
            release = base / 'relocated-checkout'
            shutil.copytree(ROOT, release, ignore=shutil.ignore_patterns('__pycache__', 'replay-output', '.venv*', '.git', 'outputs', 'receipts'))
            output = base / 'output'
            run = subprocess.run([sys.executable, str(release / 'scripts/reproduce_metrics.py'),
                                  '--output-dir', str(output)], cwd=base, text=True, capture_output=True)
            self.assertEqual(run.returncode, 0, run.stderr)
            verification = json.loads((output / 'verification.json').read_text())
            self.assertTrue(verification['metrics_exact_json_values'])
            self.assertTrue(verification['bootstrap_exact_json_values'])
            self.assertFalse(verification['private_training_and_prediction_receipt_chain_verified'])
            self.assertEqual(verification['bootstrap_draws'], 2000)
            self.assertEqual(len(verification['models']), 7)
            for relative in ('reports/results/final_metrics.json', 'data/predictions/bootstrap_draws.json'):
                self.assertEqual((output / Path(relative).name).read_bytes(), (ROOT / relative).read_bytes())
            rerun = subprocess.run([sys.executable, str(release / 'scripts/reproduce_metrics.py'),
                                    '--output-dir', str(output)], cwd=base, text=True, capture_output=True)
            self.assertNotEqual(rerun.returncode, 0)
            self.assertIn('Output directory exists', rerun.stderr)

    def test_changed_predictions_and_missing_code_fail_before_scoring(self):
        for relative, mutation in [('data/predictions/final_predictions.csv', 'change'),
                                   ('scripts/benchmark_common.py', 'delete')]:
            with self.subTest(relative=relative), tempfile.TemporaryDirectory() as directory:
                release = Path(directory) / 'release'
                shutil.copytree(ROOT, release, ignore=shutil.ignore_patterns('__pycache__', 'replay-output', '.venv*', '.git', 'outputs', 'receipts'))
                path = release / relative
                if mutation == 'change':
                    raw = path.read_bytes()
                    path.write_bytes(raw.replace(b'3.1543611517753562e-05', b'4.1543611517753562e-05', 1))
                else:
                    path.unlink()
                with self.assertRaisesRegex(ValueError, '(?i)public replay artifact'):
                    verify_artifacts(release)

    def test_manifest_change_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            release = Path(directory)
            (release / 'manifests').mkdir()
            shutil.copyfile(ROOT / 'manifests/public_scoring_replay.json',
                            release / 'manifests/public_scoring_replay.json')
            with (release / 'manifests/public_scoring_replay.json').open('a') as handle:
                handle.write(' ')
            with self.assertRaisesRegex(ValueError, 'manifest hash mismatch'):
                verify_artifacts(release)

    def test_partial_prediction_coverage_and_duplicate_id_fail(self):
        reference = [{'compound_id': 'one'}, {'compound_id': 'two'}]
        rows = [{'model': 'a', 'compound_id': 'one', 'split': 'test', 'score': '.1'},
                {'model': 'a', 'compound_id': 'two', 'split': 'test', 'score': '.9'}]
        with self.assertRaisesRegex(ValueError, 'Unequal prediction coverage'):
            validate_prediction_rows(rows[:1], reference, ['a'])
        with self.assertRaisesRegex(ValueError, 'Duplicate prediction'):
            validate_prediction_rows(rows + rows[:1], reference, ['a'])

    def test_result_comparison_rejects_real_differences(self):
        with self.assertRaisesRegex(ValueError, 'Numeric value differs'):
            compare_values({'ap': .2}, {'ap': .3})
        with self.assertRaisesRegex(ValueError, 'Keys differ'):
            compare_values({'ap': .2}, {'ap': .2, 'auc': .5})

    def test_pairing_preserves_zero_difference_for_identical_methods(self):
        report, draws = paired_cluster_bootstrap(
            [0, 1, 0, 1], {'a': [.1, .8, .2, .9], 'b': [.1, .8, .2, .9]},
            ['first', 'first', 'second', 'second'], 'a', draws=20)
        self.assertEqual(draws['a'], draws['b'])
        self.assertEqual(report['models']['b']['paired_ap_difference_vs_comparator_percentile_95'], [0., 0.])


if __name__ == '__main__':
    unittest.main()
