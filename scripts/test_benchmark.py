"""Small, synthetic contract tests; never fit or score the real cohort."""
import csv
import hashlib
import json
import tempfile
import unittest
from pathlib import Path
from unittest import mock

import numpy as np
import benchmark_common as common
import train_baselines as training


class RankingContracts(unittest.TestCase):
    def test_constant_scores_ap_and_random_topk(self):
        y = [1, 0, 1, 0, 0]
        result = common.ranking_metrics(y, [.4]*5, list('edcba'), fraction=.4)
        self.assertAlmostEqual(result['average_precision'], .4)
        self.assertAlmostEqual(result['roc_auc'], .5)
        top = result['screening']
        self.assertEqual(top['k'], 2)
        self.assertAlmostEqual(top['expected_hits'], .8)
        self.assertEqual((top['minimum_hits'], top['maximum_hits']), (0, 2))
        self.assertEqual(top['deterministic_ids'], ['a', 'b'])
        self.assertEqual(top['deterministic_hits'], 0)
        self.assertTrue(top['arbitrary_order_within_boundary_tie'])

    def test_boundary_tie_includes_expected_and_extremes(self):
        result = common.ranking_metrics([1, 1, 0, 0, 0], [1, .8, .8, .8, .1], list('abcde'), .4)
        top = result['screening']
        self.assertAlmostEqual(top['expected_hits'], 1 + 1/3)
        self.assertEqual((top['minimum_hits'], top['maximum_hits']), (1, 2))
        self.assertEqual(top['strictly_above_count'], 1)
        self.assertEqual(top['boundary_tie_count'], 3)
        self.assertEqual(top['slots_in_boundary_tie'], 1)

    def test_single_class_and_nonfinite_inputs(self):
        result = common.ranking_metrics([0, 0], [.1, .2], ['a', 'b'])
        self.assertEqual(result['average_precision'], 0)
        self.assertIsNone(result['roc_auc'])
        self.assertEqual(result['screening']['expected_hits'], 0)
        self.assertEqual(common.ranking_metrics([1, 1], [0, 0], ['a', 'b'])['average_precision'], 1)
        for bad in ([0, float('nan')], [float('inf'), 1]):
            with self.assertRaises(ValueError): common.ranking_metrics([0, 1], bad, ['a', 'b'])
        with self.assertRaises(ValueError): common.ranking_metrics([0, 1], [0, 1], ['a', 'a'])
        with self.assertRaises(ValueError): common.ranking_metrics([], [], [])

    def test_threshold_ties_prefer_higher_and_confusion_is_integer(self):
        self.assertGreater(common.choose_threshold([0, 1], [.5, .5]), .5)
        result = common.threshold_metrics([1, 1, 0, 0], [.8, .2, .6, .1], .5)
        self.assertEqual({k: result[k] for k in ['tp', 'fn', 'fp', 'tn']}, {'tp':1, 'fn':1, 'fp':1, 'tn':1})
        self.assertEqual(result['balanced_accuracy'], .5)
        for key in ['tp', 'fn', 'fp', 'tn']: self.assertIsInstance(result[key], int)
        with self.assertRaises(ValueError): common.choose_threshold([1, 1], [.2, .4])
        with self.assertRaises(ValueError): common.choose_threshold([0, 1], [.2, .4], split='test')


class FitIsolation(unittest.TestCase):
    def test_train_scaler_does_not_learn_validation_distribution(self):
        x = np.array([[0., 0.], [2., 4.], [4., 8.], [6., 12.]])
        model = common.fit_scaled_logistic(x, [0, 0, 1, 1], C=1, split='train')
        before = model.named_steps['scaler'].mean_.copy()
        model.predict_proba(np.array([[100000., 200000.]]))
        np.testing.assert_array_equal(before, [3., 6.])
        np.testing.assert_array_equal(model.named_steps['scaler'].mean_, before)
        with self.assertRaises(ValueError): common.fit_scaled_logistic(x, [0, 0, 1, 1], C=1, split='validation')

    def test_selection_uses_only_validation_with_frozen_tie_rules(self):
        predictions = [[.5, .5], [.5, .5], [.5, .5]]
        pick = common.select_validation_candidate([0, 1], predictions, [{'C':10}, {'C':.1}, {'C':1}], family='logistic')
        self.assertEqual(pick, 1)
        pick = common.select_validation_candidate([0, 1], predictions[:2], [{'min_samples_leaf':2}, {'min_samples_leaf':5}], family='forest')
        self.assertEqual(pick, 1)
        with self.assertRaises(ValueError):
            common.select_validation_candidate([0, 1], predictions[:1], [{'C':1}], family='logistic', split='test')
        with self.assertRaises(ValueError):
            common.select_validation_candidate([0, 1], predictions[:1], [{'C':100}], family='logistic')

    def test_five_neighbors_include_all_boundary_ties(self):
        # All six training fingerprints equally similar; all six must contribute.
        model = common.TieAwareKNN(k=5).fit(np.eye(6, dtype=np.uint8), [1, 0, 0, 0, 0, 0])
        result = model.predict_proba(np.ones((1, 6), dtype=np.uint8))[0, 1]
        self.assertAlmostEqual(result, 1/6)
        # Packed uint8 dot products must not overflow with 2048 set bits.
        model = common.TieAwareKNN(k=1).fit(np.vstack([np.ones(2048), np.zeros(2048)]), [1, 0])
        self.assertEqual(model.predict_proba(np.ones((1, 2048)))[0, 1], 1)

    def test_fingerprint_policy_is_matched_and_chiral_sensitivity_differs(self):
        smiles = ['C[C@H](O)C(=O)O', 'C[C@@H](O)C(=O)O']
        achiral = common.features(smiles)
        chiral = common.features(smiles, chirality=True)
        self.assertEqual(achiral.shape, (2, 2048))
        np.testing.assert_array_equal(achiral[0], achiral[1])
        self.assertFalse(np.array_equal(chiral[0], chiral[1]))

    def test_unknown_split_fails_and_test_label_is_never_parsed(self):
        rows = [{'compound_id':'a', 'split':'train', 'binary_y':'0'},
                {'compound_id':'b', 'split':'validation', 'binary_y':'1'},
                {'compound_id':'c', 'split':'test', 'binary_y':'NOT_A_LABEL'}]
        dev = common.development_rows(rows)
        self.assertEqual(set(dev), {'train', 'validation'})
        self.assertEqual([r['compound_id'] for r in dev['train']], ['a'])
        self.assertFalse(any(r['compound_id'] == 'c' for group in dev.values() for r in group))
        rows[-1]['split'] = 'future_maybe'
        with self.assertRaises(ValueError): common.development_rows(rows)

    def test_protocol_cohort_and_preflight_hash_guards(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root/'docs').mkdir(); (root/'data/processed').mkdir(parents=True)
            (root/'reports/curation').mkdir(parents=True)
            (root/'docs/PROTOCOL.md').write_text('tampered protocol')
            with self.assertRaisesRegex(ValueError, 'protocol'):
                common.load_frozen_development(root)
            protocol = b'fixed'; cohort = b'compound_id,split,binary_y\na,train,0\nb,validation,1\nc,test,UNREAD\n'
            (root/'docs/PROTOCOL.md').write_bytes(protocol)
            (root/'data/processed/cohort.csv').write_bytes(cohort)
            gate = {'status':'PASS_FOR_NARROW_SOURCE_SHIFT_BINARY_BENCHMARK', 'outputs': {'data/processed/cohort.csv':hashlib.sha256(cohort).hexdigest()}}
            (root/'reports/curation/cohort_gate.json').write_text(json.dumps(gate))
            with mock.patch.object(common, 'PROTOCOL_SHA256', hashlib.sha256(protocol).hexdigest()), mock.patch.object(common, 'COHORT_SHA256', hashlib.sha256(cohort).hexdigest()), mock.patch.object(common, 'EXPECTED_COUNTS', {'train':1,'validation':1,'test':1}):
                dev, bindings = common.load_frozen_development(root)
                self.assertEqual(len(dev['train']), 1)
                (root/'data/processed/cohort.csv').write_bytes(cohort+b'corrupt')
                with self.assertRaisesRegex(ValueError, 'cohort'): common.load_frozen_development(root)
            preflight = root/'preflight.json'
            preflight.write_text(json.dumps({'status':'PASS','protocol_sha256':common.PROTOCOL_SHA256,'cohort_sha256':'wrong'}))
            with self.assertRaisesRegex(ValueError, 'cohort'): common.require_preflight(preflight)
            preflight.write_text(json.dumps({'status':'PENDING','protocol_sha256':common.PROTOCOL_SHA256,'cohort_sha256':common.COHORT_SHA256}))
            with self.assertRaisesRegex(ValueError, 'PASS'): common.require_preflight(preflight)

    def test_cpu_preflight_requires_input_pass_and_revision_not_gpu_pass(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'preflight.json'
            receipt = {'status':'INPUT_PASS_RUNTIME_FAILED', 'protocol_sha256':common.PROTOCOL_SHA256,
                       'cohort_sha256':common.COHORT_SHA256, 'revision':'361063d0ad524ef77cf39b08469f6be770dc550f',
                       'input_gate':{'status':'PASS'}}
            path.write_text(json.dumps(receipt))
            self.assertIn('encoder_input_gate_sha256', common.require_preflight(path))
            receipt['input_gate']['status'] = 'FAIL'
            path.write_text(json.dumps(receipt))
            with self.assertRaisesRegex(ValueError, 'input'):
                common.require_preflight(path)
            receipt['input_gate']['status'] = 'PASS'
            receipt['revision'] = 'unreviewed_revision'
            path.write_text(json.dumps(receipt))
            with self.assertRaisesRegex(ValueError, 'revision'):
                common.require_preflight(path)

    def test_baseline_budget_is_exact(self):
        specs = training.candidate_specs()
        self.assertEqual(len(specs), 10)
        self.assertEqual([s['params']['C'] for s in specs if s['family']=='fingerprint_logistic'], [.1,1,10])
        self.assertEqual([s['params']['C'] for s in specs if s['family']=='fingerprint_logistic_chiral'], [.1,1,10])
        self.assertEqual([s['params']['min_samples_leaf'] for s in specs if s['family']=='random_forest'], [2,5])
        for spec in specs:
            model = training.make_model(spec)
            if spec['family'] == 'random_forest':
                self.assertEqual(model.n_estimators, 500)
                self.assertEqual(model.max_features, 'sqrt')
                self.assertEqual(model.random_state, 17)
                self.assertEqual(model.class_weight, 'balanced')
            if 'logistic' in spec['family']:
                self.assertEqual(model.max_iter, 5000)
                self.assertEqual(model.class_weight, 'balanced')

    def test_training_entry_checks_preflight_before_fit(self):
        with tempfile.TemporaryDirectory() as directory, mock.patch.object(training, 'fit_candidates') as fit:
            with self.assertRaises((FileNotFoundError, ValueError)):
                training.run(Path(directory), Path(directory)/'missing_preflight.json')
            fit.assert_not_called()


class OutputIntegrity(unittest.TestCase):
    def test_public_report_omits_local_executable_without_mutating_receipt(self):
        report = {'runtime':{'executable':'/private/environment/python','python':'3.12'},
                  'test_scored':False,'candidates':[{'candidate_id':'fixed'}]}
        published = training.public_report(report)
        self.assertNotIn('executable', published['runtime'])
        self.assertEqual(report['runtime']['executable'], '/private/environment/python')
        self.assertEqual(published['candidates'], report['candidates'])

    def test_repeated_cpu_fit_has_identical_predictions(self):
        x = np.array([[0.,0.],[1.,0.],[0.,1.],[1.,1.],[2.,0.],[0.,2.]])
        model, score, difference = training.fit_twice(common.logistic_estimator(1), x, np.array([0,0,0,1,1,1]), x[:2])
        self.assertEqual(difference, 0)
        self.assertTrue(np.isfinite(score).all())

    def test_unknown_or_test_group_never_reaches_feature_creation(self):
        with mock.patch.object(training, 'features') as feature:
            with self.assertRaises(ValueError):
                training.fit_candidates({'train':[], 'validation':[], 'test':[]}, Path('/unused'))
            feature.assert_not_called()

    def test_failed_publication_invalidates_seal(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            bindings = {'protocol_sha256':common.PROTOCOL_SHA256,'cohort_sha256':common.COHORT_SHA256,'cohort_gate_sha256':'x'}
            original_write = training.write_json
            def write_or_fail(path, value):
                if Path(path) == root/'reports/models/validation_baselines.json':
                    raise OSError('injected final publication failure')
                return original_write(path, value)
            with mock.patch.object(training, 'require_preflight', return_value={'encoder_preflight_sha256':'p'}), \
                 mock.patch.object(training, 'load_frozen_development', return_value=({'train':[],'validation':[]},bindings)), \
                 mock.patch.object(training, 'fit_candidates', return_value={'test_scored':False}), \
                 mock.patch.object(training, 'write_json', side_effect=write_or_fail):
                with self.assertRaisesRegex(OSError, 'injected'):
                    training.run(root, root/'preflight.json')
            destination = root/'data/models/baselines'
            self.assertTrue((destination/'FAILED.json').exists())
            self.assertFalse((destination/'SEALED.json').exists())
            self.assertTrue((destination/'SEAL_INVALIDATED.json').exists())


if __name__ == '__main__':
    unittest.main()
