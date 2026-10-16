import unittest
import copy
from selection_lock import validate_lock_definition, MODEL_NAMES, CODE_PATHS
from benchmark_common import PROTOCOL_SHA256, COHORT_SHA256
import numpy as np
from evaluate_locked import paired_cluster_bootstrap, validate_prediction_rows

class ClusterUncertainty(unittest.TestCase):
    def test_paired_identical_predictions_have_zero_difference(self):
        y=[0,1,0,1,0,1];score=[.1,.9,.2,.8,.3,.7]
        report,draws=paired_cluster_bootstrap(y,{'a':score,'b':score},['x','x','y','y','z','z'],'a',draws=30)
        self.assertEqual(report['models']['b']['paired_ap_difference_vs_comparator_percentile_95'],[0.,0.])
        self.assertEqual(report['valid_draws'],30)
        self.assertEqual(draws['a'],draws['b'])
    def test_invalid_single_class_draws_are_counted(self):
        report,_=paired_cluster_bootstrap([0,0,1,1],{'a':[.1,.2,.8,.9]},['n','n','p','p'],'a',draws=100,seed=2)
        self.assertGreater(report['invalid_single_class_draws'],0)
        self.assertEqual(report['valid_draws']+report['invalid_single_class_draws'],100)
    def test_no_valid_bootstrap_fails(self):
        with self.assertRaisesRegex(ValueError,'No valid'):
            paired_cluster_bootstrap([0,0],{'a':[.1,.2]},['a','b'],'a',draws=3)

class PredictionCoverage(unittest.TestCase):
    def test_order_and_exact_coverage(self):
        rows=[{'model':'a','compound_id':'2','split':'test','score':'.8'},
              {'model':'a','compound_id':'1','split':'test','score':'.1'}]
        reference=[{'compound_id':'1'},{'compound_id':'2'}]
        result=validate_prediction_rows(rows,reference,['a'])
        np.testing.assert_array_equal(result['a'],[.1,.8])
        with self.assertRaisesRegex(ValueError,'Duplicate prediction'):
            validate_prediction_rows(rows+rows[:1],reference,['a'])
        with self.assertRaisesRegex(ValueError,'Unequal prediction coverage'):
            validate_prediction_rows(rows[:1],reference,['a'])
    def test_non_test_and_nan_rejected(self):
        row={'model':'a','compound_id':'1','split':'validation','score':'.5'}
        with self.assertRaisesRegex(ValueError,'test predictions only'):
            validate_prediction_rows([row],[{'compound_id':'1'}],['a'])
        row.update(split='test',score='nan')
        with self.assertRaisesRegex(ValueError,'Invalid classifier score'):
            validate_prediction_rows([row],[{'compound_id':'1'}],['a'])


class StrictSelectionLock(unittest.TestCase):
    def setUp(self):
        self.expected={'models':{name:{'validation_threshold':.5} for name in MODEL_NAMES},
                       'classical_comparator':'random_forest',
                       'development_artifacts':{'mandatory_sealed_report':'sha'}}
        self.code={name:'digest' for name in CODE_PATHS}
        self.lock={'status':'LOCKED_BEFORE_FINAL_TEST_SCORING',
                   'protocol_sha256':PROTOCOL_SHA256,'cohort_sha256':COHORT_SHA256,
                   **copy.deepcopy(self.expected),'code_sha256':self.code.copy()}
    def test_exact_development_choices_pass(self):
        validate_lock_definition(self.lock,self.expected,self.code)
    def test_threshold_comparator_or_model_omission_fails(self):
        changes=[lambda x:x['models']['random_forest'].update(validation_threshold=.2),
                 lambda x:x.update(classical_comparator='constant_prevalence'),
                 lambda x:x['models'].pop('adapted_encoder')]
        for change in changes:
            altered=copy.deepcopy(self.lock);change(altered)
            with self.assertRaises(ValueError):validate_lock_definition(altered,self.expected,self.code)
    def test_empty_evidence_alternate_freeze_and_unbound_code_fail(self):
        changes=[lambda x:x.update(development_artifacts={}),
                 lambda x:x.update(cohort_sha256='arbitrary-new-file-hash'),
                 lambda x:x.update(code_sha256={}),
                 lambda x:x['code_sha256'].update({'scripts/evaluate_locked.py':'changed'})]
        for change in changes:
            altered=copy.deepcopy(self.lock);change(altered)
            with self.assertRaises(ValueError):validate_lock_definition(altered,self.expected,self.code)

if __name__=='__main__':unittest.main()
