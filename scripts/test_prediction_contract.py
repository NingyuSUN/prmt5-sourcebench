#!/usr/bin/env python3
"""Synthetic prediction-contract tests; no real test scores or model fitting."""
import ast
import hashlib
import json
from pathlib import Path
import tempfile
import unittest
import sys
from types import SimpleNamespace
from unittest.mock import patch
import numpy as np
import predict_locked as prediction

class LabelTrap(dict):
    def __getitem__(self, key):
        if key not in ('split','compound_id','canonical_smiles'):
            raise AssertionError('Non-input field accessed: '+key)
        return super().__getitem__(key)
    def get(self, key, *args):
        if key not in ('split','compound_id','canonical_smiles'):
            raise AssertionError('Non-input field accessed: '+key)
        return super().get(key,*args)
    def items(self):
        raise AssertionError('All columns copied')

class PredictionContract(unittest.TestCase):
    def test_projection_never_accesses_target_or_copies_extra_fields(self):
        rows=[LabelTrap(split='train',compound_id='tr',canonical_smiles='C',binary_y=object()),
              LabelTrap(split='validation',compound_id='va',canonical_smiles='CC',binary_y=object()),
              LabelTrap(split='test',compound_id='te',canonical_smiles='CO',binary_y=object())]
        projected=prediction.test_inputs(rows,{'train':1,'validation':1,'test':1})
        self.assertEqual(projected,[{'split':'test','compound_id':'te','canonical_smiles':'CO'}])
    def test_wrong_split_count_duplicate_id_and_unknown_split_block(self):
        row={'split':'test','compound_id':'a','canonical_smiles':'C'}
        for rows,expected in [([row],{'test':2}),([row,row],None),([dict(row,split='new')],None)]:
            with self.assertRaises(ValueError):prediction.test_inputs(rows,expected)
    def test_locked_hash_checked_before_deserialization(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'model.joblib';path.write_bytes(b'not-a-model')
            with self.assertRaisesRegex(ValueError,'changed'):
                prediction.verified_artifact(directory,'model.joblib',{'model.joblib':'0'*64})
            with self.assertRaisesRegex(ValueError,'absent'):
                prediction.verified_artifact(directory,'model.joblib',{})
            with self.assertRaises(ValueError):prediction.inside(directory,'../outside')
    def test_environment_mismatch_blocks(self):
        with self.assertRaisesRegex(ValueError,'version differs'):
            prediction.require_runtime({'scikit_learn':'1.7.2'},{'scikit_learn':'1.9.0'},['scikit_learn'])
    def test_all_seeds_mean_not_best_seed(self):
        scores={17:[.2,.1],29:[.8,.4],43:[.5,.7]}
        np.testing.assert_allclose(prediction.mean_seed_scores(scores,2),[.5,.4])
        with self.assertRaises(ValueError):prediction.mean_seed_scores({17:[.2,.1],29:[.8,.4]},2)
    def test_score_shape_range_finite_and_model_coverage(self):
        inputs=[{'compound_id':'a'},{'compound_id':'b'}]
        for scores in [[float('nan'),.1],[-.1,.1],[1.1,.1],[.1],[[.1],[.2]]]:
            with self.assertRaises(ValueError):prediction.prediction_rows(inputs,{'m':scores},{'m'})
        with self.assertRaises(ValueError):prediction.prediction_rows(inputs,{'m':[.1,.2]},{'other'})
        rows=prediction.prediction_rows(inputs,{'m':[.1,.2]},{'m'})
        self.assertEqual(set(rows[0]),{'model','compound_id','split','score'})
        self.assertEqual([row['compound_id'] for row in rows],['a','b'])
    def test_atomic_publication_failure_and_no_overwrite(self):
        with tempfile.TemporaryDirectory() as directory:
            destination=Path(directory)/'predictions'
            rows=[{'model':'m','compound_id':'a','split':'test','score':'.5'}]
            with patch('predict_locked.os.rename',side_effect=OSError('simulated storage failure')):
                with self.assertRaises(OSError):prediction.publish_output(destination,rows,{'status':'ready'})
            self.assertFalse(destination.exists())
            self.assertEqual(list(Path(directory).iterdir()),[])
            receipt=prediction.publish_output(destination,rows,{'status':'ready'})
            self.assertEqual(receipt['predictions_sha256'],prediction.sha256(destination/'predictions.csv'))
            with self.assertRaises(ValueError):prediction.publish_output(destination,rows,{})
    def test_end_to_end_synthetic_publication_requires_lock_and_binds_code(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);(root/'manifests').mkdir();(root/'docs').mkdir()
            (root/'data/processed').mkdir(parents=True)
            (root/'manifests/model_selection_lock.json').write_text('{}')
            (root/'docs/PROTOCOL.md').write_text('synthetic protocol')
            (root/'data/processed/cohort.csv').write_text('synthetic bytes')
            ph=prediction.sha256(root/'docs/PROTOCOL.md');ch=prediction.sha256(root/'data/processed/cohort.csv')
            lock={'models':{name:{} for name in prediction.BASELINE_NAMES},
                  'code_sha256':{},'development_artifacts':{}}
            inputs=[{'compound_id':'synthetic_a','canonical_smiles':'C','split':'test'}]
            outputs={name:[.5] for name in prediction.BASELINE_NAMES}
            module=SimpleNamespace(load_locked_model_selection=lambda unused:lock)
            with patch.dict(sys.modules,{'selection_lock':module}), \
                 patch.object(prediction,'PROTOCOL_SHA256',ph),patch.object(prediction,'COHORT_SHA256',ch), \
                 patch.object(prediction,'load_test_inputs',return_value=inputs), \
                 patch.object(prediction,'baseline_predictions',return_value=(outputs,{})):
                receipt=prediction.run(root,'baselines')
            self.assertFalse(receipt['target_columns_accessed'])
            self.assertEqual(receipt['prediction_row_count'],5)
            self.assertEqual(receipt['model_selection_lock_sha256'],prediction.sha256(root/'manifests/model_selection_lock.json'))
            self.assertTrue((root/'data/predictions/baselines/predictions.csv').exists())
    def test_missing_selection_lock_blocks_before_predictions(self):
        with tempfile.TemporaryDirectory() as directory:
            module=SimpleNamespace(load_locked_model_selection=lambda unused:(_ for _ in ()).throw(ValueError('unlocked')))
            with patch.dict(sys.modules,{'selection_lock':module}), \
                 patch.object(prediction,'baseline_predictions') as infer:
                with self.assertRaisesRegex(ValueError,'unlocked'):prediction.run(directory,'baselines')
                infer.assert_not_called()
            self.assertFalse((Path(directory)/'data/predictions').exists())
    def test_input_bytes_changed_during_inference_prevent_publication(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);(root/'manifests').mkdir();(root/'docs').mkdir()
            (root/'data/processed').mkdir(parents=True)
            (root/'manifests/model_selection_lock.json').write_text('{}')
            protocol=root/'docs/PROTOCOL.md';protocol.write_text('synthetic protocol')
            cohort=root/'data/processed/cohort.csv';cohort.write_text('synthetic bytes')
            lock={'models':{name:{} for name in prediction.BASELINE_NAMES},'code_sha256':{},'development_artifacts':{}}
            def mutate(*unused):
                protocol.write_text('changed')
                return {name:[.5] for name in prediction.BASELINE_NAMES},{}
            module=SimpleNamespace(load_locked_model_selection=lambda unused:lock)
            with patch.dict(sys.modules,{'selection_lock':module}), \
                 patch.object(prediction,'PROTOCOL_SHA256',prediction.sha256(protocol)), \
                 patch.object(prediction,'COHORT_SHA256',prediction.sha256(cohort)), \
                 patch.object(prediction,'load_test_inputs',return_value=[{'compound_id':'synthetic_a'}]), \
                 patch.object(prediction,'baseline_predictions',side_effect=mutate):
                with self.assertRaisesRegex(ValueError,'changed during prediction'):prediction.run(root,'baselines')
            self.assertFalse((root/'data/predictions').exists())
    def test_no_fitting_or_metric_calls_in_prediction_entry(self):
        tree=ast.parse(Path(prediction.__file__).read_text())
        forbidden={'fit','fit_transform','partial_fit','train_epoch','make_optimizer',
                   'ranking_metrics','threshold_metrics','average_precision_score','choose_threshold'}
        for node in ast.walk(tree):
            if isinstance(node,ast.Call):
                name=node.func.attr if isinstance(node.func,ast.Attribute) else getattr(node.func,'id','')
                self.assertNotIn(name,forbidden)

if __name__=='__main__':unittest.main()
