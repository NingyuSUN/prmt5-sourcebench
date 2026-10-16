"""Synthetic post-score checks; never load actual test predictions or labels."""
from pathlib import Path
import copy
import csv
import json
import sys
import tempfile
import unittest
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'scripts'))
import analyze_biology as biology


class BiologicalAnalysis(unittest.TestCase):
    def test_similarity_boundaries_and_reject_invalid(self):
        self.assertEqual([biology.similarity_band(x) for x in [0,.4999,.5,.6999,.7,1]],
                         ['<0.5','<0.5','0.5–<0.7','0.5–<0.7','>=0.7','>=0.7'])
        for value in [-.1,1.1,float('nan')]:
            with self.assertRaises(ValueError): biology.similarity_band(value)

    def test_single_class_and_empty_have_no_ap(self):
        for y,score in [([],[]),([0,0],[.2,.8]),([1],[.9])]:
            result=biology.subgroup_metrics(y,score,.5)
            self.assertIsNone(result['average_precision'])
            self.assertEqual(sum(result[k] for k in ['tp','fp','tn','fn']),len(y))
        result=biology.subgroup_metrics([1,0],[.5,.5],.5)
        self.assertEqual((result['tp'],result['fp'],result['fn'],result['tn']),(1,1,0,0))
        self.assertEqual(result['average_precision'],.5)

    def test_source_joint_bins_preserved_and_conflicts_rejected(self):
        evidence=json.dumps([{'values':{'a':'++','b':'+'}}])
        self.assertEqual(biology.source_bins(evidence,0),'+|++')
        with self.assertRaises(ValueError): biology.source_bins(evidence,1)
        with self.assertRaises(ValueError): biology.source_bins('[{"values":{"a":"+++","b":"+"}}]',0)
        with self.assertRaises(ValueError): biology.source_bins('[{"values":{"a":"unexpected"}}]',0)

    def test_nonisomeric_collision_and_descriptors(self):
        rows=[{'compound_id':'a','canonical_smiles':'N[C@@H](C)O','binary_y':'0'},
              {'compound_id':'b','canonical_smiles':'N[C@H](C)O','binary_y':'1'},
              {'compound_id':'c','canonical_smiles':'CCO','binary_y':'0'}]
        flags,summary=biology.collision_audit(rows)
        self.assertEqual(summary['mixed_label_groups'],1)
        self.assertEqual(summary['empirical_minimum_errors'],1)
        self.assertTrue(flags['a']['nonisomeric_mixed_label_collision'])
        self.assertFalse(flags['c']['nonisomeric_collision'])
        desc=biology.molecular_descriptors('CCO')
        self.assertAlmostEqual(desc['mw'],46.069,places=2)
        self.assertEqual((desc['hbd'],desc['hba']),(1,1))

    def test_all_seed_validation_variation_separate_from_bootstrap(self):
        report={'validation_parents':70,'partial_fine_tuning':{'selected_learning_rate':1e-5},
                'trials':[{'learning_rate':lr,'seed':seed,'best_epoch':2,
                           'best_validation_average_precision':ap}
                          for lr in [1e-5,3e-5] for seed,ap in [(17,.4),(29,.5),(43,.6)]]}
        variation=biology.seed_variation(report)
        self.assertEqual(len(variation['learning_rates']),2)
        self.assertAlmostEqual(variation['learning_rates'][0]['sample_standard_deviation'],.1)
        self.assertIn('validation',variation['metric'])
        report['trials'].pop()
        with self.assertRaises(ValueError): biology.seed_variation(report)

    def test_report_bindings_and_fixed_thresholds(self):
        lock={'cohort_sha256':'c','protocol_sha256':'p','predictions_sha256':'d',
              'models':{name:{'validation_threshold':.5} for name in biology.MODEL_NAMES}}
        final={'status':'SCORED_LOCKED_PREDICTIONS','selection_lock_sha256':'l',
               'cohort_sha256':'c','protocol_sha256':'p','predictions_sha256':'d',
               'n':246,'positive':69,'negative':177,'source_families':1,
               'metrics':{name:{'validation_selected_threshold':.5,'classification':{'threshold':.5}}
                          for name in biology.MODEL_NAMES}}
        biology.validate_final_bindings(final,lock,'l')
        for key in ['predictions_sha256','selection_lock_sha256']:
            bad=copy.deepcopy(final);bad[key]='tampered'
            with self.assertRaises(ValueError): biology.validate_final_bindings(bad,lock,'l')
        bad=copy.deepcopy(final)
        bad['metrics']['frozen_encoder']['validation_selected_threshold']=.4
        with self.assertRaises(ValueError): biology.validate_final_bindings(bad,lock,'l')

    def test_gate_failure_happens_before_cohort_csv_read_or_scoring(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            (root/'lock.json').write_text('{}')
            (root/'metrics.json').write_text('{}')
            with mock.patch.object(biology,'load_csv') as read_csv, mock.patch.object(biology,'subgroup_metrics') as score:
                with self.assertRaises(ValueError):
                    biology.analyze(root,root/'lock.json',root/'metrics.json')
                read_csv.assert_not_called();score.assert_not_called()

    def test_full_selection_gate_failure_precedes_csv_read(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            (root/'lock.json').write_text('{}');(root/'metrics.json').write_text('{}')
            with mock.patch.object(biology,'validate_final_bindings'), mock.patch.object(
                    biology,'validate_selection_lock',side_effect=ValueError('tampered lock')), mock.patch.object(biology,'load_csv') as read_csv:
                with self.assertRaisesRegex(ValueError,'tampered lock'):
                    biology.analyze(root,root/'lock.json',root/'metrics.json')
                read_csv.assert_not_called()

    def test_test_seed_probabilities_reordered_without_selection(self):
        rows=[{'compound_id':'a','binary_y':'1'},{'compound_id':'b','binary_y':'0'}]
        receipt={'compound_ids':['b','a'],'seed_probabilities':{'17':[.2,.8],'29':[.9,.1],'43':[.4,.6]}}
        result=biology.test_seed_variation(receipt,rows,[.5,.5],1e-5)
        self.assertEqual([r['seed'] for r in result['seeds']],[17,29,43])
        self.assertEqual([r['test_average_precision'] for r in result['seeds']],[1,.5,1])
        self.assertEqual(result['ensemble_average_precision'],.5)
        with self.assertRaises(ValueError): biology.test_seed_variation(receipt,rows,[.9,.1],1e-5)

    def test_neighbors_only_pinned_training_reference(self):
        rows=[{'reference_split':'train','query_split':'test','query_id':'a','nearest_id':'t','max_tanimoto':'.5'},
              {'reference_split':'validation','query_split':'test','query_id':'a','nearest_id':'v','max_tanimoto':'.9'}]
        found=biology.training_neighbors(rows,{'a'},{'t'})
        self.assertEqual(found['a']['nearest_train_tanimoto'],.5)
        with self.assertRaises(ValueError): biology.training_neighbors(rows+rows[:1],{'a'},{'t'})
        with self.assertRaises(ValueError): biology.training_neighbors(rows,{'a','b'},{'t'})

    def test_complete_synthetic_postscore_analysis_without_real_inputs(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            def write_json(relative,value):
                path=root/relative;path.parent.mkdir(parents=True,exist_ok=True)
                path.write_text(json.dumps(value));return path
            def write_csv(relative,rows):
                path=root/relative;path.parent.mkdir(parents=True,exist_ok=True)
                with path.open('w',newline='') as handle:
                    writer=csv.DictWriter(handle,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
                return path
            test=[{'compound_id':f'c{i}','canonical_smiles':'CCO','binary_y':int(i<69),
                   'split':'test','source_family':'synthetic','scaffold_group':'synthetic',
                   'source_evidence':json.dumps([{'values':{'example':'+++' if i<69 else '+'}}])}
                  for i in range(246)]
            train=dict(test[0],compound_id='train',split='train')
            cohort=write_csv('data/processed/cohort.csv',[train]+test)
            neighbors=write_csv('data/processed/nearest_neighbors.csv',[
                {'reference_split':'train','query_split':'test','query_id':r['compound_id'],
                 'nearest_id':'train','max_tanimoto':.5} for r in test])
            gate=write_json('reports/curation/cohort_gate.json',{
                'fingerprint':{'radius':2,'bits':2048,'include_chirality':True},
                'outputs':{'data/processed/cohort.csv':biology.sha(cohort),
                           'data/processed/nearest_neighbors.csv':biology.sha(neighbors)}})
            _,collision_summary=biology.collision_audit(test)
            preflight=write_json('reports/models/encoder_preflight.json',{'status':'PASS','input_gate':{
                'status':'PASS','rdkit_version':biology.rdBase.rdkitVersion,'by_split':{'test':collision_summary}}})
            encoder={'validation_parents':70,'bindings':{'cohort_gate_sha256':biology.sha(gate)},
                'partial_fine_tuning':{'selected_learning_rate':1e-5},
                'trials':[{'learning_rate':lr,'seed':seed,'best_epoch':1,'best_validation_average_precision':.5}
                          for lr in (1e-5,3e-5) for seed in (17,29,43)]}
            encoder_path=write_json('data/models/encoder/encoder_development.json',encoder)
            models={name:{'validation_threshold':.5} for name in biology.MODEL_NAMES}
            models['adapted_encoder']['learning_rate']=1e-5
            predictions=write_csv('sealed.csv',[{'model':name,'compound_id':row['compound_id'],'split':'test','score':.5}
                                               for name in biology.MODEL_NAMES for row in test])
            receipt=write_json('data/predictions/encoder/receipt.json',{
                'status':'LOCKED_TEST_PREDICTIONS_NO_METRICS','kind':'encoder','protocol_sha256':'synthetic',
                'cohort_sha256':biology.sha(cohort),'model_selection_lock_sha256':'synthetic',
                'models':{k:v for k,v in models.items() if k in ('frozen_encoder','adapted_encoder')},
                'compound_ids':[r['compound_id'] for r in test],
                'seed_probabilities':{str(seed):[.5]*246 for seed in (17,29,43)}})
            lock={'models':models,'protocol_sha256':'synthetic','cohort_sha256':biology.sha(cohort),
                'model_selection_lock_sha256':'synthetic','predictions_path':'sealed.csv',
                'predictions_sha256':biology.sha(predictions),'prediction_receipts':{
                    'data/predictions/encoder/receipt.json':biology.sha(receipt)},
                'development_artifacts':{'data/models/encoder/encoder_development.json':biology.sha(encoder_path),
                                        'reports/models/encoder_preflight.json':biology.sha(preflight)}}
            lock_path=write_json('lock.json',lock)
            classification=biology.subgroup_metrics([r['binary_y'] for r in test],[.5]*246,.5)
            final_path=write_json('metrics.json',{'status':'SCORED_LOCKED_PREDICTIONS',
                'selection_lock_sha256':biology.sha(lock_path),'cohort_sha256':biology.sha(cohort),
                'protocol_sha256':'synthetic','predictions_sha256':biology.sha(predictions),
                'n':246,'positive':69,'negative':177,'source_families':1,
                'metrics':{name:{'classification':classification,'validation_selected_threshold':.5,
                                'ranking':{'average_precision':classification['average_precision']}} for name in models}})
            with mock.patch.object(biology,'validate_selection_lock') as gate_check:
                report,rows=biology.analyze(root,lock_path,final_path)
                gate_check.assert_called_once()
            self.assertEqual(len(rows),1722)
            self.assertEqual(report['subgroups']['similarity_band'][2]['n'],0)
            self.assertIsNone(report['subgroups']['original_activity_bins'][0]['models']['frozen_encoder']['average_precision'])
            self.assertTrue(report['input_hashes_unchanged_after_analysis'])
            json.dumps(report,allow_nan=False)


if __name__=='__main__': unittest.main()
