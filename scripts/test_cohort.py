import unittest
import hashlib
import json
import os
from pathlib import Path
import tempfile
from unittest import mock
import build_cohort as cohort
from build_cohort import interval_label, later_label, old_label, collapse

class SourceLabels(unittest.TestCase):
    def test_boundary_and_straddling_intervals(self):
        self.assertEqual(interval_label(None,100),1)
        self.assertEqual(interval_label(100,1000,False),0)
        self.assertIsNone(interval_label(100,1000,True))
        self.assertIsNone(interval_label(50,150))
        self.assertIsNone(interval_label(None,None))
    def test_nonpositive_or_reversed_invalid(self):
        for args in [(0,0),(-1,100),(1000,100)]:
            with self.assertRaises(ValueError): interval_label(*args)
        with self.assertRaises(ValueError):later_label('0')
    def test_midpoint_never_drives_label(self):
        row={'lower_nM':'1001','upper_nM':'10000','bound_semantics':'[]','chembl_value_nM':'550'}
        self.assertEqual(old_label(row),0)
        self.assertEqual(later_label('+++'),1)
        self.assertEqual(later_label('++'),0)
        self.assertEqual(later_label('100'),1)
    def test_discordant_parent_quarantined_not_averaged(self):
        rows=[{'canonical_smiles':'C','binary_y':0},{'canonical_smiles':'C','binary_y':1}]
        included,conflicts=collapse(rows)
        self.assertEqual(included,[])
        self.assertEqual(len(conflicts),2)

class CohortAdmission(unittest.TestCase):
    def test_unknown_later_assay_fails(self):
        self.assertEqual(cohort.later_source('CHEMBL5739093'), ('validation', '67846363'))
        self.assertEqual(cohort.later_source('CHEMBL5737316'), ('test', '75338422'))
        with self.assertRaisesRegex(ValueError, 'Unknown later assay'):
            cohort.later_source('UNREVIEWED')

    def test_family_overlap_expands_parent_family_sets(self):
        reference = [{'canonical_smiles':'C', 'scaffold_group':'s1', 'source_family':'a;b'}]
        query = [{'canonical_smiles':'CC', 'scaffold_group':'s2', 'source_family':'b'}]
        self.assertEqual(cohort.overlap_counts(reference, query),
                         {'canonical_smiles':0, 'scaffold_group':0, 'source_family':1})

    def test_cutoff_checked_before_untrusted_later_labels(self):
        row = {'family':'50001259', 'direct_source_publication_date':'2019-09-01',
               'original_bin_mapped':'True', 'original_bin':'A',
               'binary_active_at_le100nM':'not-a-label'}
        for future_label in ('not-a-label', '0', '1'):
            row['binary_active_at_le100nM'] = future_label
            with self.assertRaisesRegex(ValueError, 'Post-cutoff training source'):
                cohort.admit_training_label(row)
        row['family'] = 'UNRESOLVED'
        self.assertEqual(cohort.admit_training_label(row),
                         (None, 'unresolved_family_identity_mapping'))

    def test_quality_flags_preserved_and_unreviewed_warning_rejected(self):
        ledger = {'1':{'data_validity_comment':'', 'potential_duplicate':'1'},
                  '2':{'data_validity_comment':'', 'potential_duplicate':'0'}}
        flags = cohort.quality_evidence('1;2', ledger)
        self.assertEqual(flags[0], {'activity_id':'1', 'data_validity_comment':'',
                                   'potential_duplicate':'1'})
        ledger['2']['data_validity_comment'] = 'Needs review'
        with self.assertRaisesRegex(ValueError, 'Unreviewed activity validity flag'):
            cohort.quality_evidence('1;2', ledger)
        with self.assertRaisesRegex(ValueError, 'Missing activity metadata'):
            cohort.quality_evidence('3', ledger)


class AtomicPublication(unittest.TestCase):
    def test_checksums_match_published_artifacts(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            output = root / 'data' / 'processed'
            report_path = root / 'reports' / 'gate.json'
            report = {'status':'PASS'}
            cohort.publish_outputs(root, output, [('a.csv',[{'x':'one'}])], report_path, report)
            saved = json.loads(report_path.read_text())
            self.assertEqual(saved['status'], 'PASS')
            self.assertEqual(saved['outputs']['data/processed/a.csv'],
                             hashlib.sha256((output/'a.csv').read_bytes()).hexdigest())

    def test_failed_replacement_cannot_leave_pass_for_mixed_artifacts(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            output = root / 'data' / 'processed'
            output.mkdir(parents=True)
            (output/'a.csv').write_text('old-a')
            (output/'b.csv').write_text('old-b')
            report_path = root / 'reports' / 'gate.json'
            report_path.parent.mkdir()
            report_path.write_text('{"status":"PASS"}')
            real_replace = os.replace
            def fail_second_artifact(source, target):
                if Path(target).name == 'b.csv':
                    raise OSError('simulated replacement failure')
                return real_replace(source, target)
            with mock.patch('build_cohort.os.replace', side_effect=fail_second_artifact):
                with self.assertRaisesRegex(OSError, 'simulated'):
                    cohort.publish_outputs(root, output,
                        [('a.csv',[{'x':'new-a'}]),('b.csv',[{'x':'new-b'}])],
                        report_path, {'status':'PASS'})
            self.assertEqual(json.loads(report_path.read_text())['status'], 'BUILD_IN_PROGRESS')
            self.assertEqual((output/'b.csv').read_text(), 'old-b')
            self.assertEqual(sorted(p.name for p in output.iterdir()), ['a.csv','b.csv'])


if __name__=='__main__':unittest.main()
