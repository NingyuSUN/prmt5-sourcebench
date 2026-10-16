"""Small semantic safeguards for exact-value inclusion and overlap quarantine."""
import unittest
from validate_sources import exact_positive_nm_ic50, strict_reasons
from cohort_audit import novel_against, parent_conflicts


class AuditSemantics(unittest.TestCase):
    def setUp(self):
        self.row = dict(standard_type='IC50', standard_relation='=', standard_units='nM',
                        standard_value='1', target_tax_id='9606')
        self.assay = dict(confidence_score=7, bao_label='protein complex format', assay_type='B')

    def test_relations_units_nonpositive_nonfinite_are_not_exact_nm(self):
        self.assertTrue(exact_positive_nm_ic50(self.row))
        for field, value in [('standard_relation', '<'), ('standard_relation', '>'),
            ('standard_units', 'uM'), ('standard_type', 'Ki'), ('standard_value', '0'),
            ('standard_value', '-1'), ('standard_value', 'NaN'), ('standard_value', 'inf')]:
            with self.subTest(field=field, value=value):
                self.assertFalse(exact_positive_nm_ic50({**self.row, field: value}))

    def test_strict_target_confidence_and_flags_are_not_silently_relaxed(self):
        self.assertEqual(strict_reasons(self.row, self.assay, ''), [])
        self.assertIn('not_direct_target', strict_reasons(self.row, {**self.assay, 'confidence_score': 6}, ''))
        self.assertIn('potential_duplicate', strict_reasons({**self.row, 'potential_duplicate': 1}, self.assay, ''))
        self.assertIn('structure_invalid', strict_reasons(self.row, self.assay, 'ValueError'))

    def test_seen_scaffolds_and_old_parent_years_are_quarantined(self):
        train = [dict(canonical_smiles='a', scaffold_group='x')]
        query = [dict(canonical_smiles='b', scaffold_group='x'),
                 dict(canonical_smiles='c', scaffold_group='y'),
                 dict(canonical_smiles='d', scaffold_group='z')]
        self.assertEqual(novel_against(train, query), query[1:])
        self.assertEqual(novel_against(train, query, {'b': 2023, 'c': 2020, 'd': 2023}, 2022), query[2:])

    def test_numeric_representation_is_not_a_value_conflict(self):
        base = dict(canonical_smiles='a', scaffold_group='x', molecule_chembl_id='M',
            assay_chembl_id='A', document_chembl_id='D')
        rows = [{**base, 'activity_id': '1', 'standard_value': '0.50'},
                {**base, 'activity_id': '2', 'standard_value': '0.5'}]
        self.assertEqual(parent_conflicts(rows)['distinct_value_parents'], 0)
        rows.append({**base, 'activity_id': '3', 'standard_value': '1510'})
        result = parent_conflicts(rows)
        self.assertEqual(result['distinct_value_parents'], 1)
        self.assertEqual(result['maximum_fold_spread'], 3020)


if __name__ == '__main__':
    unittest.main()
