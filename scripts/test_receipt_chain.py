"""Synthetic arm-receipt linkage tests; no real model inputs or labels."""
import copy
import csv
import json
from pathlib import Path
import tempfile
import unittest

from benchmark_common import sha256
from selection_lock import (BASELINE_FAMILIES, MODEL_NAMES,
                            require_merged_predictions, verified_arm_predictions)


class ReceiptChain(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.lock = {'protocol_sha256':'synthetic_protocol', 'cohort_sha256':'synthetic_cohort',
                     'code_sha256':{'synthetic_code':'code_hash'},
                     'development_artifacts':{'synthetic_development':'artifact_hash'},
                     'models':{name:{'candidate':name} for name in MODEL_NAMES}}
        path = self.root/'manifests/model_selection_lock.json'
        path.parent.mkdir(parents=True)
        path.write_text(json.dumps(self.lock))
        self.ids = [f'synthetic_{i:03d}' for i in range(246)]
        for arm, names in [('baselines', BASELINE_FAMILIES),
                           ('encoder', {'frozen_encoder','adapted_encoder'})]:
            folder = self.root/'data/predictions'/arm
            folder.mkdir(parents=True)
            rows = [{'model':name,'compound_id':cid,'split':'test','score':'0.5'}
                    for name in sorted(names) for cid in self.ids]
            self.write_rows(folder/'predictions.csv', rows)
            receipt = {key:copy.deepcopy(self.lock[key]) for key in
                       ('protocol_sha256','cohort_sha256','code_sha256','development_artifacts')}
            receipt.update(status='LOCKED_TEST_PREDICTIONS_NO_METRICS', kind=arm,
                target_columns_accessed=False, fitting_performed=False, metrics_computed=False,
                model_selection_lock_sha256=sha256(path),
                models={name:self.lock['models'][name] for name in names},
                compound_ids=self.ids, test_input_count=246, prediction_row_count=len(rows),
                predictions_sha256=sha256(folder/'predictions.csv'))
            (folder/'receipt.json').write_text(json.dumps(receipt))

    @staticmethod
    def write_rows(path, rows):
        with path.open('w', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=['model','compound_id','split','score'])
            writer.writeheader()
            writer.writerows(rows)

    def test_valid_chain_and_reordered_merge_pass(self):
        rows, receipts = verified_arm_predictions(self.root, self.lock)
        self.assertEqual(len(rows), 7*246)
        self.assertEqual(len(receipts), 2)
        merged = self.root/'merged.csv'
        self.write_rows(merged, list(reversed(rows)))
        require_merged_predictions(merged, rows)

    def test_receipt_binding_or_scope_mutation_fails(self):
        path = self.root/'data/predictions/encoder/receipt.json'
        original = json.loads(path.read_text())
        changes = [('model_selection_lock_sha256','different_lock'),
                   ('protocol_sha256','different_protocol'), ('cohort_sha256','different_cohort'),
                   ('code_sha256',{}), ('development_artifacts',{}), ('models',{}),
                   ('metrics_computed',True), ('target_columns_accessed',True),
                   ('fitting_performed',True), ('prediction_row_count',1)]
        for key, value in changes:
            with self.subTest(field=key):
                altered = copy.deepcopy(original)
                altered[key] = value
                path.write_text(json.dumps(altered))
                with self.assertRaises(ValueError):
                    verified_arm_predictions(self.root, self.lock)
        path.write_text(json.dumps(original))

    def test_prediction_bytes_tamper_fails_even_with_unchanged_receipt(self):
        path = self.root/'data/predictions/baselines/predictions.csv'
        path.write_text(path.read_text().replace('0.5','0.9',1))
        with self.assertRaisesRegex(ValueError,'bytes differ'):
            verified_arm_predictions(self.root, self.lock)

    def test_duplicate_row_fails_even_when_checksum_and_count_are_resealed(self):
        folder = self.root/'data/predictions/encoder'
        path = folder/'predictions.csv'
        with path.open(newline='') as handle:
            rows = list(csv.DictReader(handle))
        rows[-1] = dict(rows[0])
        self.write_rows(path, rows)
        receipt_path = folder/'receipt.json'
        receipt = json.loads(receipt_path.read_text())
        receipt['predictions_sha256'] = sha256(path)
        receipt_path.write_text(json.dumps(receipt))
        with self.assertRaisesRegex(ValueError,'coverage'):
            verified_arm_predictions(self.root, self.lock)

    def test_merged_score_rewrite_or_omission_fails(self):
        rows, _ = verified_arm_predictions(self.root, self.lock)
        merged = self.root/'merged.csv'
        changed = copy.deepcopy(rows)
        changed[0]['score'] = '0.9'
        for altered in (changed, rows[:-1]):
            self.write_rows(merged, altered)
            with self.assertRaisesRegex(ValueError,'differ from'):
                require_merged_predictions(merged, rows)


if __name__ == '__main__':
    unittest.main()
