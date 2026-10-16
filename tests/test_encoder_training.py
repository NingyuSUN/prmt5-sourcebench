"""Offline synthetic training tests; no PRMT5 fitting or model downloads."""
import copy
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest import mock

import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
import torch
from torch import nn

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
import train_encoder as training


class TinyEncoder(nn.Module):
    def __init__(self, blocks=4, width=4):
        super().__init__()
        self.embeddings = nn.Embedding(16, width)
        self.encoder = nn.Module()
        self.encoder.layer = nn.ModuleList([nn.Linear(width, width) for _ in range(blocks)])
        self.dropout = nn.Dropout(0.75)
        self.config = SimpleNamespace(deterministic_eval=True)

    def forward(self, input_ids, attention_mask):
        mask = attention_mask.unsqueeze(-1)
        values = (self.embeddings(input_ids) * mask).sum(1) / mask.sum(1)
        for block in self.encoder.layer:
            values = torch.tanh(block(values))
        return SimpleNamespace(pooler_output=self.dropout(values))


class TinyTokenizer:
    unk_token_id = 15

    def __call__(self, strings, **kwargs):
        if kwargs.get('truncation') is not False:
            raise AssertionError('Truncation must explicitly stay off')
        encoded = [[int(x) for x in text.split()] for text in strings]
        width = max(map(len, encoded))
        return {'input_ids':torch.tensor([row+[0]*(width-len(row)) for row in encoded]),
                'attention_mask':torch.tensor([[1]*len(row)+[0]*(width-len(row)) for row in encoded])}


def fit_head(features, labels, C, split):
    if split != 'train':
        raise AssertionError('Only training rows may fit the scaler/classifier')
    return Pipeline([('scaler', StandardScaler()),
                     ('classifier', LogisticRegression(C=C, class_weight='balanced', max_iter=5000))]).fit(features, labels)


def metric(labels, scores, ids, fraction=.1):
    return {'average_precision':float(average_precision_score(labels, scores))}


class EncoderTraining(unittest.TestCase):
    def setUp(self):
        torch.manual_seed(7)
        self.features = np.array([[1.,0,1,0],[2,1,0,0],[3,0,0,1],
                                  [4,1,1,0],[5,0,1,1],[6,1,0,1]])
        self.labels = np.array([0,0,1,1,0,1])
        self.pipeline = fit_head(self.features, self.labels, 1., 'train')

    def make_head(self):
        return training.ScaledLinearHead.from_pipeline(self.pipeline, device='cpu')

    def test_scaler_and_initial_logits_equal_frozen_classifier(self):
        head = self.make_head()
        values = torch.tensor(self.features, dtype=torch.float32)
        observed = torch.sigmoid(head(values)).detach().numpy()
        np.testing.assert_allclose(observed, self.pipeline.predict_proba(self.features)[:,1], atol=1e-6)
        np.testing.assert_allclose(head.mean.cpu().numpy(), self.features.mean(0))
        self.assertEqual(set(dict(head.named_parameters())), {'linear.weight','linear.bias'})
        self.assertFalse(head.mean.requires_grad)
        self.assertFalse(head.scale.requires_grad)

    def test_frozen_selection_scaler_uses_training_rows_only(self):
        validation = self.features + 10000.
        calls = []
        def fitting(features, labels, C, split):
            calls.append((features.copy(), labels.copy(), C, split))
            return fit_head(features, labels, C, split)
        selected = training.select_frozen_head(self.features, self.labels, validation,
            self.labels, [str(i) for i in range(len(self.labels))], fitting, metric)
        self.assertEqual([x[2] for x in calls], [.1,1.,10.])
        self.assertTrue(all(x[3]=='train' for x in calls))
        for features, labels, _, _ in calls:
            np.testing.assert_array_equal(features, self.features)
            np.testing.assert_array_equal(labels, self.labels)
        np.testing.assert_allclose(selected['pipeline'].named_steps['scaler'].mean_, self.features.mean(0))

    def test_updates_only_last_two_blocks_and_head_with_encoder_eval(self):
        model = TinyEncoder()
        head = self.make_head()
        metadata = training.configure_last_two_blocks(model, head, expected_blocks=4)
        expected = {name for name, _ in model.named_parameters()
                    if name.startswith(('encoder.layer.2.','encoder.layer.3.'))}
        self.assertEqual(set(metadata['encoder_parameter_names']), expected)
        self.assertEqual(metadata['encoder_parameter_count'], sum(p.numel() for p in model.parameters() if p.requires_grad))
        before = {name:p.detach().clone() for name,p in model.named_parameters()}
        head_before = head.linear.weight.detach().clone()
        optimizer = training.make_optimizer(model, head, learning_rate=.01)
        labels = np.array([0,1,0,1,0,1,0,1])
        loss = training.train_epoch(model, head, TinyTokenizer(),
            ['1 2','2 3','3 4','4 5','5 6','6 7','7 8','8 9'], labels, optimizer,
            training.balanced_class_weights(labels), device='cpu', micro_batch_size=2)
        self.assertTrue(np.isfinite(loss))
        self.assertFalse(model.training)
        self.assertFalse(model.dropout.training)
        for name, parameter in model.named_parameters():
            if name not in expected:
                torch.testing.assert_close(parameter, before[name], rtol=0, atol=0)
        for block in ('encoder.layer.2.','encoder.layer.3.'):
            self.assertTrue(any(not torch.equal(parameter, before[name]) for name,parameter in model.named_parameters() if name.startswith(block)))
        self.assertFalse(torch.equal(head.linear.weight, head_before))

    def test_shared_parameter_cannot_unfreeze_an_earlier_block(self):
        model = TinyEncoder()
        model.encoder.layer[3] = model.encoder.layer[0]
        with self.assertRaisesRegex(ValueError, 'shared'):
            training.configure_last_two_blocks(model, self.make_head(), expected_blocks=4)

    def test_balanced_weights_use_training_denominator(self):
        np.testing.assert_allclose(training.balanced_class_weights([0,0,0,1]), [4/6,2])
        with self.assertRaises(ValueError):
            training.balanced_class_weights([1,1])

    def test_checkpoint_is_a_copy_and_restores_best_epoch(self):
        model, head = TinyEncoder(), self.make_head()
        training.configure_last_two_blocks(model, head, expected_blocks=4)
        optimizer = training.make_optimizer(model, head, learning_rate=.01)
        training.train_epoch(model, head, TinyTokenizer(), ['1 2','2 3','3 4','4 5'],
            np.array([0,1,0,1]), optimizer, np.array([1.,1.]), 'cpu', micro_batch_size=2)
        snapshot = training.capture_training_state(model, head, optimizer, 1, {'average_precision':.8}, 'cpu')
        reference = copy.deepcopy(snapshot)
        training.train_epoch(model, head, TinyTokenizer(), ['1 2','2 3','3 4','4 5'],
            np.array([0,1,0,1]), optimizer, np.array([1.,1.]), 'cpu', micro_batch_size=2)
        for name, value in snapshot['encoder_state'].items():
            torch.testing.assert_close(value, reference['encoder_state'][name], rtol=0, atol=0)
        training.restore_training_weights(model, head, snapshot)
        for name, parameter in model.named_parameters():
            if parameter.requires_grad:
                torch.testing.assert_close(parameter.cpu(), reference['encoder_state'][name], rtol=0, atol=0)
        torch.testing.assert_close(head.linear.weight.cpu(), reference['head_state']['linear.weight'], rtol=0, atol=0)
        self.assertEqual(snapshot['epoch'], 1)
        self.assertTrue(snapshot['optimizer_state']['state'])
        self.assertIn('torch_cpu', snapshot['rng_state'])

    def test_early_stop_keeps_epoch_copy_not_last_epoch(self):
        model, head = TinyEncoder(), self.make_head()
        initial_bias = float(head.linear.bias.detach())
        def synthetic_epoch(*args, **kwargs):
            with torch.no_grad():
                head.linear.bias.add_(1.)
            return .5
        initial_scores = training.predict_probabilities(model, head, TinyTokenizer(), ['3 4','4 5'], 'cpu')
        sequence = iter([.8,.7,.6])
        def validation_metric(*args, **kwargs):
            return {'average_precision':next(sequence)}
        with tempfile.TemporaryDirectory() as directory:
            with mock.patch.object(training, 'train_epoch', side_effect=synthetic_epoch):
                result = training.run_trial(model, TinyTokenizer(), head,
                    ['1 2','2 3'], np.array([0,1]), ['3 4','4 5'], np.array([0,1]),
                    ['v1','v2'], 1e-5, 17, 'cpu', Path(directory), validation_metric,
                    max_epochs=3, patience=2, expected_blocks=4,
                    initial_validation_probabilities=initial_scores,
                    checkpoint_context={
                        'model_receipt':{'torch':torch.__version__},
                        'runtime_versions':{'torch':torch.__version__},
                        'nested':[{'count':np.int64(2),'value':np.float32(.25),
                                   'flags':(np.bool_(True), None)}]})
            self.assertEqual(result['best_epoch'], 1)
            self.assertEqual(len(result['training_trace']), 3)
            self.assertAlmostEqual(float(head.linear.bias.detach()), initial_bias+1., places=6)
            checkpoint = torch.load(Path(directory)/'best.pt', weights_only=True)
            self.assertEqual(checkpoint['epoch'], 1)
            self.assertIs(type(checkpoint['context']['model_receipt']['torch']), str)
            self.assertIs(type(checkpoint['context']['runtime_versions']['torch']), str)
            nested = checkpoint['context']['nested'][0]
            self.assertIs(type(nested['count']), int)
            self.assertIs(type(nested['value']), float)
            self.assertEqual(nested['flags'], [True, None])
            rng = checkpoint['rng_state']['numpy']
            self.assertTrue(all(type(value) is int for value in rng['keys']))
            restored = np.random.RandomState()
            restored.set_state((rng['algorithm'],np.asarray(rng['keys'],dtype=np.uint32),
                                rng['position'],rng['has_gaussian'],rng['cached_gaussian']))
            self.assertTrue(np.isfinite(restored.random_sample()))

    def test_context_rejects_unsupported_objects_and_nonfinite_scalars(self):
        for value in ({'bad':object()}, {'bad':np.float64('nan')}):
            with self.assertRaises((TypeError,ValueError)):
                training.checkpoint_primitives(value)

    def test_initialization_mismatch_rejects_before_optimizer_exists(self):
        model, head = TinyEncoder(), self.make_head()
        scores = training.predict_probabilities(model, head, TinyTokenizer(), ['3 4','4 5'], 'cpu')
        with tempfile.TemporaryDirectory() as directory:
            with mock.patch.object(training, 'make_optimizer') as optimizer:
                with self.assertRaisesRegex(ValueError, 'initialization'):
                    training.run_trial(model, TinyTokenizer(), head,
                        ['1 2','2 3'], np.array([0,1]), ['3 4','4 5'], np.array([0,1]),
                        ['v1','v2'], 1e-5, 29, 'cpu', Path(directory), metric,
                        max_epochs=1, expected_blocks=4,
                        initial_validation_probabilities=np.clip(scores+.1,0,1))
                optimizer.assert_not_called()

    def test_nonfinite_json_cannot_replace_valid_artifact(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'receipt.json'
            training.save_json(path, {'status':'VALID'})
            before = path.read_bytes()
            with self.assertRaises(ValueError):
                training.save_json(path, {'bad':float('nan')})
            self.assertEqual(path.read_bytes(), before)
            self.assertEqual(len(list(Path(directory).iterdir())), 1)

    def test_source_hash_change_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'source.txt'
            path.write_text('frozen')
            paths = {'public_source':path}
            snapshot = training.snapshot_hashes(paths)
            training.verify_snapshot(paths, snapshot)
            path.write_text('changed')
            with self.assertRaisesRegex(ValueError, 'source'):
                training.verify_snapshot(paths, snapshot)

    def test_learning_rate_selection_averages_all_seeds_and_ties_lower(self):
        trials = [{'learning_rate':lr,'seed':seed,'best_validation_average_precision':score}
                  for lr,scores in [(1e-5,[.99,.1,.1]),(3e-5,[.5,.5,.5])]
                  for seed,score in zip([17,29,43],scores)]
        selected = training.select_learning_rate(trials)
        self.assertEqual(selected['selected_learning_rate'], 3e-5)
        for trial in trials:
            trial['best_validation_average_precision'] = .5
        self.assertEqual(training.select_learning_rate(trials)['selected_learning_rate'], 1e-5)
        with self.assertRaisesRegex(ValueError, 'seeds'):
            training.select_learning_rate(trials[:-1])


if __name__ == '__main__':
    unittest.main()
