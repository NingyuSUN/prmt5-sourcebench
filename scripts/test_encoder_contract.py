#!/usr/bin/env python3
"""Meaningful input-contract tests; do not download weights or execute remote code."""
import tempfile,unittest
from pathlib import Path
from prepare_encoder import audit_input_rows,require_hash,save_json
from unittest.mock import patch

class FakeEncoding:
 def __init__(self,s,n=None,unknown=False):
  self.tokens=['<bos>',s,'<eos>'];self.special_tokens_mask=[1,0,1]
  self.ids=[0]+([2361] if unknown else [5])*(n-2 if n else 1)+[1]
class FakeTokenizer:
 def __init__(self,n=None,unknown=False):self.n=n;self.unknown=unknown
 def token_to_id(self,s):return 2361
 def encode(self,s,add_special_tokens=True):return FakeEncoding(s,self.n,self.unknown)
def row(cid,smiles,y,split='train'):
 return {'compound_id':cid,'canonical_smiles':smiles,'binary_y':str(y),'split':split}

class EncoderInputContract(unittest.TestCase):
 def test_stereo_collision_mixed_labels_reports_unavoidable_errors(self):
  rows=[row('a','N[C@@H](C)C(=O)O',0),row('b','N[C@H](C)C(=O)O',1),row('c','CCO',1)]
  report,prepared,collisions=audit_input_rows(rows,FakeTokenizer())
  self.assertEqual(report['status'],'PASS');self.assertEqual(report['by_split']['train']['mixed_label_groups'],1)
  self.assertEqual(report['by_split']['train']['empirical_minimum_errors'],1)
  self.assertAlmostEqual(report['by_split']['train']['empirical_minimum_error_fraction'],1/3)
 def test_cross_split_stereoisomer_collision_blocks(self):
  rows=[row('a','N[C@@H](C)C(=O)O',0),row('b','N[C@H](C)C(=O)O',1,'test')]
  report,_,_=audit_input_rows(rows,FakeTokenizer())
  self.assertEqual(report['status'],'FAIL');self.assertEqual(report['cross_split_collision_groups'],1)
 def test_202_token_boundary_and_unknown_block(self):
  rows=[row('a','CCO',1)]
  self.assertEqual(audit_input_rows(rows,FakeTokenizer(202))[0]['status'],'PASS')
  self.assertEqual(audit_input_rows(rows,FakeTokenizer(203))[0]['status'],'FAIL')
  self.assertEqual(audit_input_rows(rows,FakeTokenizer(3,True))[0]['status'],'FAIL')
 def test_corrupt_input_hash_blocks(self):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d)/'data';p.write_bytes(b'changed')
   with self.assertRaises(ValueError):require_hash(p,'0'*64)
 def test_failed_atomic_replace_keeps_previous_receipt(self):
  with tempfile.TemporaryDirectory() as d:
   p=Path(d)/'receipt.json';p.write_text('{"status":"old"}')
   with patch('prepare_encoder.os.replace',side_effect=OSError('simulated write failure')):
    with self.assertRaises(OSError):save_json(p,{'status':'new'})
   self.assertEqual(p.read_text(),'{"status":"old"}')
 def test_invalid_structure_blocks(self):
  with self.assertRaises(ValueError):audit_input_rows([row('a','definitely not smiles',1)],FakeTokenizer())

if __name__=='__main__':unittest.main()
