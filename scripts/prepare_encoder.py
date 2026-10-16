#!/usr/bin/env python3
"""Pinned, audited MoLFormer preparation. No fitting. All data remain in ignored paths."""
from __future__ import annotations
import argparse, collections, csv, hashlib, json, os, platform, sys, time, urllib.request
from pathlib import Path

MODEL_ID = 'ibm-research/MoLFormer-XL-both-10pct'
REVISION = '361063d0ad524ef77cf39b08469f6be770dc550f'
PROTOCOL_SHA = '283d38f275ef4c91d17dc87f1049213936543353d873198ba219b9a0b5890d1f'
COHORT_SHA = '17772cc5f360464abb2d5b913f1055ca60921d71315dfc03aa93baa1d8d89c96'
EXPECTED_FILES = {
 'README.md':'fff331c6a973fa6995662ac31354c702ac0f4805d494e7506b979a3a6cee6f44',
 'config.json':'3ef9eaac8c7ca6282fd6256ed038d151bd4ff42a4ff855367e0d7197bbc1c284',
 'configuration_molformer.py':'b88ea8d4b7b5e54f4f186cc7a230eff308928020a030f153a038fd12c05e3bed',
 'modeling_molformer.py':'6f1ef72022de2c69e95661899422a7bb39a40a2cc5a6cb6216f14e9b7d84559c',
 'tokenizer.json':'3df1f2219653c44fac9fa03b7f788b372eb2544ecc176737bb9aca8411b471a5',
 'tokenizer_config.json':'56f93dc43e4383fcfc6342e09d782ef38d48ef227215a4a92ff8762fd5b8ff3e',
}

def sha256(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as f:
  for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
 return h.hexdigest()

def save_json(path,data):
 path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
 temporary=path.with_name(path.name+f'.tmp-{os.getpid()}')
 temporary.write_text(json.dumps(data,indent=2,sort_keys=True)+'\n')
 os.replace(temporary,path)

def require_hash(path,expected):
 observed=sha256(path)
 if observed!=expected:raise ValueError(f'SHA256 mismatch for {Path(path).name}: {observed}')
 return observed

def binding(root):
 return {'protocol_sha256':require_hash(root/'docs/PROTOCOL.md',PROTOCOL_SHA),
         'cohort_sha256':require_hash(root/'data/processed/cohort.csv',COHORT_SHA)}

def verify_reviewed_files(model_dir):
 return {name:require_hash(Path(model_dir)/name,digest) for name,digest in EXPECTED_FILES.items()}

def fetch_file(url,path,expected=None):
 path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
 if path.exists():
  if expected:require_hash(path,expected)
  return
 temporary=path.with_suffix(path.suffix+'.part')
 with urllib.request.urlopen(url,timeout=180) as response,temporary.open('wb') as output:
  while True:
   chunk=response.read(1024*1024)
   if not chunk:break
   output.write(chunk)
 if expected:require_hash(temporary,expected)
 temporary.replace(path)

def prepare_assets(model_dir,include_weights=False):
 model_dir=Path(model_dir);model_dir.mkdir(parents=True,exist_ok=True)
 base=f'https://huggingface.co/{MODEL_ID}/resolve/{REVISION}/'
 for name,digest in EXPECTED_FILES.items():fetch_file(base+name,model_dir/name,digest)
 manifest={'model_id':MODEL_ID,'revision':REVISION,'files':dict(EXPECTED_FILES),
  'weights_sha256':None,'weights_downloaded':False,'format':'safetensors',
  'publisher_license':'Apache-2.0 model card and code headers',
  'pretraining':'masked SMILES modeling, 10% ZINC15 and 10% PubChem',
  'pretraining_membership':'NOT_VERIFIED','supervised_PRMT5_pretraining':'No use identified in inspected publisher documentation; not a complete corpus audit',
  'review':'reports/curation/PRETRAINED_MODEL_REVIEW.md','source_base_url':base}
 if include_weights:
  url=f'https://huggingface.co/api/models/{MODEL_ID}/tree/{REVISION}?recursive=false&expand=true'
  with urllib.request.urlopen(url,timeout=90) as response:entries=json.load(response)
  record=next(x for x in entries if x.get('path')=='model.safetensors')
  expected=record['lfs']['oid']
  fetch_file(base+'model.safetensors',model_dir/'model.safetensors',expected)
  manifest['files']['model.safetensors']=require_hash(model_dir/'model.safetensors',expected)
  manifest.update(weights_sha256=expected,weights_downloaded=True,weight_bytes=(model_dir/'model.safetensors').stat().st_size,weight_hash_basis='Pinned HF tree LFS object SHA256 and downloaded bytes agree')
 save_json(model_dir/'asset_manifest.json',manifest)
 return manifest

def audit_input_rows(rows,tokenizer):
 from rdkit import Chem,rdBase
 groups=collections.defaultdict(list);prepared=[];bad=[]
 for row in rows:
  mol=Chem.MolFromSmiles(row['canonical_smiles'])
  if mol is None:raise ValueError('Invalid SMILES: '+row['compound_id'])
  smiles=Chem.MolToSmiles(mol,canonical=True,isomericSmiles=False)
  enc=tokenizer.encode(smiles,add_special_tokens=True)
  unknown=sum(i==tokenizer.token_to_id('<unk>') for i in enc.ids)
  recovered=''.join(token for token,special in zip(enc.tokens,enc.special_tokens_mask) if not special)
  item={'compound_id':row['compound_id'],'split':row['split'],'binary_y':int(row['binary_y']),
        'canonical_smiles':row['canonical_smiles'],'nonisomeric_smiles':smiles,
        'token_count':len(enc.ids),'unknown_tokens':unknown,'reconstruction_equal':recovered==smiles}
  prepared.append(item);groups[smiles].append(item)
  if unknown or len(enc.ids)>202 or recovered!=smiles:bad.append(row['compound_id'])
 cross=[g for g in groups.values() if len({x['split'] for x in g})>1]
 stats={};collision_rows=[]
 for split in sorted({r['split'] for r in rows}):
  subset=[r for r in prepared if r['split']==split]
  gs=[g for g in groups.values() if g[0]['split']==split and len({x['split'] for x in g})==1]
  multi=[g for g in gs if len(g)>1];mixed=[g for g in multi if len({x['binary_y'] for x in g})>1]
  minimum_errors=sum(min(sum(x['binary_y']==1 for x in g),sum(x['binary_y']==0 for x in g)) for g in gs)
  stats[split]={'rows':len(subset),'distinct_inputs':len(gs),'collision_groups':len(multi),
   'collision_members':sum(map(len,multi)),'mixed_label_groups':len(mixed),'mixed_label_members':sum(map(len,mixed)),
   'empirical_minimum_errors':minimum_errors,'empirical_minimum_error_fraction':minimum_errors/len(subset),
   'token_min':min(x['token_count'] for x in subset),'token_max':max(x['token_count'] for x in subset)}
  for g in multi:collision_rows.append({'split':split,'nonisomeric_smiles':g[0]['nonisomeric_smiles'],
   'compound_ids':[x['compound_id'] for x in g],'labels':[x['binary_y'] for x in g],
   'mixed_labels':len({x['binary_y'] for x in g})>1})
 report={'status':'PASS' if not bad and not cross else 'FAIL','rows':len(rows),'rdkit_version':rdBase.rdkitVersion,
  'unknown_token_count':sum(x['unknown_tokens'] for x in prepared),'token_min':min(x['token_count'] for x in prepared),
  'token_max':max(x['token_count'] for x in prepared),'overlength_rows':sum(x['token_count']>202 for x in prepared),
  'reconstruction_failures':sum(not x['reconstruction_equal'] for x in prepared),'truncation':False,
  'unsupported_compound_ids':bad,'cross_split_collision_groups':len(cross),
  'cross_split_collision_ids':[[x['compound_id'] for x in g] for g in cross], 'by_split':stats,
  'collision_bound_interpretation':'Empirical minimum integer mistakes for a deterministic classifier using only this shared representation on these rows; not a clinical error rate, population bound, or model performance.'}
 return report,prepared,collision_rows

def audit_inputs(root,model_dir):
 from tokenizers import Tokenizer,__version__ as tokenizer_version
 verify_reviewed_files(model_dir)
 tokenizer=Tokenizer.from_file(str(model_dir/'tokenizer.json'));tokenizer.no_truncation();tokenizer.no_padding()
 with (root/'data/processed/cohort.csv').open() as f:rows=list(csv.DictReader(f))
 if len(rows)!=920 or len({x['compound_id'] for x in rows})!=920:raise ValueError('Frozen cohort must have 920 unique compound IDs')
 report,prepared,collisions=audit_input_rows(rows,tokenizer);report['tokenizers_version']=tokenizer_version
 output=root/'data/model/encoder_inputs.csv';output.parent.mkdir(parents=True,exist_ok=True)
 with output.open('w',newline='') as f:
  writer=csv.DictWriter(f,fieldnames=list(prepared[0]));writer.writeheader();writer.writerows(prepared)
 report['prepared_inputs_sha256']=sha256(output)
 save_json(root/'data/model/encoder_collisions.json',collisions)
 return report,prepared

def verify_asset_manifest(model_dir):
 model_dir=Path(model_dir);manifest=json.loads((model_dir/'asset_manifest.json').read_text())
 if manifest.get('revision')!=REVISION or manifest.get('model_id')!=MODEL_ID:raise ValueError('Unreviewed model identity/revision')
 if not manifest.get('weights_downloaded') or not manifest.get('weights_sha256'):raise ValueError('Missing verified safetensors weights')
 verify_reviewed_files(model_dir)
 if manifest.get('files',{}).get('model.safetensors')!=manifest['weights_sha256']:raise ValueError('Inconsistent weight digest')
 require_hash(model_dir/'model.safetensors',manifest['weights_sha256'])
 for name,digest in EXPECTED_FILES.items():
  if manifest['files'].get(name)!=digest:raise ValueError('Unreviewed metadata digest: '+name)
 return manifest

def verify_feature_map_buffers(model,model_dir):
 """Every persistent random-feature matrix must equal the pinned checkpoint, not initialization."""
 import torch
 from safetensors import safe_open
 buffers={name:tensor for name,tensor in model.named_buffers() if name.endswith('feature_map.weight')}
 if len(buffers)!=12:raise RuntimeError('Expected 12 persistent feature-map buffers')
 checks=[]
 with safe_open(str(Path(model_dir)/'model.safetensors'),framework='pt',device='cpu') as archive:
  keys=set(archive.keys())
  for name,actual in sorted(buffers.items()):
   checkpoint_key='molformer.'+name
   if checkpoint_key not in keys:raise RuntimeError('Checkpoint is missing a feature-map buffer')
   expected=archive.get_tensor(checkpoint_key);actual_cpu=actual.detach().cpu()
   if not torch.equal(actual_cpu,expected):raise RuntimeError('Loaded feature-map buffer differs from checkpoint')
   digest=hashlib.sha256(expected.contiguous().numpy().tobytes()).hexdigest()
   checks.append({'buffer_name':name,'checkpoint_key':checkpoint_key,'shape':list(expected.shape),'tensor_sha256':digest,'exactly_loaded':True})
 aggregate=hashlib.sha256(json.dumps(checks,sort_keys=True).encode()).hexdigest()
 return {'status':'PASS','count':12,'aggregate_sha256':aggregate,'buffers':checks}

def load_encoder_bundle(model_dir,device='cpu'):
 """Audited loader used by later fitting code. Does not fit or change data."""
 model_dir=Path(model_dir).resolve();manifest=verify_asset_manifest(model_dir)
 os.environ.setdefault('HF_HOME',str(model_dir.parent/'hf_cache'))
 os.environ['HF_HUB_OFFLINE']='1';os.environ['TRANSFORMERS_OFFLINE']='1'
 import torch,transformers
 if transformers.__version__!='5.12.1':raise RuntimeError('This reviewed revision requires the pinned Transformers 5.12.1 environment')
 model,loading_info=transformers.AutoModel.from_pretrained(str(model_dir),trust_remote_code=True,local_files_only=True,
       use_safetensors=True,deterministic_eval=True,output_loading_info=True)
 loading_info={key:sorted(value) if isinstance(value,set) else value for key,value in loading_info.items()}
 if loading_info.get('missing_keys') or loading_info.get('mismatched_keys') or loading_info.get('error_msgs'):
  raise RuntimeError('Encoder checkpoint has missing, mismatched, or erroneous weights')
 unexpected=loading_info.get('unexpected_keys',[])
 if any(not name.startswith('lm_head.') for name in unexpected):raise RuntimeError('Unexpected non-MLM checkpoint keys')
 tokenizer=transformers.AutoTokenizer.from_pretrained(str(model_dir),trust_remote_code=True,local_files_only=True)
 model=model.to(device).eval()
 if not model.config.deterministic_eval or len(model.encoder.layer)!=12:raise RuntimeError('Unexpected encoder contract')
 receipt={'revision':REVISION,'model_id':MODEL_ID,'manifest_sha256':sha256(model_dir/'asset_manifest.json'),'asset_manifest_sha256':sha256(model_dir/'asset_manifest.json'),
  'weights_sha256':manifest['weights_sha256'],'code_sha256':{k:v for k,v in EXPECTED_FILES.items() if k.endswith('.py')},
  'deterministic_eval':True,'transformers':transformers.__version__,'torch':torch.__version__,'device':str(device),
  'pooling':'model.pooler_output: attention-mask-weighted mean including BOS/EOS','encoder_blocks':12,'feature_dimension':768,
  'loaded_parameter_count':sum(p.numel() for p in model.parameters()),'loading_info':loading_info,
  'feature_map_buffers':verify_feature_map_buffers(model,model_dir)}
 return model,tokenizer,receipt

def extract_features(model,tokenizer,smiles,device,batch_size=32):
 import numpy as np,torch
 batches=[]
 with torch.no_grad():
  for start in range(0,len(smiles),batch_size):
   inputs=tokenizer(smiles[start:start+batch_size],padding=True,truncation=False,return_tensors='pt')
   if inputs['input_ids'].shape[1]>202:raise RuntimeError('Token length contract violated at inference')
   if (inputs['input_ids']==tokenizer.unk_token_id).any():raise RuntimeError('Unknown token during inference')
   output=model(**{k:v.to(device) for k,v in inputs.items()})
   pooled=output.pooler_output
   if pooled.shape!=(len(smiles[start:start+batch_size]),768) or not torch.isfinite(pooled).all():raise RuntimeError('Nonfinite or wrong-size features')
   batches.append(pooled.detach().cpu().numpy())
 return np.concatenate(batches)

def runtime_preflight(root,model_dir,prepared):
 import numpy as np,torch
 if not torch.cuda.is_available():raise RuntimeError('CUDA unavailable; GPU proof required')
 torch.manual_seed(17);torch.cuda.manual_seed_all(17)
 torch.backends.cuda.matmul.allow_tf32=False;torch.backends.cudnn.allow_tf32=False
 a=torch.arange(16,dtype=torch.float32,device='cuda').reshape(4,4)
 product=a@a.T;torch.cuda.synchronize()
 if not torch.isfinite(product).all():raise RuntimeError('CUDA tensor probe failed')
 model,tokenizer,loader=load_encoder_bundle(model_dir,'cuda')
 inputs=[r['nonisomeric_smiles'] for r in prepared]
 started=time.monotonic();features=extract_features(model,tokenizer,inputs,'cuda',32)
 # Cover shortest/longest molecules plus all three splits; every primary feature is checked finite.
 chosen=set()
 for split in ['train','validation','test']:
  ix=[i for i,r in enumerate(prepared) if r['split']==split]
  ix.sort(key=lambda i:prepared[i]['token_count']);chosen.update(ix[:8]+ix[-8:])
 indices=sorted(chosen);sample=[inputs[i] for i in indices]
 repeat=extract_features(model,tokenizer,sample,'cuda',32)
 other=extract_features(model,tokenizer,sample,'cuda',7)
 baseline=features[indices]
 repeat_ok=np.allclose(baseline,repeat,atol=1e-5,rtol=1e-5);batch_ok=np.allclose(repeat,other,atol=1e-5,rtol=1e-5)
 for name,ok in [('repeat',repeat_ok),('different_batch',batch_ok)]:
  if not ok:raise RuntimeError(name+' feature reproducibility failed')
 seed_checks=[]
 seed_reference=None
 for seed in [17,29,43]:
  torch.manual_seed(seed);torch.cuda.manual_seed_all(seed)
  fresh,fresh_tokenizer,fresh_receipt=load_encoder_bundle(model_dir,'cuda')
  seed_features=extract_features(fresh,fresh_tokenizer,inputs[:2],'cuda',2)
  if seed_reference is None:seed_reference=seed_features
  if not np.allclose(seed_reference,seed_features,atol=1e-5,rtol=1e-5):raise RuntimeError('Initialization seed changes fixed encoder predictions')
  if fresh_receipt['feature_map_buffers']['aggregate_sha256']!=loader['feature_map_buffers']['aggregate_sha256']:raise RuntimeError('Initialization seed changes feature-map buffers')
  seed_checks.append({'seed':seed,'feature_map_aggregate_sha256':fresh_receipt['feature_map_buffers']['aggregate_sha256'],'max_abs_difference_vs_seed17':float(np.max(np.abs(seed_reference-seed_features)))})
  del fresh
 output=root/'data/model/encoder_features.npz'
 temporary=output.with_name(output.stem+f'.tmp-{os.getpid()}.npz')
 np.savez_compressed(temporary,features=features,compound_ids=np.array([r['compound_id'] for r in prepared]),splits=np.array([r['split'] for r in prepared]))
 os.replace(temporary,output)
 return {'status':'PASS','loader':loader,'cuda_tensor_probe':{'status':'PASS','checksum':float(product.sum().cpu())},
  'cuda_version':torch.version.cuda,'device_capability':list(torch.cuda.get_device_capability()),'compiled_cuda_architectures':torch.cuda.get_arch_list(),
  'all_features_shape':list(features.shape),'all_features_finite':bool(np.isfinite(features).all()),
  'determinism_sample_count':len(indices),'determinism_sample_ids':[prepared[i]['compound_id'] for i in indices],
  'atol':1e-5,'rtol':1e-5,'repeat_pass':bool(repeat_ok),'different_batch_pass':bool(batch_ok),
  'repeat_max_abs_difference':float(np.max(np.abs(baseline-repeat))),
  'different_batch_max_abs_difference':float(np.max(np.abs(repeat-other))),
  'initialization_seed_checks':{'status':'PASS','sample_count':2,'results':seed_checks},
  'features_sha256':sha256(output),'feature_cache':'data/model/encoder_features.npz',
  'elapsed_seconds':round(time.monotonic()-started,3),'python':platform.python_version()}

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--project-root',type=Path,default=Path(__file__).resolve().parents[1]);ap.add_argument('--model-dir',type=Path)
 ap.add_argument('--stage',choices=['inputs','assets','runtime','all'],default='inputs');ap.add_argument('--with-weights',action='store_true');args=ap.parse_args()
 root=args.project_root.resolve();model_dir=args.model_dir or root/'data/model'/('molformer-'+REVISION)
 bound=binding(root);receipt_path=root/'reports/models/encoder_preflight.json';receipt={'status':'PENDING','revision':REVISION,'model_revision':REVISION,'bindings':bound,**bound}
 failure_stage=args.stage
 try:
  if args.stage in ['assets','all']:
   manifest=prepare_assets(model_dir,include_weights=args.with_weights or args.stage=='all')
   save_json(root/'manifests/pretrained_model.json',manifest)
  if args.stage in ['inputs','runtime','all']:
   input_gate,prepared=audit_inputs(root,model_dir);receipt['input_gate']=input_gate
   receipt['status']='INPUT_PASS_RUNTIME_PENDING' if input_gate['status']=='PASS' else 'FAIL'
   save_json(receipt_path,receipt)
   if input_gate['status']!='PASS':raise RuntimeError('Input gate failed')
  if args.stage in ['runtime','all']:
   failure_stage='runtime'
   runtime=runtime_preflight(root,model_dir,prepared);receipt['runtime_gate']=runtime
   receipt['model_asset_manifest_sha256']=runtime['loader']['manifest_sha256']
   receipt['weights_sha256']=runtime['loader']['weights_sha256'];receipt['status']='PASS'
  if args.stage!='assets':save_json(receipt_path,receipt)
  print(json.dumps(receipt if args.stage!='assets' else manifest,indent=2))
 except Exception as error:
  import traceback
  private_log=root/'receipts/encoder_preparation_failure.log';private_log.parent.mkdir(parents=True,exist_ok=True)
  with private_log.open('a') as log:log.write(traceback.format_exc()+'\n')
  public_message=str(error).replace(str(root),'<project>').replace(str(model_dir),'<model>')
  if failure_stage=='runtime':receipt['runtime_gate']={'status':'FAIL','error_type':type(error).__name__,'error':public_message}
  receipt.update(status='FAIL',failure_stage=failure_stage,error_type=type(error).__name__,error=public_message)
  save_json(receipt_path,receipt);raise

if __name__=='__main__':main()
