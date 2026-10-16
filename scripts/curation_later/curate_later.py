"""Later PRMT5 feasibility only; preserves source measurements, no model fitting."""
import csv,json,re,html,collections,hashlib,datetime
from pathlib import Path
import argparse
_parser = argparse.ArgumentParser(description=__doc__)
_parser.add_argument('--root', type=Path, required=True, help='Canonical project root')
_parser.add_argument('--work-dir', type=Path, required=True, help='Local source and output directory')
_args = _parser.parse_args()
_ROOT = _args.root.resolve()
_WORK = _args.work_dir.resolve()
def _verify_inputs():
    manifest = json.loads((_WORK / 'reproduction_inputs.json').read_text())
    for entry in manifest:
        base = _ROOT if entry['base'] == 'root' else _WORK
        path = base / entry['path']
        if hashlib.sha256(path.read_bytes()).hexdigest() != entry['sha256']:
            raise RuntimeError('Input hash mismatch: ' + str(path))
_verify_inputs()

from rdkit import Chem,DataStructs
from rdkit.Chem import rdFingerprintGenerator
ROOT=_ROOT; OUT=_WORK
ASES=['CHEMBL5737316','CHEMBL5739093','CHEMBL5544428','CHEMBL6088129','CHEMBL5509028']
allrows=list(csv.DictReader(open(ROOT/'reports/activity_ledger.csv')))
# Match original candidate inclusion only, retaining confidence-6 and duplicate annotations for source audit.
def cand(x):return set(x['exclusion_reasons'].split(';'))-{''} <= {'not_direct_target','potential_duplicate'}
rows=[x for x in allrows if cand(x)]
cs={x['assay_id']:x for x in json.load(open(ROOT/'reports/review_candidates.json'))}
# Older source scope is explicit, it is not a frozen training dataset.
ep=[x for x in rows if x['assay_description'].find('75 nM')>=0 and 'SAM' in x['assay_description']]
mi=[x for x in rows if x['assay_chembl_id'] in ['CHEMBL5738640','CHEMBL5738641','CHEMBL5738688','CHEMBL5738689','CHEMBL5048053','CHEMBL5048052']]
fpgen=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048,includeChirality=True)
def structures(rs):return {x['canonical_smiles']:x for x in rs}
def fprints(rs):return [fpgen.GetFingerprint(Chem.MolFromSmiles(s)) for s in structures(rs)]
bfps={'epizyme_75nM_text_candidates':fprints(ep),'mirati_selected_assays':fprints(mi)}
out={'verified_at_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'scope':'Data feasibility only. Cohorts are candidates, not admitted model labels. Baseline overlap uses candidate rows and text-selected 75 nM scope, not frozen training sets.','cohorts':[]}
for a in ASES:
 rr=[x for x in rows if x['assay_chembl_id']==a];dd=structures(rr); ps=set(dd);ss={x['scaffold_group'] for x in rr}
 by=collections.defaultdict(set)
 for x in rr:by[x['canonical_smiles']].add(float(x['standard_value']))
 item={**cs[a],'unique_values':len({x['standard_value'] for x in rr}),'value_counts':dict(collections.Counter(x['standard_value'] for x in rr)),'parents_multiple_values':sum(len(v)>1 for v in by.values()),'parent_classes_100nM':dict(collections.Counter('ambiguous' if min(v)<=100<max(v) else 'at_most_100nM' if max(v)<=100 else 'above_100nM' for v in by.values())),'overlap':{}}
 for name,base in [('epizyme_75nM_text_candidates',ep),('mirati_selected_assays',mi)]:
  bps={x['canonical_smiles'] for x in base};bss={x['scaffold_group'] for x in base}
  fresh=[s for s in ps if s not in bps];sims=[max(DataStructs.BulkTanimotoSimilarity(fpgen.GetFingerprint(Chem.MolFromSmiles(s)),bfps[name]),default=0) for s in fresh]
  item['overlap'][name]={'baseline_rows':len(base),'baseline_parents':len(bps),'exact_parent_overlap':len(ps&bps),'parents_no_exact_overlap':len(fresh),'parents_no_scaffold_overlap':sum(dd[s]['scaffold_group'] not in bss for s in fresh),'fresh_parent_max_similarity_counts':{str(t):sum(v>=t for v in sims) for t in [0.5,0.7,0.8,0.9]},'similarity':'Morgan radius2 2048bit with chirality, Tanimoto; no claim of scaffold novelty from a cutoff alone'}
 out['cohorts'].append(item)
(OUT/'later_capacity.json').write_text(json.dumps(out,indent=2)+'\n')
for x in out['cohorts']:print(x['assay_id'],x['rows'],x['parents'],x['scaffolds'],x['parent_classes_100nM'],x['parents_multiple_values'],{n:(v['exact_parent_overlap'],v['parents_no_scaffold_overlap']) for n,v in x['overlap'].items()})

assert [(x['rows'],x['parents']) for x in out['cohorts']] == [(256,250),(73,73),(27,27),(30,30),(20,20)]
