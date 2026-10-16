import csv,json,collections,hashlib,datetime
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

O=_WORK; rows=list(csv.DictReader(open(O/'patent_crosswalk.csv'))); out=[]
for a in ['CHEMBL5737316','CHEMBL5739093']:
 by=collections.defaultdict(list)
 for x in rows:
  if x['assay_id']==a:by[x['canonical_smiles']].append(x)
 for sm,rr in by.items():
  vv={float(x['database_value_nM']) for x in rr};labels={int(v<=100) for v in vv}; controls=a=='CHEMBL5739093' and any(x['source_examples'] in ['C1','C2','C3'] for x in rr)
  out.append({'assay_id':a,'canonical_smiles':sm,'scaffold_group':rr[0]['scaffold_group'],'binary_y':next(iter(labels)) if len(labels)==1 else '', 'decision':'quarantine_threshold_conflict' if len(labels)>1 else 'reference_control_exclude_temporal_novel' if controls else 'candidate_include','activity_ids':';'.join(x['activity_id'] for x in rr),'source_examples':';'.join(sorted(set(x['source_examples'] for x in rr))),'source_values_original':';'.join(sorted(set(x['source_values'] for x in rr))),'database_values_nM_preserved_only':';'.join(str(v) for v in sorted(vv)),'source_date_verified':rr[0]['first_verified_source_date'],'source_date_type':'activity_table_publication_verified_not_exhaustive_first_pair_date','numeric_precision':'interval_only_do_not_regress' if a=='CHEMBL5737316' else 'source_numeric','catalytic_endpoint':'PRMT5/MEP50 histoneH4 methylation_inhibition','source_protocol_held_separate':'True'})
with open(O/'binary_parent_candidates.csv','w') as f:
 w=csv.DictWriter(f,fieldnames=out[0]);w.writeheader();w.writerows(out)
stats={}
for a in ['CHEMBL5737316','CHEMBL5739093']:
 rr=[x for x in out if x['assay_id']==a]; inc=[x for x in rr if x['decision']=='candidate_include'];stats[a]={'all_parents':len(rr),'decision_counts':dict(collections.Counter(x['decision'] for x in rr)),'included_parents':len(inc),'positive_at_most_100nM':sum(x['binary_y']==1 for x in inc),'negative_above_100nM':sum(x['binary_y']==0 for x in inc),'scaffolds':len({x['scaffold_group'] for x in inc}),'labels_are_source_conditioned':True}
(O/'binary_candidate_summary.json').write_text(json.dumps(stats,indent=2)+'\n');print(json.dumps(stats,indent=2))
# Every temporary artifact has a digest; original download times in source_manifest.json.
m=[]
for p in sorted(O.iterdir()):
 if p.is_file() and p.name!='artifact_hashes.json':m.append({'path':str(p.relative_to(O)),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
(O/'artifact_hashes.json').write_text(json.dumps(m,indent=2)+'\n')

assert len(rows)==329 and len(out)==323
assert [(stats[a]['included_parents'],stats[a]['positive_at_most_100nM'],stats[a]['negative_above_100nM']) for a in ['CHEMBL5737316','CHEMBL5739093']]==[(246,69,177),(70,39,31)]
