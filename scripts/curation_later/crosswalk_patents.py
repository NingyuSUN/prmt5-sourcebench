"""Original patent table/record-name crosswalk. Does not infer values from positions."""
import re,json,csv,html,collections,hashlib
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

R=_ROOT;O=_WORK
def clean(s):return re.sub(r'\s+',' ',html.unescape(re.sub('<[^>]+>',' ',s))).strip()
def table(doc,tid):
 s=(O/(doc+'.html')).read_text();t=re.search(r'<tables\b[^>]*id="'+tid+r'"[^>]*>(.*?)</tables>',s,re.S).group(1)
 return [[clean(c) for c in re.findall(r'<td\b[^>]*>(.*?)</td>',r,re.S)] for r in re.findall(r'<tr\b[^>]*>(.*?)</tr>',t,re.S)]
skt=table('US11111237B2','TABLE-US-00002');ska=table('US20210101892A1','TABLE-US-00002')
assert skt==ska
skt={x[1]:x[2] for x in skt if len(x)==3 and x[1].startswith('Example ')}
pbt={x[1]:x[2] for x in table('US11591326B2','TABLE-US-00003') if len(x)==3 and re.fullmatch(r'C?\d+[ab]?',x[1])}
allraw={}
for p in (R/'data/raw').glob('*activity_[0-9]*.json'):
 for x in json.load(open(p))['activities']:
  if x['assay_chembl_id'] in ['CHEMBL5737316','CHEMBL5739093'] and x['standard_type']=='IC50' and x['standard_value'] is not None:allraw[x['activity_id']]=x
ledger={int(x['activity_id']):x for x in csv.DictReader(open(R/'reports/activity_ledger.csv'))}
result=[]
for aid,doc,pat,tbl,pattern in [('CHEMBL5737316','CHEMBL5728719','US11111237',skt,r'US11111237, (Example \d+)'),('CHEMBL5739093','CHEMBL5729615','US11591326',pbt,r'US11591326, Compound (C?\d+[ab]?)')]:
 recs={x['record_id']:x for p in (R/'data/raw').glob(f'compound_records_{doc}_*.json') for x in json.load(open(p))['compound_records']}
 rr=[x for x in allraw.values() if x['assay_chembl_id']==aid]
 for x in rr:
  record=recs[x['record_id']]; examples=re.findall(pattern,record.get('compound_name') or '');vals={v:tbl[v] for v in examples if v in tbl}
  v=float(x['standard_value'])
  surrogate={'+':5500.,'++':550.,'+++':50.}
  expected=[surrogate[t] if pat=='US11111237' else float(t) for t in vals.values() if t!='—']
  item={'assay_id':aid,'activity_id':x['activity_id'],'record_id':x['record_id'],'molecule_id':x['molecule_chembl_id'],'compound_name':record['compound_name'],'bindingdb_key':record['compound_key'],'source_examples':';'.join(examples),'source_values':json.dumps(vals),'database_value_nM':v,'has_matching_source_value_or_surrogate':v in expected,'exact_activity_to_example_unique':len(examples)==1,'canonical_smiles':ledger[x['activity_id']]['canonical_smiles'],'scaffold_group':ledger[x['activity_id']]['scaffold_group'],'interpretation':'INTERVAL_ONLY' if pat=='US11111237' else 'NUMERIC_SOURCE'}
  result.append(item)
  if not expected or v not in expected: print('UNMATCHED',item)
 print(aid,len(rr),'matched',sum(x['has_matching_source_value_or_surrogate'] for x in result if x['assay_id']==aid),'records',len(recs),'source_table_values',len(tbl))
assert len(result)==329 and all(x['has_matching_source_value_or_surrogate'] for x in result)
manual=json.loads((O/'wo2019_table2_manual.json').read_text())['values_nM']
assert manual=={k:None if v=='—' else float(v) for k,v in pbt.items()}
for item in result:
 item['earlier_source_verified']='True'
 item['first_verified_source_date']='2021-04-08' if item['assay_id']=='CHEMBL5737316' else '2019-09-12'
 item['date_evidence']='US20210101892A1 Table2 exact text match plus Table1 structure names match' if item['assay_id']=='CHEMBL5737316' else 'WO2019173804A1 originalTable2 all75cells manually checked againstUS11591326B2; exactstructure onlysource-familycontext noteachimage'
with open(O/'patent_crosswalk.csv','w') as f:
 w=csv.DictWriter(f,fieldnames=list(result[0]));w.writeheader();w.writerows(result)
by=collections.defaultdict(list)
for x in result:
 if x['assay_id']=='CHEMBL5737316':by[x['canonical_smiles']].append(x)
conflicts=[]
for sm,rows in by.items():
 vals=sorted({x['database_value_nM'] for x in rows})
 if len(vals)>1:conflicts.append({'canonical_smiles':sm,'values_surrogates_nM':vals,'straddles_100nM':min(vals)<=100<max(vals),'source_rows':rows})
(O/'sk_bin_conflicts.json').write_text(json.dumps(conflicts,indent=2)+'\n')
print('SK conflicting parents',len(conflicts),'binary conflicting',sum(x['straddles_100nM'] for x in conflicts))

assert len(conflicts)==5 and sum(x['straddles_100nM'] for x in conflicts)==4
