from pathlib import Path
import json,csv,re,collections,hashlib,datetime
import argparse
parser = argparse.ArgumentParser(description="Recover source intervals and conservative Epizyme candidates; never treat bin representatives as exact IC50.")
parser.add_argument('--project-root', type=Path, default=Path(__file__).resolve().parents[2])
parser.add_argument('--patent-dir', type=Path, required=True)
parser.add_argument('--output-dir', type=Path, required=True)
args = parser.parse_args()
P = args.project_root.resolve()
PATENTS = args.patent_dir.resolve()
OUT = args.output_dir.resolve()
OUT.mkdir(parents=True, exist_ok=True)
ASSAYS={'CHEMBL5732693':('US9777008B2','WO2014100695A1','US20140228343A1'),'CHEMBL5732555':('US9745291B2','WO2014100730A1','US20140221345A1'),'CHEMBL5733925':('US10307413B2','WO2014100719A2','US20140329794A1'),'CHEMBL5734324':('US10391089B2','WO2014100719A2','US20140329794A1'),'CHEMBL5736675':('US10980794B2','WO2014100719A2','US20140329794A1'),'CHEMBL5731839':('US10150758B2','WO2014100716A1','US20140213582A1')}
DOCS={p.stem:json.load(open(p)) for p in PATENTS.glob('*.json') if re.fullmatch(r'(?:US|WO)\d+[AB]\d',p.stem)}
def table(pid):
 ts=DOCS[pid]['tables'];inds=[i for i,t in enumerate(ts) if 'Biochemical IC' in t['text']];rows=[]
 if inds:
  for t in ts[inds[0]:]:
   for row in t['rows']:
    cells=[c['text'] for c in row]
    while len(cells)>4 and not cells[0]:cells.pop(0)
    if len(cells)==4 and re.fullmatch(r'\d+[A-Za-z]?',cells[0]):rows.append(cells)
 (OUT/(pid+'_table2.json')).write_text(json.dumps(rows,indent=2));return {r[0]:r[1:] for r in rows}
TABLES={i:table(i) for v in ASSAYS.values() for i in (v[0],v[2])}
records={}
for p in sorted((P/'data/raw').glob('compound_records_*.json')):
 for r in json.load(open(p)).get('compound_records',[]):records[r['record_id']]=r
raw={}
for p in sorted((P/'data/raw').glob('*activity_*.json')):
 if 'probe' in p.name:continue
 for r in json.load(open(p))['activities']:
  if r['assay_chembl_id'] in ASSAYS:raw[str(r['activity_id'])]=r
ledger={r['activity_id']:r for r in csv.DictReader(open(P/'reports/activity_ledger.csv'))}
rows=[]
for aid,a in sorted(raw.items(), key=lambda item: int(item[0])):
 if a['standard_type']!='IC50':continue
 l=ledger[aid];pat,wo,early=ASSAYS[a['assay_chembl_id']];rec=records.get(a['record_id'],{});name=rec.get('compound_name','');prefix=pat[:-2]
 nos=list(dict.fromkeys(re.findall(re.escape(prefix)+r',\s*(?:Compound|Cmpd No)\s+(\d+[A-Za-z]?)\b',name,re.I)))
 candidates={no:TABLES[pat].get(no) for no in nos};cats={v[0] for v in candidates.values() if v};complete=bool(nos) and all(candidates.values());category=next(iter(cats)) if complete and len(cats)==1 else None
 top=40000 if pat=='US9777008B2' else 50000
 bins={'A':(None,100,'<',100),'B':(101,1000,'[]',550),'C':(1001,10000,'[]',5500),'D':(10001,top,'[]',25000 if top==40000 else 30000),'E':(top,None,'>',top),'F':(1000,None,'>',1000),'G':(5000,None,'>',5000),'*':(10000,None,'>',10000),'**':(20000,None,'>',20000),'—':(None,None,'missing',None)}
 low,up,closure,rep=bins.get(category,(None,None,'unmapped',None))
 earlymatch=complete and all(TABLES[early].get(no)==v for no,v in candidates.items())
 reason=[]
 if not nos:reason.append('no_patent_specific_compound_name')
 if any(v is None for v in candidates.values()):reason.append('compound_number_not_in_B2_table')
 if len(cats)>1:reason.append('multiple_aliases_have_conflicting_original_bins')
 if l['structure_error'] or not l['canonical_smiles']:reason.append('structure_invalid')
 if category in [None,'','—']:reason.append('no_unique_source_label')
 label=1 if category == 'A' else (0 if category in ['B','C','D','E','F','G','*','**'] else None)
 # Same-number A1 matches are locators only, not identity or date evidence.
 # The directly indexed B2 public date is a conservative upper bound, never first disclosure.
 bound=DOCS[pat]['meta']['DC.date'] # conservative indexed source, not same-number cross-publication identity
 r={'activity_id':aid,'record_id':a['record_id'],'assay_chembl_id':a['assay_chembl_id'],'document_chembl_id':a['document_chembl_id'],'patent':pat,'family':DOCS[pat]['family'],'compound_name':name,'source_compound_numbers':';'.join(nos),'candidate_original_rows':json.dumps(candidates),'original_bin':category,'lower_nM':low,'upper_nM':up,'bound_semantics':closure,'binary_active_at_le100nM':label,'chembl_relation':a['standard_relation'],'chembl_value_nM':a['standard_value'],'chembl_expected_representative_nM':rep,'chembl_representative_matches_source_bin':rep is not None and float(a['standard_value'])==rep,'value_is_exact':False,'original_bin_mapped':category is not None,'mapping_reasons':';'.join(reason),'early_A1':early,'early_A1_date':DOCS[early]['meta']['DC.date'],'early_A1_same_number_all3categories_match':earlymatch,'WO_candidate':wo,'WO_publication_date':'2014-06-26','WO_row_verified':False,'direct_source_publication_date':DOCS[pat]['meta']['DC.date'],'record_disclosed_by_date_bound':bound,'earliest_exact_date_verified':False,'canonical_smiles':l['canonical_smiles'],'scaffold_group':l['scaffold_group'],'structure_error':l['structure_error'],'data_validity_comment':l['data_validity_comment'],'potential_duplicate':l['potential_duplicate'],'parent_molecule_chembl_id':l['parent_molecule_chembl_id'],'molecule_chembl_id':l['molecule_chembl_id'],'structure_identity_in_original_figure':'not_independently_checked'}
 rows.append(r)
parents=collections.defaultdict(set)
for r in rows:
 if r['binary_active_at_le100nM'] is not None:parents[r['canonical_smiles']].add(r['binary_active_at_le100nM'])
for r in rows:
 r['parent_has_conflicting_source_binary_labels']=len(parents[r['canonical_smiles']])>1
 r['eligible_source_label_only']=bool(not r['mapping_reasons'] and not r['parent_has_conflicting_source_binary_labels'])
with (OUT/'epizyme_activity_curation.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
summary={}
def count(rr):return {'rows':len(rr),'parents':len({r['canonical_smiles'] for r in rr}),'scaffolds':len({r['scaffold_group'] for r in rr}),'row_labels':dict(collections.Counter(r['binary_active_at_le100nM'] for r in rr))}
for assay,(pat,wo,early) in ASSAYS.items():
 rr=[r for r in rows if r['assay_chembl_id']==assay];ok=[r for r in rr if r['eligible_source_label_only']];summary[assay]={'patent':pat,'family':DOCS[pat]['family'],'all_IC50':count(rr),'table_rows':len(TABLES[pat]),'source_bin_mapped_rows':sum(r['original_bin_mapped'] for r in rr),'categories':dict(collections.Counter(r['original_bin'] for r in rr)),'source_label_eligible':count(ok),'early_A1':early,'early_A1_date':DOCS[early]['meta']['DC.date'],'same_number_all3categories_in_early_A1':sum(r['early_A1_same_number_all3categories_match'] for r in rr),'exceptions':[{'activity_id':r['activity_id'],'compound_name':r['compound_name'],'mapping_reasons':r['mapping_reasons'],'source_bin':r['original_bin'],'chembl_value':r['chembl_value_nM'],'parent_binary_conflict':r['parent_has_conflicting_source_binary_labels']} for r in rr if r['mapping_reasons'] or not r['chembl_representative_matches_source_bin'] or r['parent_has_conflicting_source_binary_labels']]}
summary['overall']={'all_IC50':count(rows),'source_label_eligible':count([r for r in rows if r['eligible_source_label_only']]),'exact_continuous_values':0,'earliest_dates_verified':0,'source_family_count':4,'training_readiness':'NOT_ESTABLISHED_BY_THIS_CURATION_ALONE'}

whitelist = {'50001259', '50029211', '50001260'}
candidates = [r for r in rows if r['family'] in whitelist and r['eligible_source_label_only'] and r['record_disclosed_by_date_bound'] < '2019-01-01']
parent_records = {}
for row in candidates:
 parent_records.setdefault(row['canonical_smiles'], row)
unique_parents = list(parent_records.values())
summary['conservative_candidates'] = {
 'filter': "family in {50001259,50029211,50001260} AND eligible_source_label_only=true AND record_disclosed_by_date_bound<2019-01-01; deduplicate by canonical_smiles for parent counts",
 'counts': count(candidates),
 'parent_labels': dict(collections.Counter(row['binary_active_at_le100nM'] for row in unique_parents)),
 'threshold_nM': 100,
 'positive_rule': 'source A strictly below100nM; no source reports the unobserved exact100nM boundary',
 'quarantined_family': '49958717',
 'quarantine_reason': 'source numbering differs after276; same database aliases can connect inconsistent original categories, including seven parents with binary conflicts across22rows; original structure-figure identity not independently confirmed',
 'candidate_status': 'source-label and conservative-date capacity only; assay comparability and final split validation handled separately'
}
for filename, subset in [('epizyme_conservative_rows.csv', candidates), ('epizyme_conservative_parents.csv', unique_parents)]:
 with (OUT/filename).open('w') as handle:
  writer=csv.DictWriter(handle,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(subset)
(OUT/'epizyme_curation_summary.json').write_text(json.dumps(summary,indent=2))
print(json.dumps(summary['conservative_candidates'],indent=2))
