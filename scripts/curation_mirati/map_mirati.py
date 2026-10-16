#!/usr/bin/env python3
"""Match activities to source tables by document-specific compound aliases + endpoint, MTA condition and exact relation/value.
No value-only compound identification. Ambiguous source occurrences remain unresolved.
"""
import argparse,csv,json,re,collections,hashlib
from pathlib import Path
from curate_mirati import write_csv,sha

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--source-dir',type=Path,default=Path(__file__).resolve().parents[2] / 'data' / 'curation' / 'mirati');ap.add_argument('--raw-dir',type=Path,required=True);a=ap.parse_args();d=a.source_dir
 sources=list(csv.DictReader((d/'mirati_source_table_rows.csv').open()))
 assay={'CHEMBL5738640':('US11479551B2','MTA_2uM'),'CHEMBL5738641':('US11479551B2','MTA_absent'),'CHEMBL5738688':('US11492351B2','MTA_2uM'),'CHEMBL5738689':('US11492351B2','MTA_absent')}
 records={};input_manifest=[]
 for p in sorted(a.raw_dir.glob('compound_records_*.json')):
  if not any(doc in p.name for doc in ['CHEMBL5729406','CHEMBL5729433']):continue
  data=json.loads(p.read_text())
  for r in data.get('compound_records',[]):records[r['record_id']]=(r,p)
  input_manifest.append(dict(file=p.name,sha256=sha(p)))
 activities={}
 for p in sorted(a.raw_dir.glob('*activity_0*.json')):
  for r in json.loads(p.read_text()).get('activities',[]):
   if r['assay_chembl_id'] in assay:activities[r['activity_id']]=(r,p)
  input_manifest.append(dict(file=p.name,sha256=sha(p)))
 result=[];candidates=[];used=collections.Counter();alias_rows=[]
 for aid,(r,p) in activities.items():
  patent,mta=assay[r['assay_chembl_id']];record,rp=records.get(r['record_id'],({},None))
  name=record.get('compound_name','');aliases=sorted(set(re.findall(re.escape(re.sub('[AB]\\d$','',patent))+r',\s*Example\s+(\d+-\d+)',name)))
  identity_ok=bool(record and record['document_chembl_id']==r['document_chembl_id'] and record['molecule_chembl_id']==r['molecule_chembl_id'])
  possible=[s for s in sources if s['patent']==patent and s['mta_condition']==mta and s['source_example'] in aliases and s['parse_status']=='numeric'] if identity_ok and r['standard_type']=='IC50' else []
  matches=[s for s in possible if s['relation']==r['standard_relation'] and float(s['value_nM'])==float(r['standard_value']) and r['standard_units']=='nM']
  status='excluded_non_IC50' if r['standard_type']!='IC50' else 'missing_record_metadata' if not record else 'record_identity_mismatch' if not identity_ok else 'no_document_specific_example_alias' if not aliases else 'unmapped_value_or_condition' if not matches else 'unique' if len(matches)==1 else 'ambiguous_source_occurrences'
  s=matches[0] if len(matches)==1 else {}
  row=dict(activity_id=aid,record_id=r['record_id'],molecule_chembl_id=r['molecule_chembl_id'],parent_molecule_chembl_id=r['parent_molecule_chembl_id'],document_chembl_id=r['document_chembl_id'],assay_chembl_id=r['assay_chembl_id'],patent=patent,endpoint=r['standard_type'],mta_condition=mta,relation=r['standard_relation'],value_nM=r['standard_value'] if r['standard_units']=='nM' else '',record_aliases='|'.join(aliases),mapping_status=status,match_count=len(matches),table=s.get('table',''),source_example=s.get('source_example',''),example_occurrence=s.get('example_occurrence',''),source_value_nM=s.get('value_nM',''),protocol=s.get('protocol',''),source_locator=s.get('html_locator',''),source_sha256=s.get('source_sha256',''),earliest_verified_table_publication='2021-03-18' if s else '',family_id=s.get('family_id',''),activity_raw_file=p.name,compound_record_file=rp.name if rp else '',mapping_basis='document_specific_example_alias+record_molecule_identity+endpoint+MTA+relation+value' if matches else '',structure_validation='ChEMBL record_id/molecule crosswalk; patent structure images not independently decoded')
  result.append(row)
  for s in matches:
   candidates.append(dict(activity_id=aid,patent=patent,table=s['table'],source_example=s['source_example'],example_occurrence=s['example_occurrence'],mta_condition=mta,value_nM=s['value_nM'],source_locator=s['html_locator'],source_sha256=s['source_sha256']))
   used[(s['patent'],s['html_locator'])]+=1
  if record:alias_rows.append(dict(record_id=record['record_id'],document_chembl_id=record['document_chembl_id'],molecule_chembl_id=record['molecule_chembl_id'],compound_name=name,patent=patent,document_specific_examples='|'.join(aliases),file=rp.name))
 write_csv(d/'mirati_activity_mapping.csv',result);write_csv(d/'mirati_mapping_candidates.csv',candidates)
 dedup={r['record_id']:r for r in alias_rows};write_csv(d/'mirati_compound_aliases.csv',list(dedup.values()))
 uncovered=[s for s in sources if s['patent'] in ['US11479551B2','US11492351B2'] and s['parse_status']=='numeric' and not used[(s['patent'],s['html_locator'])]]
 write_csv(d/'mirati_uncovered_numeric_source_cells.csv',uncovered)
 summary=dict(activity_count=len(result),mapping_status_counts=dict(collections.Counter(x['mapping_status'] for x in result)),mapped_by_protocol_condition=[dict(protocol=k[0],mta_condition=k[1],activities=n) for k,n in collections.Counter((x['protocol'],x['mta_condition']) for x in result if x['mapping_status']=='unique').items()],uncovered_numeric_source_cells=len(uncovered),multiply_used_source_cells=sum(n>1 for n in used.values()),input_manifest=input_manifest,limitations=['Structure identity follows ChEMBL compound_record and molecule crosswalk; patent structure images not independently decoded.','Exact unique table assignment is not pharmacological reproducibility or target specificity validation.','2021-03-18 is earliest verified row disclosure in inspected US A1s, not proven worldwide earliest disclosure.'])
 (d/'mirati_mapping_summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps({k:v for k,v in summary.items() if k!='input_manifest'},indent=2))
if __name__=='__main__':main()
