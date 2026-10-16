#!/usr/bin/env python3
"""Validate parsed coverage, raw-metadata joins, explicit ambiguous retention, and earliest-disclosure comparisons."""
import argparse,csv,json,collections
from pathlib import Path
ap=argparse.ArgumentParser();ap.add_argument('--source-dir',type=Path,default=Path(__file__).resolve().parents[2] / 'data' / 'curation' / 'mirati');D=ap.parse_args().source_dir
read=lambda n:list(csv.DictReader((D/n).open()))
s=read('mirati_source_table_rows.csv');m=read('mirati_activity_mapping.csv');c=read('mirati_disclosure_comparison.csv');f=read('mirati_family_alias_candidates.csv')
assert len(s)==2848 and not any(x['parse_status']=='unparsed' for x in s)
assert collections.Counter(x['patent'] for x in s)==dict.fromkeys(['US11479551B2','US11492351B2','US20210078994A1','US20210079003A1'],712)
assert len(c)==2848 and all(x['status']=='identical' for x in c)
assert len(m)==3378 and len({x['activity_id'] for x in m})==3378
assert collections.Counter(x['mapping_status'] for x in m)=={'unique':1108,'excluded_non_IC50':2252,'ambiguous_source_occurrences':6,'unmapped_value_or_condition':12}
assert len([x for x in m if x['endpoint']=='IC50' and float(x['value_nM'])==0])==8
assert len(set(x['activity_id'] for x in f))==12
assert len([x for x in f if x['candidate_count']=='1'])==10
index={(x['patent'],x['html_locator']):x for x in s}
for x in m:
 if x['mapping_status']=='unique':
  y=index[x['patent'],x['source_locator']]
  assert x['source_example'] in x['record_aliases'].split('|')
  assert float(x['value_nM'])==float(y['value_nM']) and x['relation']==y['relation'] and x['mta_condition']==y['mta_condition']
  assert y['parse_status']=='numeric'
known={(x['table'],x['mta_condition']):x['value_nM'] for x in s if x['patent']=='US11479551B2' and x['source_example']=='4-100'}
assert known=={('9a','MTA_2uM'):'0.5',('9a','MTA_absent'):'80.0',('9b','MTA_2uM'):'1510.0',('9b','MTA_absent'):'10000.0'}
result={'status':'PASS','scope':'Engineering extraction and metadata joins only; no scientific or external validation','source_cells':len(s),'comparisons':len(c),'activity_rows':len(m),'direct_unique':1108,'direct_ambiguous':6,'direct_unmapped':12,'excluded_non_IC50':2252,'nonpositive_IC50_records':8,'indirect_unique_candidates':10,'indirect_ambiguous_activities':2}
(D/'mirati_validation.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
