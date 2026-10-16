#!/usr/bin/env python3
"""Report indirect same-molecule, same-family alias candidates separately from direct mappings."""
import argparse,csv,json,re,collections
from pathlib import Path
from curate_mirati import write_csv
p=argparse.ArgumentParser();p.add_argument('--source-dir',type=Path,default=Path(__file__).resolve().parents[2] / 'data' / 'curation' / 'mirati');a=p.parse_args();d=a.source_dir
rows=list(csv.DictReader((d/'mirati_activity_mapping.csv').open()));sources=list(csv.DictReader((d/'mirati_source_table_rows.csv').open()));aliases=list(csv.DictReader((d/'mirati_compound_aliases.csv').open()))
by_molecule=collections.defaultdict(set)
for r in aliases:
 by_molecule[r['molecule_chembl_id']].update(re.findall(r'US(?:11479551|11492351),\s*Example\s+(\d+-\d+)',r['compound_name']))
result=[]
for r in rows:
 if r['mapping_status']!='unmapped_value_or_condition':continue
 family_examples=by_molecule[r['molecule_chembl_id']]
 candidates=[s for s in sources if s['patent']==r['patent'] and s['source_example'] in family_examples and s['mta_condition']==r['mta_condition'] and s['parse_status']=='numeric' and s['relation']==r['relation'] and float(s['value_nM'])==float(r['value_nM'])]
 for s in candidates:
  result.append(dict(activity_id=r['activity_id'],molecule_chembl_id=r['molecule_chembl_id'],patent=r['patent'],direct_aliases=r['record_aliases'],family_molecule_aliases='|'.join(sorted(family_examples)),candidate_count=len(candidates),table=s['table'],source_example=s['source_example'],example_occurrence=s['example_occurrence'],mta_condition=s['mta_condition'],relation=s['relation'],value_nM=s['value_nM'],source_locator=s['html_locator'],source_sha256=s['source_sha256'],evidence_status='INFERRED_family_same_molecule_alias_transfer',basis='Four US original tables verified identical; same ChEMBL molecule has example aliases in sibling grant; exact endpoint condition relation value match',limitation='Not a direct document-specific alias; source structure images not independently decoded'))
write_csv(d/'mirati_family_alias_candidates.csv',result)
print(dict(unmatched_activities=sum(r['mapping_status']=='unmapped_value_or_condition' for r in rows),activities_with_candidates=len(set(r['activity_id'] for r in result)),unique_indirect=sum(r['candidate_count']==1 for r in result),ambiguous_indirect_activity_ids=sorted(set(r['activity_id'] for r in result if r['candidate_count']>1))))
