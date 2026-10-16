#!/usr/bin/env python3
"""Parse saved Google Patents HTML; standard-library only. Never use values alone to identify compounds."""
import argparse,csv,hashlib,html,json,re,collections,datetime
from pathlib import Path

def text(s):return ' '.join(html.unescape(re.sub('<[^>]+>','',s)).split())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def write_csv(p,rows):
 if rows:
  with p.open('w',newline='') as f:
   w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def parse(p):
 raw=p.read_text();pat=p.stem;out=[];table_info=[];source_sha=sha(p)
 dates=re.findall(r'<meta name="DC.date" content="([^"]+)"([^>]*)>',raw)
 date=next((d for d,a in dates if 'dateSubmitted' not in a),'')
 family=re.search(r'<h2>ID=(\d+)</h2>',raw)
 for wrapper in re.finditer(r'<tables\b([^>]*)>(.*?)</tables>',raw,re.S):
  content=wrapper.group(2);m=re.search(r'TABLE\s+(9[ab])',content)
  if not m:continue
  table=m.group(1);wid=re.search(r'id="([^"]+)"',wrapper.group(1))[1]
  blocks=list(re.finditer(r'<table\b[^>]*>(.*?)</table>',content,re.S))
  table_info.append(dict(patent=pat,table=table,html_id=wid,header=text(blocks[0][1]),publication_date=date))
  for ti,b in enumerate(blocks,1):
   for ri,row in enumerate(re.finditer(r'<tr\b[^>]*>(.*?)</tr>',b[1],re.S),1):
    cells=list(re.finditer(r'<t[dh]\b[^>]*>(.*?)</t[dh]>',row[1],re.S)); vals=[text(c[1]) for c in cells]
    for offset in (0,3):
     if len(vals)<=offset or not re.fullmatch(r'\d+-\d+',vals[offset]):continue
     assert len(vals)>=offset+3,(pat,table,vals)
     for delta,mta in [(1,'MTA_2uM'),(2,'MTA_absent')]:
      v=vals[offset+delta];vnorm=v.replace(',','').replace(' ','')
      vm=re.fullmatch(r'([<>≤≥]=?)?(\d+(?:\.\d+)?)',vnorm)
      status='numeric' if vm else 'not_determined' if vnorm in ['N.D.','ND','N.D',''] else 'not_reported_dash' if vnorm=='—' else 'unparsed'
      line=raw.count('\n',0,wrapper.start(2)+b.start(1)+row.start())+1
      out.append(dict(patent=pat,family_id=family[1] if family else '',publication_date=date,table=table,source_example=vals[offset],endpoint='IC50',mta_condition=mta,relation=(vm[1] or '=') if vm else '',value_nM=str(float(vm[2])) if vm else '',raw_value=v,parse_status=status,protocol='FlashPlate_H4' if table=='9a' else 'HotSpot_H2A',enzyme_nM=3 if table=='9a' else 1,substrate='H4_1-15_biotin_40nM' if table=='9a' else 'full_length_H2A_5uM',SAM_uM=1,preincubation_min=20,reaction_min=20 if table=='9a' else 60,source_url=f'https://patents.google.com/patent/{pat}/en',source_file=p.name,source_sha256=source_sha,html_locator=f'#{wid} table:nth-of-type({ti}) tr:nth-of-type({ri}) td:nth-of-type({offset+delta+1})',html_line=line))
 counts=collections.Counter()
 for r in out:
  key=(r['table'],r['source_example'],r['mta_condition']);counts[key]+=1;r['example_occurrence']=counts[key]
 return out,table_info,dict(patent=pat,source_file=p.name,sha256=sha(p),bytes=p.stat().st_size,retrieved_utc=datetime.datetime.fromtimestamp(p.stat().st_mtime,datetime.timezone.utc).isoformat(),publication_date=date,family_id=family[1] if family else '',structured_table_count=len(table_info))

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--source-dir',type=Path,default=Path(__file__).resolve().parents[2] / 'data' / 'curation' / 'mirati');args=ap.parse_args();d=args.source_dir
 rows=[];infos=[];manifest=[]
 for p in sorted(d.glob('*.html')):
  r,i,m=parse(p);rows+=r;infos+=i;manifest.append(m)
 write_csv(d/'mirati_source_table_rows.csv',rows)
 write_csv(d/'mirati_table_headers.csv',infos)
 comparisons=[]
 pairs=[('US11479551B2','US20210078994A1'),('US11492351B2','US20210079003A1'),('US11479551B2','US11492351B2'),('US20210078994A1','US20210079003A1')]
 for later,earlier in pairs:
  a={(x['table'],x['source_example'],x['mta_condition'],x['example_occurrence']):x for x in rows if x['patent']==later};b={(x['table'],x['source_example'],x['mta_condition'],x['example_occurrence']):x for x in rows if x['patent']==earlier}
  assert len(a)==sum(x['patent']==later for x in rows);assert len(b)==sum(x['patent']==earlier for x in rows)
  for k in sorted(a.keys()|b.keys()):
   x=a.get(k);y=b.get(k)
   equal=bool(x and y and all(x[f]==y[f] for f in ['relation','value_nM','parse_status','raw_value']))
   comparisons.append(dict(later_patent=later,earlier_patent=earlier,table=k[0],source_example=k[1],mta_condition=k[2],example_occurrence=k[3],status='identical' if equal else 'only_later' if not y else 'only_earlier' if not x else 'different',later_raw=x['raw_value'] if x else '',earlier_raw=y['raw_value'] if y else '',later_locator=x['html_locator'] if x else '',earlier_locator=y['html_locator'] if y else ''))
 write_csv(d/'mirati_disclosure_comparison.csv',comparisons)
 summary={'source_manifest':manifest,'parse_counts':[dict(patent=k[0],table=k[1],status=k[2],cells=n) for k,n in collections.Counter((r['patent'],r['table'],r['parse_status']) for r in rows).items()],'comparison_counts':[dict(later=k[0],earlier=k[1],status=k[2],cells=n) for k,n in collections.Counter((r['later_patent'],r['earlier_patent'],r['status']) for r in comparisons).items()],'limits':['US table9b titles say FlashPlate; preceding H2A protocol supports HotSpot. WO title says HotSpot.','WO tables are images and not transcribed; row comparisons use US A1 structured tables.','2021-03-18 is earliest verified table disclosure among inspected documents; priority date is not public disclosure; worldwide earlier absence is not proven.']}
 (d/'mirati_parse_summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
