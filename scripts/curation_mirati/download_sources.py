#!/usr/bin/env python3
"""Download only public patent HTML. Existing snapshots are preserved unless --refresh requested."""
import argparse,json,hashlib,datetime,urllib.request
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--source-dir','--output-dir',dest='output_dir',type=Path,default=Path(__file__).resolve().parents[2] / 'data' / 'curation' / 'mirati');p.add_argument('--refresh',action='store_true');a=p.parse_args();a.output_dir.mkdir(parents=True,exist_ok=True)
manifest=[]
for patent in ['US11479551B2','US11492351B2','US20210078994A1','US20210079003A1','WO2021050915A1']:
 path=a.output_dir/(patent+'.html');url=f'https://patents.google.com/patent/{patent}/en'
 if not path.exists() or a.refresh:
  with urllib.request.urlopen(url,timeout=90) as r:path.write_bytes(r.read())
 manifest.append(dict(patent=patent,url=url,file=path.name,sha256=hashlib.sha256(path.read_bytes()).hexdigest(),bytes=path.stat().st_size,retrieved_utc=datetime.datetime.fromtimestamp(path.stat().st_mtime,datetime.timezone.utc).isoformat()))
(a.output_dir/'download_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
