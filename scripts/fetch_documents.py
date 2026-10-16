import argparse, concurrent.futures, datetime, hashlib, importlib.util, json, time
from pathlib import Path
p = argparse.ArgumentParser()
p.add_argument('--root', type=Path, required=True)
p.add_argument('--client', type=Path, required=True)
a = p.parse_args()
s = importlib.util.spec_from_file_location('client', a.client)
m = importlib.util.module_from_spec(s)
s.loader.exec_module(m)
ids = sorted({r['document_chembl_id'] for f in (a.root / 'data/raw').glob('*_activity_0*.json') for r in json.loads(f.read_text())['activities']})

def fetch(pair):
    i, batch = pair
    path = a.root / f'data/raw/documents_{i:03d}.json'
    params = {'document_chembl_id__in': ','.join(batch), 'limit': 100}
    for attempt in range(3):
        out = m.execute(dict(base_url='https://www.ebi.ac.uk/chembl/api/data', path='document.json', params=params, record_path='documents', max_items=1, timeout_sec=50, save_raw=True, raw_output_path=str(path)))
        if out['ok']:
            break
        time.sleep(2)
    if not out['ok']:
        raise RuntimeError(out)
    d = json.loads(path.read_text())
    print(i, len(d['documents']), flush=True)
    return dict(path=str(path.relative_to(a.root)), url='https://www.ebi.ac.uk/chembl/api/data/document.json', params=params, rows=len(d['documents']), sha256=hashlib.sha256(path.read_bytes()).hexdigest(), accessed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as pool:
    results = list(pool.map(fetch, enumerate([ids[i:i + 40] for i in range(0, len(ids), 40)])))
(a.root / 'manifests/document_downloads.json').write_text(json.dumps(results, indent=2) + '\n')
print('DONE', len(ids), sum((x['rows'] for x in results)))
