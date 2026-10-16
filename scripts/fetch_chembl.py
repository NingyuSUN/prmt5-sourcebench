"""Snapshot public ChEMBL PRMT5 records through the ChEMBL skill client."""
import argparse, concurrent.futures, datetime, hashlib, importlib.util, json, time
from pathlib import Path

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--root', type=Path, required=True)
    p.add_argument('--client', type=Path, required=True)
    a = p.parse_args()
    spec = importlib.util.spec_from_file_location('chembl_rest_client', a.client)
    client = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(client)
    raw = a.root / 'data/raw'
    raw.mkdir(parents=True, exist_ok=True)

    def fetch(kind, target, offset=0, limit=1000):
        path = raw / f'{target}_{kind}_{offset:06d}.json'
        params = {'target_chembl_id': target, 'limit': limit, 'offset': offset, 'order_by': {'activity': 'activity_id', 'assay': 'assay_chembl_id'}[kind]}
        for attempt in range(3):
            out = client.execute({'base_url': 'https://www.ebi.ac.uk/chembl/api/data', 'path': kind + '.json', 'params': params, 'record_path': {'activity': 'activities', 'assay': 'assays'}[kind], 'max_items': 1, 'timeout_sec': 60, 'save_raw': True, 'raw_output_path': str(path)})
            if out['ok']:
                break
            time.sleep(2 * (attempt + 1))
        if not out['ok']:
            raise RuntimeError(out)
        d = json.loads(path.read_text())
        n = len(d[{'activity': 'activities', 'assay': 'assays'}[kind]])
        rec = {'source': 'ChEMBL', 'path': str(path.relative_to(a.root)), 'url': 'https://www.ebi.ac.uk/chembl/api/data/' + kind + '.json', 'params': params, 'accessed_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(), 'sha256': hashlib.sha256(path.read_bytes()).hexdigest(), 'bytes': path.stat().st_size, 'rows': n, 'total_count': d['page_meta']['total_count']}
        print(json.dumps({'file': path.name, 'rows': n, 'total': rec['total_count']}), flush=True)
        return rec
    jobs = [(k, t) for t in ['CHEMBL1795116', 'CHEMBL3137261'] for k in ['activity', 'assay']]
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        first = list(pool.map(lambda j: fetch(*j), jobs))
        rest = [(k, t, o) for (k, t), r in zip(jobs, first) for o in range(1000, r['total_count'], 1000)]
        records = first + list(pool.map(lambda j: fetch(*j), rest))
    (a.root / 'manifests/chembl_downloads.json').write_text(json.dumps(records, indent=2) + '\n')
    print('DONE', len(records), sum((x['rows'] for x in records)))
if __name__ == '__main__':
    main()
