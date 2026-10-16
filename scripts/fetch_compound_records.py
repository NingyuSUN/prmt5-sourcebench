#!/usr/bin/env python3
"""Retrieve public source-example crosswalks for explicitly selected assays."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
from pathlib import Path
import time

SELECTED_ASSAYS = {
    'CHEMBL5731839', 'CHEMBL5732555', 'CHEMBL5732693',
    'CHEMBL5733925', 'CHEMBL5734324', 'CHEMBL5736675',
    'CHEMBL5738640', 'CHEMBL5738641', 'CHEMBL5738688', 'CHEMBL5738689',
    'CHEMBL5737316', 'CHEMBL5739093', 'CHEMBL5544428', 'CHEMBL6088129',
}

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--client', type=Path, required=True)
    args = parser.parse_args()
    spec = importlib.util.spec_from_file_location('chembl_client', args.client)
    client = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(client)
    docs = sorted({row['document_chembl_id']
                   for file in (args.root / 'data/raw').glob('*_assay_0*.json')
                   for row in json.loads(file.read_text())['assays']
                   if row['assay_chembl_id'] in SELECTED_ASSAYS})
    def fetch_document(document):
        output = []
        offset = 0
        while True:
            path = args.root / f'data/raw/compound_records_{document}_{offset:06d}.json'
            params = {'document_chembl_id': document, 'limit': 1000, 'offset': offset}
            for attempt in range(3):
                response = client.execute({
                    'base_url': 'https://www.ebi.ac.uk/chembl/api/data',
                    'path': 'compound_record.json', 'params': params,
                    'record_path': 'compound_records', 'max_items': 1,
                    'timeout_sec': 60, 'save_raw': True, 'raw_output_path': str(path),
                })
                if response['ok']:
                    break
                time.sleep(attempt + 1)
            if not response['ok']:
                raise RuntimeError(response)
            data = json.loads(path.read_text())
            output.append({
                'path': str(path.relative_to(args.root)),
                'url': 'https://www.ebi.ac.uk/chembl/api/data/compound_record.json',
                'params': params, 'rows': len(data['compound_records']),
                'total_count': data['page_meta']['total_count'],
                'sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
                'bytes': path.stat().st_size,
                'retrieved_utc': datetime.now(timezone.utc).isoformat(),
            })
            print(document, offset, len(data['compound_records']), flush=True)
            if data['page_meta']['next'] is None:
                return output
            offset += 1000
    with ThreadPoolExecutor(max_workers=3) as pool:
        results = [row for batch in pool.map(fetch_document, docs) for row in batch]
    (args.root / 'manifests/compound_record_downloads.json').write_text(
        json.dumps(results, indent=2) + '\n')
    print('Completed', len(docs), 'documents', sum(row['rows'] for row in results), 'source records')

if __name__ == '__main__':
    main()
