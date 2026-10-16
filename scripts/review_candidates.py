#!/usr/bin/env python3
"""List candidate assays requiring primary-source review, without admitting them."""
import argparse
import csv
import itertools
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    args = parser.parse_args()
    root = args.root
    with (root / 'reports/activity_ledger.csv').open() as stream:
        rows = list(csv.DictReader(stream))
    documents = {
        record['document_chembl_id']: record
        for path in (root / 'data/raw').glob('documents_*.json')
        for record in json.loads(path.read_text())['documents']
    }
    output = []
    grouped = itertools.groupby(sorted(rows, key=lambda r: r['assay_chembl_id']),
                                key=lambda r: r['assay_chembl_id'])
    for assay_id, group in grouped:
        group = list(group)
        candidates = [row for row in group
                      if set(row['exclusion_reasons'].split(';')) - {''}
                      <= {'not_direct_target', 'potential_duplicate'}]
        if len(candidates) < 10:
            continue
        unflagged = [row for row in candidates
                     if 'potential_duplicate' not in row['exclusion_reasons']]
        document = documents[group[0]['document_chembl_id']]
        output.append({
            'assay_id': assay_id,
            'rows': len(candidates),
            'nonflagged': len(unflagged),
            'parents': len({row['canonical_smiles'] for row in candidates}),
            'scaffolds': len({row['scaffold_group'] for row in candidates}),
            'confidence': group[0]['confidence'],
            'year': document['year'],
            'document_id': document['document_chembl_id'],
            'patent_id': document['patent_id'],
            'doi': document['doi'],
            'title': document['title'],
            'description': group[0]['assay_description'],
        })
    output.sort(key=lambda row: row['parents'], reverse=True)
    (root / 'reports/review_candidates.json').write_text(json.dumps(output, indent=2) + '\n')
    print(f'Saved {len(output)} assay-review candidates; none newly admitted for training.')


if __name__ == '__main__':
    main()
