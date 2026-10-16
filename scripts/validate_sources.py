#!/usr/bin/env python3
"""Independent stdlib source/pagination/ledger audit; no chemistry recomputation."""
import argparse
import collections
import csv
import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path


def file_record(path, root=None):
    path = Path(path)
    with path.open('rb') as stream:
        digest = hashlib.file_digest(stream, 'sha256').hexdigest()
    return {'path': str(path.relative_to(root)) if root else path.name,
            'bytes': path.stat().st_size, 'sha256': digest}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def safe_path(root, relative):
    path = (root / relative).resolve()
    require(not Path(relative).is_absolute() and path.is_relative_to(root),
            f'Unsafe manifest path: {relative}')
    return path


def exact_positive_nm_ic50(row):
    try:
        value = float(row.get('standard_value'))
    except (ValueError, TypeError):
        return False
    return (row.get('standard_type') == 'IC50'
            and row.get('standard_relation') == '='
            and row.get('standard_units') == 'nM'
            and math.isfinite(value) and value > 0)


def strict_reasons(activity, assay, structure_error):
    reasons = []
    if structure_error:
        reasons.append('structure_invalid')
    for field, expected, reason in [
        ('standard_type', 'IC50', 'not_IC50'),
        ('standard_relation', '=', 'not_exact'),
        ('standard_units', 'nM', 'not_nM')]:
        if activity.get(field) != expected:
            reasons.append(reason)
    try:
        value = float(activity.get('standard_value'))
        positive = math.isfinite(value) and value > 0
    except (ValueError, TypeError):
        positive = False
    if not positive:
        reasons.append('nonpositive_or_missing')
    if activity.get('data_validity_comment'):
        reasons.append('validity_flag')
    if activity.get('potential_duplicate'):
        reasons.append('potential_duplicate')
    if assay.get('confidence_score') not in [7, 9]:
        reasons.append('not_direct_target')
    if (assay.get('bao_label') not in ['single protein format', 'protein complex format']
            or assay.get('assay_type') != 'B'):
        reasons.append('not_biochemical_format')
    if str(activity.get('target_tax_id')) != '9606':
        reasons.append('target_not_human')
    return reasons


def counts(rows):
    def unique(field):
        return len({r[field] for r in rows if r.get(field) not in [None, '']})
    return {'rows': len(rows), 'molecules': unique('molecule_chembl_id'),
            'parents': unique('canonical_smiles'), 'scaffolds': unique('scaffold_group'),
            'assays': unique('assay_chembl_id'), 'documents': unique('document_chembl_id')}


def load_verified_sources(root):
    root = Path(root).resolve()
    source_manifest = root / 'manifests/chembl_downloads.json'
    doc_manifest = root / 'manifests/document_downloads.json'
    records = json.loads(source_manifest.read_text())
    doc_records = json.loads(doc_manifest.read_text())
    pages = collections.defaultdict(list)
    activities, assays, documents = [], {}, {}
    inputs = [file_record(source_manifest, root), file_record(doc_manifest, root)]
    seen_paths = set()
    for rec in [*records, *doc_records]:
        require(rec['path'] not in seen_paths, f'Duplicate manifest path: {rec["path"]}')
        seen_paths.add(rec['path'])
        path = safe_path(root, rec['path'])
        info = file_record(path, root)
        require(info['sha256'] == rec['sha256'], f'Hash mismatch: {rec["path"]}')
        require('bytes' not in rec or info['bytes'] == rec['bytes'], f'Size mismatch: {rec["path"]}')
        inputs.append(info)
        data = json.loads(path.read_text())
        resource = next((name for name in ['activities', 'assays', 'documents'] if name in data), None)
        require(resource is not None, f'Unknown source response: {rec["path"]}')
        rows = data[resource]
        require(len(rows) == rec['rows'], f'Manifest row mismatch: {rec["path"]}')
        meta = data['page_meta']
        if resource == 'documents':
            requested = set(rec['params']['document_chembl_id__in'].split(','))
            returned = [r['document_chembl_id'] for r in rows]
            require(len(returned) == len(set(returned)) and set(returned) == requested,
                    f'Incomplete/overlapping document batch: {rec["path"]}')
            require(meta['total_count'] == len(rows), f'Incomplete document pagination: {rec["path"]}')
            for row in rows:
                require(row['document_chembl_id'] not in documents, 'Document batches overlap')
                documents[row['document_chembl_id']] = row
        else:
            target = rec['params']['target_chembl_id']
            require(meta['total_count'] == rec['total_count'], f'Page total mismatch: {rec["path"]}')
            require(meta['offset'] == rec['params']['offset'], f'Offset mismatch: {rec["path"]}')
            require(meta['limit'] == rec['params']['limit'], f'Limit mismatch: {rec["path"]}')
            require(all(r['target_chembl_id'] == target for r in rows), f'Target mismatch: {rec["path"]}')
            pages[(resource, target)].append((rec, rows))
            if resource == 'activities':
                activities.extend(rows)
            else:
                for row in rows:
                    require(row['assay_chembl_id'] not in assays, 'Assay pages overlap')
                    assays[row['assay_chembl_id']] = row
    pagination = {}
    for (resource, target), grouped in sorted(pages.items()):
        grouped.sort(key=lambda pair: pair[0]['params']['offset'])
        totals = {rec['total_count'] for rec, _ in grouped}
        limits = {rec['params']['limit'] for rec, _ in grouped}
        require(len(totals) == len(limits) == 1, f'Inconsistent pagination: {resource}/{target}')
        total, limit = next(iter(totals)), next(iter(limits))
        offsets = [rec['params']['offset'] for rec, _ in grouped]
        require(offsets == list(range(0, total, limit)), f'Missing/duplicate page offsets: {resource}/{target}')
        for rec, rows in grouped:
            require(len(rows) == min(limit, total - rec['params']['offset']), f'Partial page: {rec["path"]}')
        require(sum(len(rows) for _, rows in grouped) == total, 'Incomplete total')
        pagination[f'{resource}/{target}'] = {'pages': len(grouped), 'rows': total, 'offsets': offsets}
    ids = [str(r['activity_id']) for r in activities]
    require(len(ids) == len(set(ids)), 'Activity pages overlap')
    for row in activities:
        require(row['assay_chembl_id'] in assays, f'Missing assay join: {row["activity_id"]}')
        require(row['document_chembl_id'] in documents, f'Missing document join: {row["activity_id"]}')
        assay = assays[row['assay_chembl_id']]
        require(assay['target_chembl_id'] == row['target_chembl_id'], 'Activity/assay target mismatch')
        require(assay['document_chembl_id'] == row['document_chembl_id'], 'Activity/assay document mismatch')
    return activities, assays, documents, {'pagination': pagination, 'input_files': inputs,
        'activities': len(activities), 'assays': len(assays), 'documents': len(documents),
        'duplicate_activity_ids': 0, 'missing_assay_joins': 0, 'missing_document_joins': 0,
        'activity_assay_target_or_document_mismatches': 0}


def validate_ledger(root, activities, assays):
    root = Path(root).resolve()
    ledger_path = root / 'reports/activity_ledger.csv'
    rows = list(csv.DictReader(ledger_path.open()))
    by_id = {str(r['activity_id']): r for r in activities}
    ledger_ids = [r['activity_id'] for r in rows]
    require(len(ledger_ids) == len(set(ledger_ids)) and set(ledger_ids) == set(by_id),
            'Ledger activity coverage/uniqueness mismatch')
    raw_fields = ['molecule_chembl_id', 'parent_molecule_chembl_id', 'target_chembl_id',
        'assay_chembl_id', 'document_chembl_id', 'document_year', 'standard_type',
        'standard_relation', 'standard_units', 'standard_value', 'data_validity_comment',
        'potential_duplicate', 'src_id']
    selected = []
    for row in rows:
        original = by_id[row['activity_id']]
        for field in raw_fields:
            expected = '' if original.get(field) is None else str(original[field])
            require(row[field] == expected, f'Ledger field changed: {row["activity_id"]}/{field}')
        valid_structure = not row['structure_error']
        require(valid_structure == bool(row['canonical_smiles'] and row['scaffold_group']),
                f'Inconsistent structure audit: {row["activity_id"]}')
        reasons = strict_reasons(original, assays[row['assay_chembl_id']], row['structure_error'])
        require(row['exclusion_reasons'].split(';') == (reasons or ['']), 'Strict-filter reasons mismatch')
        require(row['candidate_exact_biochemical'] == str(not reasons), 'Strict inclusion mismatch')
        if not reasons:
            selected.append(row)
    inventory_path = root / 'reports/inventory.json'
    inventory = json.loads(inventory_path.read_text())
    require(counts(rows) == inventory['all'], 'Inventory all-count mismatch')
    require(counts(selected) == inventory['strict_candidate'], 'Inventory strict-count mismatch')
    return rows, {'all': counts(rows), 'strict_metadata_candidates': counts(selected),
        'raw_field_mismatches': 0, 'filter_mismatches': 0,
        'structure_standardization_recomputed': False,
        'limitation': 'Canonical parents/scaffolds are trusted from existing ledger; no independent RDKit standardization rerun.',
        'input_files': [file_record(ledger_path, root), file_record(inventory_path, root)]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    activities, assays, documents, source_report = load_verified_sources(args.root)
    _, ledger_report = validate_ledger(args.root, activities, assays)
    report = {'status': 'PASS_ENGINEERING_SOURCE_AND_FILTER_AUDIT',
        'verified_utc': datetime.now(timezone.utc).isoformat(),
        'source_validation': source_report, 'ledger_validation': ledger_report,
        'script': file_record(Path(__file__)), 'models_trained': False,
        'scientific_comparability_confirmed': False}
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / 'validation.json').write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')
    print(json.dumps({k: v for k, v in report.items() if k not in ['source_validation', 'ledger_validation']}, indent=2))


if __name__ == '__main__':
    main()
