#!/usr/bin/env python3
"""Audit predefined candidate assay counts/conflicts; no fitting or new chemistry."""
import argparse
import collections
import csv
import json
from datetime import datetime, timezone
from pathlib import Path
from validate_sources import (counts, exact_positive_nm_ic50, file_record,
                              load_verified_sources, validate_ledger, require)

# Explicit metadata-reviewed list: the separate 75 nM *peptide*/15 uM SAM assay
# CHEMBL5737229 is intentionally excluded. This is NOT a trained/accepted cohort.
SAM75_ASSAYS = ('CHEMBL5731839', 'CHEMBL5732555', 'CHEMBL5732693',
                'CHEMBL5733925', 'CHEMBL5734324', 'CHEMBL5736675')
MTA_ASSAYS = ('CHEMBL5738640', 'CHEMBL5738688', 'CHEMBL5738641', 'CHEMBL5738689')


def parent_conflicts(rows):
    grouped = collections.defaultdict(list)
    for row in rows:
        grouped[row['canonical_smiles']].append(row)
    conflicts = []
    for parent, entries in sorted(grouped.items()):
        values = sorted({float(row['standard_value']) for row in entries})
        if len(values) > 1:
            conflicts.append({'canonical_smiles': parent,
                'activity_ids': sorted(row['activity_id'] for row in entries),
                'assay_ids': sorted({row['assay_chembl_id'] for row in entries}),
                'document_ids': sorted({row['document_chembl_id'] for row in entries}),
                'values_nM': values, 'fold_spread': max(values) / min(values)})
    return {'counts': counts(rows),
        'repeated_parents': sum(len(entries) > 1 for entries in grouped.values()),
        'distinct_value_parents': len(conflicts),
        'maximum_fold_spread': max((r['fold_spread'] for r in conflicts), default=None),
        'conflicts': conflicts}


def novel_against(train, query, first_year=None, cutoff=None):
    train_parents = {r['canonical_smiles'] for r in train}
    train_groups = {r['scaffold_group'] for r in train}
    return [r for r in query if r['canonical_smiles'] not in train_parents
            and r['scaffold_group'] not in train_groups
            and (cutoff is None or first_year.get(r['canonical_smiles'], 0) > cutoff)]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    root = args.root.resolve()
    activities, assays, documents, source_check = load_verified_sources(root)
    ledger, ledger_check = validate_ledger(root, activities, assays)
    raw = {str(r['activity_id']): r for r in activities}
    eligible = [r for r in ledger if r['canonical_smiles'] and exact_positive_nm_ic50(r)]
    def quality_ok(row):
        source = raw[row['activity_id']]
        return not source.get('data_validity_comment') and not source.get('potential_duplicate')
    mta = {}
    for assay_id in MTA_ASSAYS:
        selected = [r for r in eligible if r['assay_chembl_id'] == assay_id]
        mta[assay_id] = parent_conflicts(selected)
        mta[assay_id]['quality_unflagged_counts'] = counts([r for r in selected if quality_ok(r)])
        validity_unflagged = [r for r in selected if not raw[r['activity_id']].get('data_validity_comment')]
        mta[assay_id]['validity_unflagged_retaining_duplicate_flags'] = parent_conflicts(validity_unflagged)
        mta[assay_id]['source_assay_description'] = assays[assay_id]['description']
    cluster = [r for r in eligible if r['assay_chembl_id'] in SAM75_ASSAYS]
    require(len({r['assay_chembl_id'] for r in cluster}) == 6, 'Missing predefined candidate assay')
    for assay_id in SAM75_ASSAYS:
        desc = assays[assay_id]['description']
        require('3H-SAM was 75 nM' in desc and 'PRMT5/MEP50 was 4 nM' in desc,
                f'Candidate protocol text changed: {assay_id}')
    per_assay = []
    for assay_id in SAM75_ASSAYS:
        selected = [r for r in cluster if r['assay_chembl_id'] == assay_id]
        assay = assays[assay_id]
        doc = documents[assay['document_chembl_id']]
        result = parent_conflicts(selected)
        result.update(assay_id=assay_id, document_id=doc['document_chembl_id'],
            document_year=doc['year'], patent_id=doc['patent_id'],
            confidence_score=assay['confidence_score'],
            quality_unflagged_counts=counts([r for r in selected if quality_ok(r)]))
        per_assay.append(result)
    # This compares local document IDs, not legally verified patent families.
    source_holdouts = []
    for doc_id in sorted({r['document_chembl_id'] for r in cluster}):
        train = [r for r in cluster if r['document_chembl_id'] != doc_id]
        query = [r for r in cluster if r['document_chembl_id'] == doc_id]
        source_holdouts.append({'document_id': doc_id, 'train': counts(train),
            'query': counts(query), 'query_parent_and_scaffold_absent_from_train':
                counts(novel_against(train, query))})
    first_year = {}
    for row in ledger:
        if row['canonical_smiles'] and row['document_year']:
            parent = row['canonical_smiles']
            first_year[parent] = min(first_year.get(parent, 9999), int(row['document_year']))
    temporal = []
    for cutoff in sorted({int(r['document_year']) for r in cluster}):
        earlier = [r for r in cluster if int(r['document_year']) <= cutoff]
        later = [r for r in cluster if int(r['document_year']) > cutoff]
        temporal.append({'cutoff_year': cutoff, 'earlier': counts(earlier), 'later': counts(later),
            'later_new_dated_parent_and_training_scaffold':
                counts(novel_against(earlier, later, first_year, cutoff))})
    combined = parent_conflicts(cluster)
    within_count = sum(r['distinct_value_parents'] for r in per_assay)
    report = {'status': 'NOT_READY_FOR_TRAINING_PROTOCOL_AND_PROVENANCE_REVIEW_REQUIRED',
        'verified_utc': datetime.now(timezone.utc).isoformat(),
        'source_and_filter_checks': {'status': 'PASS_ENGINEERING_ONLY',
            'source_counts': {k: source_check[k] for k in ['activities', 'assays', 'documents']},
            'ledger_counts': ledger_check['all'], 'strict_candidates': ledger_check['strict_metadata_candidates']},
        'mta_assays': mta,
        'sam75_metadata_protocol_candidate': {'assay_ids': list(SAM75_ASSAYS),
            'exact_nm_ic50_valid_structure_counts': counts(cluster),
            'quality_unflagged_counts': counts([r for r in cluster if quality_ok(r)]),
            'confidence_scores': dict(collections.Counter(str(assays[a]['confidence_score']) for a in SAM75_ASSAYS)),
            'per_assay': per_assay, 'pooled_conflict_inventory': combined,
            'within_assay_conflicting_parent_keys': within_count,
            'source_holdout_capacity_only': source_holdouts,
            'document_year_proxy_capacity_only': temporal,
            'biological_comparability_confirmed': False,
            'patent_family_independence_confirmed': False},
        'limitations': ['Exact-value disagreement is a provenance-review trigger; not automatically an assay error.',
            'No averaging or selecting among discordant values was performed.',
            'Scaffolds/canonical parents reused from ledger; standardization not independently recomputed.',
            'Source documents are not independent patent families.',
            'Document years are not earliest disclosure or patent priority dates.',
            'Source/scaffold capacities are exploratory; no split or outcome threshold selected.',
            'Strict 7/9 confidence gate excludes the six confidence-6 candidate assays; manual target review is required before relaxing it.',
            'Matching description fragments do not establish matched full assay conditions.',
            'No training, prediction, test-set performance or biological validation performed.'],
        'input_manifest_paths': ['manifests/chembl_downloads.json', 'manifests/document_downloads.json'],
        'input_files': source_check['input_files'] + ledger_check['input_files'],
        'script_files': [file_record(Path(__file__)), file_record(Path(__file__).with_name('validate_sources.py'))]}
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / 'cohort_audit.json').write_text(json.dumps(report, indent=2, allow_nan=False) + '\n')
    lines = ['# PRMT5 independent cohort audit', '', '**NOT READY for training.** Engineering checks passed; assay comparability and patent-family independence remain unresolved.', '',
        f"Verified {len(activities):,} source activities, {len(assays)} assays, {len(documents)} documents; {ledger_check['strict_metadata_candidates']['rows']} strict metadata candidate rows.", '',
        '## Within-assay exact IC50 disagreements', '',
        '| Assay | Selection denominator | Exact rows | Parents | Repeated parents | Parents with distinct values | Maximum spread |',
        '|---|---|---:|---:|---:|---:|---:|']
    for assay_id, result in mta.items():
        c = result['counts']
        lines.append(f"| {assay_id} | All exact positive nM IC50 with valid structure | {c['rows']} | {c['parents']} | {result['repeated_parents']} | {result['distinct_value_parents']} | {result['maximum_fold_spread']} |")
        checked = result['validity_unflagged_retaining_duplicate_flags']
        v = checked['counts']
        lines.append(f"| {assay_id} | No data-validity comment; potential-duplicate flags retained | {v['rows']} | {v['parents']} | {checked['repeated_parents']} | {checked['distinct_value_parents']} | {checked['maximum_fold_spread']} |")
    c, q = counts(cluster), counts([r for r in cluster if quality_ok(r)])
    lines += ['', '## 75 nM SAM candidate protocol cluster', '',
        f"Six explicitly selected assays: {c['rows']} exact positive nM IC50 rows, {c['parents']} canonical parents, {c['scaffolds']} scaffold groups. Quality-unflagged subset: {q['rows']} rows / {q['parents']} parents / {q['scaffolds']} groups.", '',
        f"All six have database confidence 6; none is included by the original 7/9 strict gate. Across this pooled cluster, {combined['distinct_value_parents']} parents have distinct exact values; {within_count} parent-assay keys disagree within one assay. Pooled differences may reflect different contexts and must not be labeled errors automatically.", '',
        'The separate CHEMBL5737229 uses 75 nM peptide and 15 uM SAM; it is excluded from this cluster.', '',
        'Source and document-year holdout tables in cohort_audit.json are capacity inventories. No comparability-confirmed, patent-family-independent split has been established.', '',
        '## Scope', '', *['- ' + x for x in report['limitations']], '',
        'Inputs are recorded as project-relative paths and hashes in cohort_audit.json. Existing source bytes and ledger were read only.']
    (args.output_dir / 'cohort_audit.md').write_text('\n'.join(lines) + '\n')
    print(json.dumps({'status': report['status'], 'sam75_counts': c,
                      'sam75_quality_unflagged': q, 'mta_conflicting_parents':
                        {a: r['distinct_value_parents'] for a, r in mta.items()}}, indent=2))


if __name__ == '__main__':
    main()
