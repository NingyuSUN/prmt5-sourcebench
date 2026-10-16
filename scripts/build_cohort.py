#!/usr/bin/env python3
"""Freeze source-defined PRMT5 binary labels and audit source/structure separation."""
import argparse
import collections
import csv
import hashlib
import importlib.util
import json
import os
import tempfile
from datetime import date
from pathlib import Path
import numpy as np
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator

UPSTREAM_SHA = '0a059a1ec5528ff140ad5634a277f6fb96022acbfec4ab9aa8f3bb450698627b'
TRAIN_FAMILIES = {'50001259', '50029211', '50001260'}
CUTOFF = '2019-09-01'


def require(condition, message):
    if not condition:
        raise ValueError(message)


def interval_label(lower, upper, lower_inclusive=True, threshold=100):
    """Positive only if all positive possible values <= threshold; negative if all >."""
    require(lower is None or lower >= 0, 'Negative interval boundary')
    require(upper is None or upper > 0, 'Nonpositive interval')
    require(lower is None or upper is None or lower <= upper, 'Reversed interval')
    if upper is not None and upper <= threshold:
        return 1
    if lower is not None and (lower > threshold or (lower == threshold and not lower_inclusive)):
        return 0
    return None


def old_label(row):
    lower = float(row['lower_nM']) if row['lower_nM'] else None
    upper = float(row['upper_nM']) if row['upper_nM'] else None
    return interval_label(lower, upper, row['bound_semantics'] not in ('>', '()','(]'))


def later_label(value):
    if value == '+++':
        return interval_label(1, 100)
    if value == '++':
        return interval_label(100, 1000, False)
    if value == '+':
        return interval_label(1000, 10000, False)
    value = float(value)
    require(value > 0, 'Nonpositive exact value')
    return interval_label(value, value)



def later_source(assay_id):
    sources = {'CHEMBL5739093': ('validation', '67846363'),
               'CHEMBL5737316': ('test', '75338422')}
    require(assay_id in sources, 'Unknown later assay: ' + assay_id)
    return sources[assay_id]


def admit_training_label(row):
    # Source/date admission happens before inspecting its activity label.
    if row['family'] not in TRAIN_FAMILIES:
        return None, 'unresolved_family_identity_mapping'
    require(date.fromisoformat(row['direct_source_publication_date']) < date.fromisoformat(CUTOFF),
            'Post-cutoff training source')
    if row['original_bin_mapped'] != 'True' or not row['original_bin']:
        return None, 'ambiguous_original_bin'
    label = old_label(row)
    require(label == int(row['binary_active_at_le100nM']), 'Interval-derived label disagreement')
    require(not row['structure_error'] and bool(row['canonical_smiles']), 'Invalid training structure')
    return label, None


def quality_evidence(activity_ids, ledger):
    flags = []
    for activity_id in sorted(set(activity_ids.split(';'))):
        require(activity_id in ledger, 'Missing activity metadata: ' + activity_id)
        source = ledger[activity_id]
        require(not source['data_validity_comment'], 'Unreviewed activity validity flag: ' + activity_id)
        flags.append({'activity_id': activity_id,
                      'data_validity_comment': source['data_validity_comment'],
                      'potential_duplicate': source['potential_duplicate']})
    return flags


def overlap_counts(reference, query):
    counts = {}
    for field in ('canonical_smiles', 'scaffold_group', 'source_family'):
        def values(rows):
            return {part for row in rows for part in
                    (row[field].split(';') if field == 'source_family' else [row[field]])}
        counts[field] = len(values(reference) & values(query))
    return counts


def atomic_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode='w', encoding='utf-8', dir=path.parent,
                                         prefix='.' + path.name + '.', delete=False) as handle:
            temporary = Path(handle.name)
            handle.write(json.dumps(value, indent=2) + '\n')
            handle.flush()
            os.fsync(handle.fileno())
        os.replace(temporary, path)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def publish_outputs(root, output, artifacts, report_path, report):
    # All cohort validation must pass before calling this function.
    # Each file is replaced atomically. During multi-file promotion no PASS is exposed.
    output.mkdir(parents=True, exist_ok=True)
    staged = []
    try:
        for name, rows in artifacts:
            require(bool(rows), 'Empty artifact: ' + name)
            with tempfile.NamedTemporaryFile(mode='w', newline='', encoding='utf-8', dir=output,
                                             prefix='.' + name + '.', delete=False) as handle:
                temporary = Path(handle.name)
                staged.append((temporary, output / name))
                writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
                writer.writeheader()
                writer.writerows(rows)
                handle.flush()
                os.fsync(handle.fileno())
        report['outputs'] = {str(target.relative_to(root)): sha(temporary)
                             for temporary, target in staged}
        atomic_json(report_path, {'status': 'BUILD_IN_PROGRESS',
                                  'reason': 'Artifact promotion; PASS is published last.'})
        for temporary, target in staged:
            os.replace(temporary, target)
        for _, target in staged:
            require(sha(target) == report['outputs'][str(target.relative_to(root))],
                    'Published artifact hash mismatch: ' + target.name)
        atomic_json(report_path, report)
    finally:
        for temporary, _ in staged:
            temporary.unlink(missing_ok=True)


def collapse(rows):
    groups = collections.defaultdict(list)
    for row in rows:
        groups[row['canonical_smiles']].append(row)
    included, conflicts = [], []
    for parent, group in sorted(groups.items()):
        if len({row['binary_y'] for row in group}) != 1:
            conflicts.extend(group)
            continue
        first = group[0]
        require(len({r['scaffold_group'] for r in group}) == 1, 'Inconsistent parent scaffold')
        included.append({
            'compound_id': hashlib.sha256(parent.encode()).hexdigest()[:16],
            'canonical_smiles': parent, 'scaffold_group': first['scaffold_group'],
            'binary_y': first['binary_y'], 'split': first['split'],
            'source_family': ';'.join(sorted({r['source_family'] for r in group})),
            'assay_ids': ';'.join(sorted({r['assay_ids'] for r in group})),
            'source_date_bound': min(r['source_date_bound'] for r in group),
            'activity_ids': ';'.join(sorted({x for r in group for x in r['activity_ids'].split(';')})),
            'source_evidence': json.dumps([r['source_evidence'] for r in group], sort_keys=True),
        })
    return included, conflicts


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def readcsv(path):
    with path.open(newline='') as handle:
        return list(csv.DictReader(handle))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', required=True, type=Path)
    parser.add_argument('--chemistry', required=True, type=Path)
    parser.add_argument('--epizyme', type=Path)
    args = parser.parse_args()
    root = args.root.resolve()
    require(sha(args.chemistry) == UPSTREAM_SHA, 'Upstream chemistry hash mismatch')
    spec = importlib.util.spec_from_file_location('pinned_chemistry', args.chemistry)
    chemistry = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(chemistry)
    early_path = args.epizyme or root / 'data/curation/epizyme/results/epizyme_activity_curation.csv'
    later_path = root / 'data/curation/later/binary_parent_candidates.csv'
    early, later = readcsv(early_path), readcsv(later_path)
    ledger_path = root / 'reports/activity_ledger.csv'
    ledger = {row['activity_id']: row for row in readcsv(ledger_path)}
    train_rows, exclusions = [], []
    for row in early:
        label, reason = admit_training_label(row)
        if reason:
            exclusions.append({'activity_id': row['activity_id'], 'reason': reason})
            continue
        train_rows.append({
            'canonical_smiles': row['canonical_smiles'], 'scaffold_group': row['scaffold_group'],
            'binary_y': label, 'split': 'train', 'source_family': row['family'],
            'assay_ids': row['assay_chembl_id'], 'source_date_bound': row['direct_source_publication_date'],
            'activity_ids': row['activity_id'], 'source_evidence': {
                'patent': row['patent'], 'examples': row['source_compound_numbers'],
                'quality_flags': quality_evidence(row['activity_id'], ledger),
                'original_bin': row['original_bin'], 'lower_nM': row['lower_nM'],
                'upper_nM': row['upper_nM'], 'bound_semantics': row['bound_semantics'],
                'family_first_publication': row['WO_publication_date'],
                'pair_first_publication_verified': False,
                'structure_provenance': 'ChEMBL compound-record crosswalk; original figure not independently decoded'},
        })
    train, conflicts = collapse(train_rows)
    require(not conflicts, 'Unexpected pre-cutoff parent label conflict')
    split_rows = {'train': train, 'validation': [], 'test': []}
    for row in later:
        split, family = later_source(row['assay_id'])
        if row['decision'] != 'candidate_include':
            exclusions.append({'activity_id': row['activity_ids'], 'reason': row['decision']})
            continue
        values = json.loads(row['source_values_original'])
        labels = {later_label(v) for v in values.values()}
        require(labels == {int(row['binary_y'])}, 'Unresolved later source-bin conflict')
        require(date.fromisoformat(row['source_date_verified']) >= date.fromisoformat(CUTOFF),
                'Unexpected later date')
        split_rows[split].append({
            'canonical_smiles': row['canonical_smiles'], 'scaffold_group': row['scaffold_group'],
            'binary_y': int(row['binary_y']), 'split': split, 'source_family': family,
            'assay_ids': row['assay_id'], 'source_date_bound': row['source_date_verified'],
            'activity_ids': row['activity_ids'], 'source_evidence': {
                'examples': row['source_examples'], 'values': values,
                'quality_flags': quality_evidence(row['activity_ids'], ledger),
                'precision': row['numeric_precision'], 'date_type': row['source_date_type']},
        })
    for split in ('validation', 'test'):
        split_rows[split], conflicts = collapse(split_rows[split])
        require(not conflicts, 'Unexpected later conflict')
    all_rows = [r for group in split_rows.values() for r in group]
    for row in all_rows:
        checked = chemistry.standardize_smiles(row['canonical_smiles'])
        require(checked['canonical_smiles'] == row['canonical_smiles'], 'Parent not stable')
        require(checked['scaffold_group'] == row['scaffold_group'], 'Scaffold differs from pinned upstream')
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048, includeChirality=True)
    fps = {k: [generator.GetFingerprint(Chem.MolFromSmiles(r['canonical_smiles'])) for r in v]
           for k, v in split_rows.items()}
    overlaps, neighbors = {}, []
    for ref, query in [('train','validation'), ('train','test'), ('validation','test')]:
        key = ref + '_vs_' + query
        counts = overlap_counts(split_rows[ref], split_rows[query])
        require(not any(counts.values()), f'Leakage: {key}: {counts}')
        maxima = []
        for row, fp in zip(split_rows[query], fps[query]):
            sims = np.array(DataStructs.BulkTanimotoSimilarity(fp, fps[ref]))
            best = int(np.argmax(sims)); maximum = float(sims[best]); maxima.append(maximum)
            neighbors.append({'reference_split': ref, 'query_split': query,
                'query_id': row['compound_id'], 'nearest_id': split_rows[ref][best]['compound_id'],
                'max_tanimoto': maximum})
        overlaps[key] = {'overlap_counts': counts, 'query_count': len(maxima),
            'maximum_similarity_quantiles': dict(zip(['min','p25','median','p75','max'],
                np.quantile(maxima,[0,.25,.5,.75,1]).tolist())),
            'query_ge_0_8': sum(x >= .8 for x in maxima), 'query_ge_0_9': sum(x >= .9 for x in maxima)}
    output = root / 'data/processed'
    counts = {k: {'parents':len(v), 'positive':sum(r['binary_y'] for r in v),
                  'negative':sum(1-r['binary_y'] for r in v),
                  'scaffolds':len({r['scaffold_group'] for r in v}),
                  'families':sorted({family for r in v for family in r['source_family'].split(';')}),
                  'source_date_bounds': sorted({r['source_date_bound'] for r in v})}
              for k,v in split_rows.items()}
    expected={'train':(604,164), 'validation':(70,39), 'test':(246,69)}
    for split,(n,positive) in expected.items():
        require((counts[split]['parents'],counts[split]['positive']) == (n,positive), 'Frozen capacity changed')
    report={'status':'PASS_FOR_NARROW_SOURCE_SHIFT_BINARY_BENCHMARK',
        'continuous_potency_regression':'NO_GO_original_intervals_and_protocol_incompatibility',
        'threshold_nM':100, 'train_source_publication_cutoff':CUTOFF,
        'counts':counts, 'overlap_and_similarity':overlaps,
        'excluded_records_by_reason':dict(collections.Counter(r['reason'] for r in exclusions)),
        'training_source_rows':len(train_rows), 'final_parents':len(all_rows),
        'source_record_dedup_policy': {
            'unit': 'One standardized parent per split, never activity-row counts as sample counts.',
            'same_parent_same_label': 'Collapse while preserving all activity IDs, source evidence and quality flags.',
            'same_parent_conflicting_labels': 'Quarantine, never average or select a favorable label.',
            'potential_duplicate': 'Retain only with reviewed original-table label; preserve flag; parent collapse prevents row inflation.',
            'validity_warning': 'Reject until explicitly reviewed; no included warnings in this frozen cohort.'},
        'training_potential_duplicate_source_rows':sum(
            any(flag['potential_duplicate'] == '1' for flag in row['source_evidence']['quality_flags'])
            for row in train_rows),
        'fingerprint':{'radius':2,'bits':2048,'include_chirality':True},
        'limitations':['Not exact first-ever molecule-label disclosure dates; documented public-date bounds.',
            'Source, protocol, time and chemical series are confounded; one final test family.',
            'Shared catalytic target does not establish equal assay scale.',
            'Original chemical figures not all independently decoded; structures rely on curated database crosswalk.',
            'No prospective, experimental, cellular or clinical validation.',
            'No performance result or model fitting is part of this data gate.'],
        'inputs':{str(early_path.relative_to(root)) if early_path.is_relative_to(root) else 'epizyme_activity_curation.csv':sha(early_path),
                  str(later_path.relative_to(root)):sha(later_path),
                  str(ledger_path.relative_to(root)):sha(ledger_path)},
        'reused_chemistry_sha256':UPSTREAM_SHA}
    publish_outputs(root, output,
        [('cohort.csv', all_rows), ('nearest_neighbors.csv', neighbors), ('exclusions.csv', exclusions)],
        root / 'reports/curation/cohort_gate.json', report)
    print(json.dumps({'status':report['status'],'counts':counts,'overlap_and_similarity':overlaps},indent=2))

if __name__=='__main__':
    main()
