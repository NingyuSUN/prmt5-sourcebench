"""Activity-independent metadata/structure feasibility audit; no model fitting."""
import argparse, collections, csv, hashlib, json, math, re, sys
from pathlib import Path
from datetime import datetime, timezone

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--root', type=Path, required=True)
    p.add_argument('--reuse-src', type=Path, required=True)
    a = p.parse_args()
    sys.path.insert(0, str(a.reuse_src))
    import rdkit
    from rdkit import Chem, RDLogger
    from drugperturb.chemistry import standardize_smiles
    RDLogger.DisableLog('rdApp.warning')
    RDLogger.DisableLog('rdApp.error')
    expected = '0a059a1ec5528ff140ad5634a277f6fb96022acbfec4ab9aa8f3bb450698627b'
    if hashlib.sha256((a.reuse_src / 'drugperturb/chemistry.py').read_bytes()).hexdigest() != expected:
        raise RuntimeError('Pinned chemistry hash mismatch')
    records = json.loads((a.root / 'manifests/chembl_downloads.json').read_text())
    rows = []
    assays = {}
    for rec in records:
        path = a.root / rec['path']
        data = path.read_bytes()
        if hashlib.sha256(data).hexdigest() != rec['sha256']:
            raise RuntimeError('Source hash mismatch')
        d = json.loads(data)
        rows.extend(d.get('activities', []))
        assays.update({x['assay_chembl_id']: x for x in d.get('assays', [])})
    if len(rows) != len({x['activity_id'] for x in rows}):
        raise RuntimeError('Overlapping activity pages')
    for t in {r['target_chembl_id'] for r in rows}:
        n = sum((r['target_chembl_id'] == t for r in rows))
        expected_n = next((x['total_count'] for x in records if x['params']['target_chembl_id'] == t and '/activity.' in x['url']))
        if n != expected_n:
            raise RuntimeError('Incomplete pagination')
    cache = {}
    output = []
    for r in rows:
        assay = assays[r['assay_chembl_id']]
        sm = r.get('canonical_smiles')
        if sm not in cache:
            try:
                cache[sm] = standardize_smiles(sm or '')
            except Exception as ex:
                cache[sm] = {'error': type(ex).__name__ + ': ' + str(ex)}
        st = cache[sm]
        desc = assay['description']
        flags = []
        if 'error' in st:
            flags.append('structure_invalid')
        if r.get('standard_type') != 'IC50':
            flags.append('not_IC50')
        if r.get('standard_relation') != '=':
            flags.append('not_exact')
        if r.get('standard_units') != 'nM':
            flags.append('not_nM')
        try:
            val = float(r.get('standard_value'))
            good = math.isfinite(val) and val > 0
        except (ValueError, TypeError):
            good = False
        if not good:
            flags.append('nonpositive_or_missing')
        if r.get('data_validity_comment'):
            flags.append('validity_flag')
        if r.get('potential_duplicate'):
            flags.append('potential_duplicate')
        if assay.get('confidence_score') not in [7, 9]:
            flags.append('not_direct_target')
        if assay.get('bao_label') not in ['single protein format', 'protein complex format'] or assay.get('assay_type') != 'B':
            flags.append('not_biochemical_format')
        if r.get('target_tax_id') != '9606' and str(r.get('target_tax_id')) != '9606':
            flags.append('target_not_human')
        o = {k: r.get(k) for k in ['activity_id', 'molecule_chembl_id', 'parent_molecule_chembl_id', 'target_chembl_id', 'assay_chembl_id', 'document_chembl_id', 'document_year', 'standard_type', 'standard_relation', 'standard_units', 'standard_value', 'data_validity_comment', 'potential_duplicate', 'src_id']}
        o.update(canonical_smiles=st.get('canonical_smiles', ''), scaffold_group=st.get('scaffold_group', ''), structure_error=st.get('error', ''), exclusion_reasons=';'.join(flags), candidate_exact_biochemical=not flags, assay_format=assay['bao_label'], confidence=assay['confidence_score'], assay_description=desc, mentions_MTA=bool(re.search('\\bMTA\\b|methylthioadenosine', desc, re.I)), mentions_SAM=bool(re.search('\\bSAM\\b|adenosyl', desc, re.I)))
        output.append(o)

    def counts(rs):
        return dict(rows=len(rs), molecules=len({r['molecule_chembl_id'] for r in rs}), parents=len({r['canonical_smiles'] for r in rs if r['canonical_smiles']}), scaffolds=len({r['scaffold_group'] for r in rs if r['scaffold_group']}), assays=len({r['assay_chembl_id'] for r in rs}), documents=len({r['document_chembl_id'] for r in rs}))

    def hist(rs, key):
        return dict(sorted(collections.Counter((str(r.get(key)) for r in rs)).items()))
    selected = [r for r in output if r['candidate_exact_biochemical']]
    assay_summ = []
    for aid, rs in __import__('itertools').groupby(sorted(selected, key=lambda r: r['assay_chembl_id']), lambda r: r['assay_chembl_id']):
        rs = list(rs)
        assay_summ.append(dict(assay_id=aid, **counts(rs), document_id=rs[0]['document_chembl_id'], year=rs[0]['document_year'], mentions_MTA=rs[0]['mentions_MTA'], description=rs[0]['assay_description']))
    assay_summ.sort(key=lambda r: r['parents'], reverse=True)
    first = {}
    for r in output:
        if r['canonical_smiles'] and r['document_year']:
            key = r['canonical_smiles']
            first[key] = min(first.get(key, 9999), int(r['document_year']))
    temporal = []
    for cutoff in [2018, 2019, 2020, 2021, 2022, 2023, 2024]:
        earlier = [r for r in selected if r['document_year'] and int(r['document_year']) <= cutoff]
        later = [r for r in selected if r['document_year'] and int(r['document_year']) > cutoff]
        novel = [r for r in later if first[r['canonical_smiles']] > cutoff]
        seen_sc = {r['scaffold_group'] for r in earlier}
        both = [r for r in novel if r['scaffold_group'] not in seen_sc]
        temporal.append(dict(cutoff_year=cutoff, earlier=counts(earlier), later=counts(later), later_new_parent=counts(novel), later_new_parent_and_scaffold=counts(both)))
    result = dict(created_utc=datetime.now(timezone.utc).isoformat(), chemistry_sha256=expected, rdkit_version=rdkit.__version__, all=counts(output), endpoint_counts=hist(output, 'standard_type'), format_counts=hist(output, 'assay_format'), years=hist(output, 'document_year'), strict_candidate=counts(selected), strict_years=hist(selected, 'document_year'), strict_source_counts=hist(selected, 'src_id'), strict_MTA_mention=hist(selected, 'mentions_MTA'), exclusions=dict(collections.Counter((x for r in output for x in r['exclusion_reasons'].split(';') if x))), temporal_year_proxies=temporal, largest_assays=assay_summ, limitations=['Candidate biochemical != confirmed comparable conditions', 'Year is document year not earliest public disclosure', 'No cross-source or patent-family de-duplication yet', 'MTA keyword presence is not known concentration and no mention does not mean absence', 'No training or predictive evaluation performed'])
    (a.root / 'reports/inventory.json').write_text(json.dumps(result, indent=2) + '\n')
    with (a.root / 'reports/activity_ledger.csv').open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(output[0]))
        w.writeheader()
        w.writerows(output)
    with (a.root / 'reports/assay_inventory.csv').open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(assay_summ[0]))
        w.writeheader()
        w.writerows(assay_summ)
    print(json.dumps({k: v for k, v in result.items() if k != 'largest_assays'}, indent=2))
    print('TOP_ASSAYS')
    print(json.dumps(assay_summ[:15], indent=2))
if __name__ == '__main__':
    main()
