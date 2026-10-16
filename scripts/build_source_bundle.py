#!/usr/bin/env python3
"""Build or verify/replay the scoped offline PRMT5 source snapshot; no network or fitting."""
from __future__ import annotations
import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import io
import json
import os
from pathlib import Path, PurePosixPath
import subprocess
import sys
import tarfile
import tempfile

PROTOCOL_SHA = '283d38f275ef4c91d17dc87f1049213936543353d873198ba219b9a0b5890d1f'
COHORT_SHA = '17772cc5f360464abb2d5b913f1055ca60921d71315dfc03aa93baa1d8d89c96'
CHEMISTRY_SHA = '0a059a1ec5528ff140ad5634a277f6fb96022acbfec4ab9aa8f3bb450698627b'
LICENSE_SHA = '522f52ca66991e1411a22d3baf46df8f34743428e0933337e91c8d275b9e1af8'
DATA_SCRIPTS = ['audit_inventory.py', 'review_candidates.py', 'build_cohort.py',
                'build_source_bundle.py', 'fetch_chembl.py', 'fetch_documents.py',
                'fetch_compound_records.py', 'validate_sources.py', 'cohort_audit.py',
                'test_audit.py', 'test_cohort.py']
MANIFESTS = ['chembl_downloads.json', 'document_downloads.json',
             'compound_record_downloads.json', 'epizyme_patent_sources.json',
             'later_patent_sources.json', 'source_metadata.json']
REPORTS = ['activity_ledger.csv', 'assay_inventory.csv', 'inventory.json',
           'review_candidates.json', 'curation/cohort_gate.json',
           'curation/epizyme_summary.json', 'curation/later_summary.json',
           'curation/EPIZYME_SOURCE_REVIEW.md', 'curation/LATER_SOURCE_REVIEW.md',
           'curation/FEASIBILITY_UPDATE_ZH.md']
EXACT_TARGETS = ['reports/activity_ledger.csv', 'reports/assay_inventory.csv',
                 'reports/review_candidates.json', 'reports/curation/cohort_gate.json',
                 'data/processed/cohort.csv', 'data/processed/exclusions.csv',
                 'data/processed/nearest_neighbors.csv',
                 'data/curation/later/patent_crosswalk.csv',
                 'data/curation/later/binary_parent_candidates.csv',
                 'data/curation/later/binary_candidate_summary.json'] + [
    'data/curation/epizyme/results/'+name for name in [
        'epizyme_activity_curation.csv', 'epizyme_conservative_rows.csv',
        'epizyme_conservative_parents.csv', 'epizyme_curation_summary.json',
        'epizyme_protocols.json']]
PRIVATE_MARKERS = [b'/'+name+b'/' for name in (b'mnt',b'home')] + [b'/'+b'tmp/prmt5']


def require(condition, message):
    if not condition:
        raise ValueError(message)


def digest(data):
    return hashlib.sha256(data).hexdigest()


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024), b''): h.update(block)
    return h.hexdigest()


def json_bytes(value):
    return (json.dumps(value, indent=2, sort_keys=True, ensure_ascii=False)+'\n').encode()


def write_json(path, value):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_bytes(json_bytes(value))


def selected_files(root, chemistry, license_path):
    files={}
    def add(path, name=None):
        require(path.is_file() and not path.is_symlink(), 'Missing or linked input: '+str(path))
        relative=name or path.relative_to(root).as_posix()
        require(relative not in files, 'Duplicate member: '+relative)
        files[relative]=path
    for source_manifest in MANIFESTS[:3]:
        for entry in json.loads((root/'manifests'/source_manifest).read_text()):
            add(root/entry['path'])
    for directory in ['data/curation/epizyme/raw',
                      'data/curation/epizyme/results', 'data/processed']:
        for path in sorted((root/directory).rglob('*')):
            if path.is_file() and '__pycache__' not in path.parts and path.suffix not in {'.py','.pyc'}:
                add(path)
    for path in sorted((root/'data/curation/later').iterdir()):
        if path.is_file() and path.suffix in {'.html','.png','.json','.csv'} and path.name != 'artifact_hashes.json':
            add(path)
    for name in DATA_SCRIPTS: add(root/'scripts'/name)
    for directory in ['scripts/curation_epizyme','scripts/curation_later']:
        for path in sorted((root/directory).glob('*.py')): add(path)
    for name in MANIFESTS: add(root/'manifests'/name)
    for name in REPORTS: add(root/'reports'/name)
    add(root/'docs/PROTOCOL.md');add(root/'requirements/baseline.txt')
    add(chemistry, 'data/upstream/drugperturb/chemistry.py')
    add(license_path, 'data/upstream/LICENSE.drugperturb')
    require(sha(chemistry)==CHEMISTRY_SHA, 'Pinned chemistry hash mismatch')
    require(sha(license_path)==LICENSE_SHA, 'Pinned MIT notice hash mismatch')
    require(sha(root/'docs/PROTOCOL.md')==PROTOCOL_SHA, 'Protocol hash mismatch')
    require(sha(root/'data/processed/cohort.csv')==COHORT_SHA, 'Cohort hash mismatch')
    return files


def build(args):
    root=args.root.resolve()
    output=(args.archive or root/'data/releases/prmt5-sourcebench-source-snapshot.tar.gz').resolve()
    manifest=(args.manifest or root/'manifests/source_snapshot.json').resolve()
    require(not output.exists() and not manifest.exists(), 'Refusing to overwrite archive or manifest')
    files=selected_files(root,args.chemistry.resolve(),args.chemistry_license.resolve())
    members=[]
    for name,path in sorted(files.items()):
        data=path.read_bytes()
        if path.suffix.lower() in {'.json','.csv','.txt','.md','.html','.py'}:
            require(not any(marker in data for marker in PRIVATE_MARKERS), 'Private location marker in '+name)
        members.append({'path':name,'bytes':len(data),'sha256':digest(data)})
    contents={'schema_version':1,'scope':'Offline source/data reconstruction; no model weights or private runtime receipts',
              'protocol_sha256':PROTOCOL_SHA,'cohort_sha256':COHORT_SHA,
              'chemistry_commit':'cfa256c201ef45ad54efdc78eb2cef529f211d51',
              'chemistry_sha256':CHEMISTRY_SHA,'members':members,
              'excluded':'Model artifacts, private receipts, optional OCR, auxiliary data-directory scripts, BindingDB and Mirati auxiliary sources not admitted to this cohort, mutable project status/operational notes',
              'format':'Sorted regular files, uid/gid 0, mode 0644, tar mtime 0, gzip mtime 0'}
    embedded=json_bytes(contents)
    all_members=members+[{'path':'SOURCE_BUNDLE_CONTENTS.json','bytes':len(embedded),'sha256':digest(embedded)}]
    output.parent.mkdir(parents=True,exist_ok=True)
    fd,tmp=tempfile.mkstemp(prefix='.source_snapshot_',dir=output.parent);os.close(fd)
    try:
        with open(tmp,'wb') as raw, gzip.GzipFile(filename='',mode='wb',fileobj=raw,mtime=0) as zipped, tarfile.open(fileobj=zipped,mode='w',format=tarfile.PAX_FORMAT) as tar:
            for entry in sorted(all_members,key=lambda x:x['path']):
                data=embedded if entry['path']=='SOURCE_BUNDLE_CONTENTS.json' else files[entry['path']].read_bytes()
                require(digest(data)==entry['sha256'], 'Input changed during archive creation')
                info=tarfile.TarInfo(entry['path']);info.size=len(data);info.mtime=0;info.mode=0o644
                tar.addfile(info,io.BytesIO(data))
        os.replace(tmp,output)
        result={**contents,'members':all_members,'status':'SOURCE_ARCHIVE_BUILT_HASH_VERIFIED',
                'archive_path':output.relative_to(root).as_posix() if output.is_relative_to(root) else output.name,
                'archive_bytes':output.stat().st_size,'archive_sha256':sha(output),
                'member_count':len(all_members),'uncompressed_member_bytes':sum(x['bytes'] for x in all_members),
                'publication':'Local artifact only; not uploaded','fresh_environment_install_verified':False}
        verify_archive(output,result)
        write_json(manifest,result)
        print(json.dumps({key:result[key] for key in ['status','archive_sha256','archive_bytes','member_count']}))
    finally:
        if os.path.exists(tmp):os.unlink(tmp)


def verify_archive(archive, manifest, destination=None):
    require(sha(archive)==manifest['archive_sha256'],'Archive hash mismatch')
    expected={entry['path']:entry for entry in manifest['members']}
    require(len(expected)==len(manifest['members']),'Duplicate manifest names')
    seen=set()
    with tarfile.open(archive,'r:gz') as tar:
        for member in tar:
            name=member.name;parts=PurePosixPath(name)
            require(member.isfile() and not parts.is_absolute() and '..' not in parts.parts,'Unsafe member')
            require(name in expected and name not in seen,'Unexpected or duplicate member')
            data=tar.extractfile(member).read();entry=expected[name]
            require(len(data)==entry['bytes'] and digest(data)==entry['sha256'],'Member hash mismatch: '+name)
            if destination is not None:
                path=destination/name;path.parent.mkdir(parents=True,exist_ok=True);path.write_bytes(data)
            seen.add(name)
    require(seen==set(expected),'Missing archive members')


def canonical_csv_bytes(data, unordered_id_field=None):
    rows=list(csv.DictReader(io.StringIO(data.decode())))
    for row in rows:
        if unordered_id_field:
            row[unordered_id_field]=';'.join(sorted(row[unordered_id_field].split(';')))
    return Counter(json.dumps(row,sort_keys=True) for row in rows)


def auxiliary_semantics(directory, archive):
    policies={'data/curation/later/patent_crosswalk.csv':None,
              'data/curation/later/binary_parent_candidates.csv':'activity_ids',
              'data/processed/exclusions.csv':'activity_id'}
    results=[]
    with tarfile.open(archive,'r:gz') as tar:
        for name,field in policies.items():
            expected=tar.extractfile(name).read();observed=(directory/name).read_bytes()
            equal=canonical_csv_bytes(expected,field)==canonical_csv_bytes(observed,field)
            results.append({'path':name,'equal':equal,'comparison':'Exact row multiset after sorting semicolon-delimited '+field if field else 'Exact row multiset; row order ignored'})
        name='reports/curation/cohort_gate.json'
        expected=json.loads(tar.extractfile(name).read());observed=json.loads((directory/name).read_bytes())
        for role in ['inputs','outputs']:
            require(set(expected[role])==set(observed[role]),'Gate hash target set changed')
            for relative,value in observed[role].items():
                require(sha(directory/relative)==value,'Rebuilt gate has stale hashes')
            expected.pop(role);observed.pop(role)
        results.append({'path':name,'equal':expected==observed,'comparison':'All gate fields equal except byte hashes of regenerated files; all new hash targets independently verified'})
    return results


def replay(args):
    require(args.archive and args.manifest and args.replay_dir and args.report,'Replay requires archive, manifest, fresh directory and report')
    directory=args.replay_dir.resolve();require(not directory.exists(),'Replay directory must not already exist')
    manifest=json.loads(args.manifest.read_text());directory.mkdir(parents=True)
    result={'status':'RUNNING','archive_sha256':manifest['archive_sha256'],
            'fresh_environment_install_verified':False,'network_used':False,'model_fitting':False,
            'python':sys.version.split()[0],'started_utc':datetime.now(timezone.utc).isoformat(),'commands':[]}
    try:
        verify_archive(args.archive,manifest,directory)
        expected={x['path']:x['sha256'] for x in manifest['members']}
        result['verified_extracted_members']=len(expected)
        upstream='data/upstream'
        commands=[['scripts/audit_inventory.py','--root','.','--reuse-src',upstream],
                  ['scripts/review_candidates.py','--root','.'],
                  ['scripts/curation_epizyme/reproduce.py','--project-root','.'],
                  ['scripts/curation_later/curate_later.py','--root','.','--work-dir','data/curation/later'],
                  ['scripts/curation_later/crosswalk_patents.py','--root','.','--work-dir','data/curation/later'],
                  ['scripts/curation_later/make_binary_candidates.py','--root','.','--work-dir','data/curation/later'],
                  ['scripts/build_cohort.py','--root','.','--chemistry',upstream+'/drugperturb/chemistry.py']]
        logs=directory/'replay_logs';logs.mkdir()
        for i,command in enumerate(commands,1):
            log=logs/f'{i:02d}.log'
            with log.open('w') as handle:
                process=subprocess.run([sys.executable,*command],cwd=directory,stdout=handle,stderr=subprocess.STDOUT)
            result['commands'].append({'argv':['python',*command],'exit_code':process.returncode,'log_sha256':sha(log)})
            require(process.returncode==0,'Offline rebuild failed: '+command[0])
        comparisons=[{'path':name,'expected_sha256':expected[name],'observed_sha256':sha(directory/name),
                      'exact_bytes_match':expected[name]==sha(directory/name)} for name in EXACT_TARGETS]
        result['exact_comparisons']=comparisons
        semantic=auxiliary_semantics(directory,args.archive)
        result['auxiliary_semantic_comparisons']=semantic
        semantic_names={x['path'] for x in semantic}
        require(all(x['exact_bytes_match'] for x in comparisons if x['path'] not in semantic_names),'A required exact output differs')
        require(all(x['equal'] for x in semantic),'An auxiliary output changed semantically')
        gate=json.loads((directory/'reports/curation/cohort_gate.json').read_text())
        require(gate['status']=='PASS_FOR_NARROW_SOURCE_SHIFT_BINARY_BENCHMARK','Cohort gate failed')
        require([(gate['counts'][s]['parents'],gate['counts'][s]['scaffolds']) for s in ['train','validation','test']]==[(604,285),(70,44),(246,89)],'Counts changed')
        require(all(not any(x['overlap_counts'].values()) for x in gate['overlap_and_similarity'].values()),'Split overlap')
        require(sha(directory/'data/processed/cohort.csv')==COHORT_SHA,'Frozen cohort changed')
        result.update(status='PASS_EXACT_COHORT_AND_AUXILIARY_SEMANTICS',counts=gate['counts'],
                      overlap_and_similarity=gate['overlap_and_similarity'],cohort_sha256=COHORT_SHA,
                      exact_comparison_count=len(comparisons),exact_byte_matches=sum(x['exact_bytes_match'] for x in comparisons),
                      auxiliary_semantic_comparison_count=len(semantic),
                      observed_order_variance='Filesystem/glob input order changes crosswalk row order and one semicolon-separated conflicting activity-ID pair; gate hashes follow changed auxiliary bytes. No labels, structures, exclusion reasons, counts or cohort bytes change.',
                      non_exact_targets='Timestamped audit/capacity receipts and directory-content inventories are not byte-identical targets; original snapshots remain intact in archive.')
    except Exception as error:
        result.update(status='FAIL',error_type=type(error).__name__,error=str(error))
        raise
    finally:
        result['completed_utc']=datetime.now(timezone.utc).isoformat();write_json(args.report,result)
    print(json.dumps({'status':result['status'],'exact_comparison_count':result['exact_comparison_count'],'cohort_sha256':COHORT_SHA}))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode',choices=['build','replay'])
    parser.add_argument('--root',type=Path,default=Path('.'))
    parser.add_argument('--chemistry',type=Path)
    parser.add_argument('--chemistry-license',type=Path)
    parser.add_argument('--archive',type=Path)
    parser.add_argument('--manifest',type=Path)
    parser.add_argument('--replay-dir',type=Path)
    parser.add_argument('--report',type=Path)
    args=parser.parse_args()
    if args.mode=='build':
        require(args.chemistry and args.chemistry_license,'Build requires pinned chemistry and MIT notice')
        build(args)
    else:replay(args)

if __name__=='__main__':main()
