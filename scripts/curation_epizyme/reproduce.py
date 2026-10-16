"""Offline Epizyme source curation. Python standard library only; no OCR needed."""
import argparse
import collections
import hashlib
import json
import subprocess
import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--project-root', type=Path, default=SCRIPT_DIR.parents[1])
parser.add_argument('--raw-dir', type=Path)
parser.add_argument('--output-dir', type=Path)
parser.add_argument('--verify-against', type=Path, help='Compare deterministic CSV/summary hashes to an existing output directory.')
args = parser.parse_args()
project = args.project_root.resolve()
raw = (args.raw_dir or project/'data/curation/epizyme/raw').resolve()
out = (args.output_dir or project/'data/curation/epizyme/results').resolve()
out.mkdir(parents=True, exist_ok=True)
subprocess.run([sys.executable, str(SCRIPT_DIR/'extract_patents.py'), '--raw-dir', str(raw), '--output-dir', str(out)], check=True)
subprocess.run([sys.executable, str(SCRIPT_DIR/'curate.py'), '--project-root', str(project), '--patent-dir', str(out), '--output-dir', str(out)], check=True)
summary = json.loads((out/'epizyme_curation_summary.json').read_text())
protocols = []
for value in summary.values():
    if 'patent' not in value:
        continue
    patent = value['patent']
    doc = json.loads((out/(patent+'.json')).read_text())
    text = (out/(patent+'.txt')).read_text()
    start = text.find('Biochemical Assay')
    end = text.find('Z-138 Methylation Assay', start)
    assert start >= 0 and end > start, patent
    method = text[start:end]
    legend_start = text.find('For Table 2')
    legend_end = text.lower().find('other embodiments', legend_start)
    assert legend_start >= 0 and legend_end > legend_start, patent
    legend = text[legend_start:legend_end]
    (out/(patent+'_biochemical_protocol.txt')).write_text(method+'\n\n'+legend)
    protocols.append({
        'patent': patent, 'source_url': doc['source'], 'source_sha256': doc['sha256'],
        'family_id': doc['family'], 'direct_publication_date': doc['meta']['DC.date'],
        'organism': 'human', 'construct': 'full-length human FLAG-PRMT5 (NM_006109.3 variant1) /6His-MEP50 (NM_024102)',
        'expression': 'coexpression in Sf9 insect cells',
        'substrate': 'histone H4 residues1-15, C-terminal biotin tag; Ac-SGRGKGGKGLGKGGA[K-Biot]-amide',
        'enzyme_nM': 4, 'SAM_nM': 75, 'peptide_nM': 40,
        'preincubation_minutes': 30, 'preincubation_temperature_C': 25,
        'reaction_duration_minutes': None, 'reaction_duration_status': 'not specified in inspected original biochemical assay section',
        'buffer': '20mM Bicine pH7.6,1mM TCEP,0.005%BSG,0.002%Tween20',
        'readout': 'tritium methyl transfer to H4 peptide, streptavidin Flashplate/TopCount',
        'stop': 'coldSAM final600uM', 'MTA': 'not reported as added',
        'DMSO_source_reported_percent': 1,
        'DMSO_volume_inconsistency': '1uL100%DMSO in51uL implies approximately1.96%; preserve source1% without correction',
        'interval_legend': legend, 'protocol_text_file': patent+'_biochemical_protocol.txt'
    })
(out/'epizyme_protocols.json').write_text(json.dumps(protocols, ensure_ascii=False, indent=2))
products = ['epizyme_activity_curation.csv', 'epizyme_conservative_rows.csv', 'epizyme_conservative_parents.csv', 'epizyme_curation_summary.json', 'epizyme_protocols.json']
hashes = {name: hashlib.sha256((out/name).read_bytes()).hexdigest() for name in products}
expected = {'rows': 607, 'parents': 604, 'scaffolds': 285}
actual = summary['conservative_candidates']['counts']
checks = {key: actual[key] == value for key,value in expected.items()}
checks['parent_labels_164_positive_440_negative'] = summary['conservative_candidates']['parent_labels'] == {'0': 440, '1': 164}
if args.verify_against:
    for name,digest in hashes.items():
        checks['reproduction_sha256_'+name] = hashlib.sha256((args.verify_against/name).read_bytes()).hexdigest() == digest
receipt = {'status': 'PASS_ENGINEERING_ONLY' if all(checks.values()) else 'FAIL', 'checks': checks, 'output_sha256': hashes, 'scientific_validation': False, 'earliest_date_validation': False}
(out/'reproduction_receipt.json').write_text(json.dumps(receipt, indent=2))
print(json.dumps(receipt, indent=2))
if not all(checks.values()):
    raise SystemExit(1)
