# Later PRMT5 source curation

Run from the canonical project root with Python 3.12 and RDKit 2026.09.1 available. The first script needs RDKit; the other scripts use the standard library. No model training occurs.

Order matters:

```bash
python scripts/curation_later/curate_later.py --root . --work-dir data/curation/later
python scripts/curation_later/crosswalk_patents.py --root . --work-dir data/curation/later
python scripts/curation_later/make_binary_candidates.py --root . --work-dir data/curation/later
```

`--root` points to the project containing the preserved ChEMBL raw snapshots, activity ledger and review candidates. `--work-dir` contains original patent HTML, original table images, manual WO transcription, `reproduction_inputs.json`, and generated outputs. It must exist before execution. All arguments are explicit; no machine-specific path is embedded.

Every script validates preserved source inputs against `reproduction_inputs.json`; changed/missing inputs fail. The crosswalk script checks 329 activities (256 SK +73 Pharmablock) against real compound-record aliases and primary table labels, and verifies all75 original WO cells against the B2 table. The final script checks323 assay-parent rows and316 candidate includes (SK246=69positive+177negative; Pharmablock70=39positive+31negative), with4 threshold-conflict parents and3 old reference controls excluded. The329 activities must not be described as329 independent compounds.

Source dates mean verified activity-table publication dates, not an exhaustive first structure-label disclosure search. SK database numeric surrogates are preserved for traceability, never promoted to exact IC50. Original table bins and conflicts remain available. Original source pages and curation snapshots are not distributed. The public release includes selected derived cohort, exclusion and prediction tables; see data/README.md for their scope.
