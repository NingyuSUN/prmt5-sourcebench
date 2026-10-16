# Mirati original-table curation — 2026-10-16, Asia/Tokyo

## Verified result

The two ChEMBL grant documents are the same patent family (`72659905`) and the two original biochemical tables have already been disclosed in the corresponding US A1 publications dated **2021-03-18**. Each US version contains **712 IC50 cells**, comprising **563 numeric/censored values**, **148 em-dash missing values**, and **1 N.D. value**. All 712 cells, including raw strings, example numbers, condition columns and source occurrence order, agree in each of four comparisons: each B2 versus its A1; B2 versus B2; and A1 versus A1. This is the earliest **verified table disclosure among the inspected documents**, not proof that no earlier worldwide disclosure exists. Priority date 2019-09-12 is not used as a public disclosure date.

Sources: [US11479551B2](https://patents.google.com/patent/US11479551B2/en), [US20210078994A1](https://patents.google.com/patent/US20210078994A1/en), [US11492351B2](https://patents.google.com/patent/US11492351B2/en), [US20210079003A1](https://patents.google.com/patent/US20210079003A1/en). [WO2021050915A1](https://patents.google.com/patent/WO2021050915A1/en) confirms the family/date/protocol text; its image tables were **not** independently transcribed.

## Why assay names cannot determine the endpoint cohort

Table 9a: 197 source example rows, 394 cells; FlashPlate, 3 nM PRMT5/MEP50, biotinylated H4 amino acids 1–15 at 40 nM, SAM 1 µM, 20-minute preincubation and 20-minute reaction. Table 9b: 159 source example rows, 318 cells; HotSpot/H2A protocol, 1 nM PRMT5/MEP50, full-length H2A 5 µM, SAM 1 µM, 20-minute preincubation and 60-minute reaction. Both have MTA 2 µM and absent-MTA columns. The US Table 9b caption incorrectly says FlashPlate; the adjacent protocol establishes the H2A/filter-paper method, and the WO caption explicitly says HotSpot. Retain this source inconsistency rather than silently correcting the raw caption.

ChEMBL assays 5738640/41 and 5738688/89 combine both source tables under FlashPlate descriptions. These assay identifiers therefore do not define a homogeneous biochemical endpoint. Row-level table reassignment is essential.

## Activity mapping

Across the four assays, 3,378 activity records include 1,126 IC50 records and 2,252 non-IC50 records (`kon`, `k_off`), excluded from this mapping scope.

- **1,108 direct unique matches**: grant-specific example alias from compound_record + exact record/molecule/document identity + IC50 endpoint + MTA column + relation + numeric value. No compound is identified by value alone.
- **6 direct ambiguous records**: repeated source values for compatible aliases, including example 9-2 present in both tables; retain candidate rows and do not arbitrarily select one.
- **12 direct unmatched records**, all from US11492351B2: document-specific alias names omit some examples. A separate **INFERRED** family/same-molecule alias transfer gives 10 unique candidates and two ambiguous candidates. These are not promoted to direct matches.
- Of the 1,108 direct matches, 1,004 have positive exact values, 97 are positive `>` bounds, and seven are zero. There are **eight zero-valued IC50 database records overall** (four source cells duplicated across two grants). Original source tables literally contain the zeros; matching them does not make them valid logarithmic potency targets. Exclude nonpositive values and keep censoring separate.

Original Table 9b repeats example **4-152** with MTA+ values 978 and 1,910 nM. Both occurrences are retained with occurrence numbers; no unsupported correction of the example number is made.

Example **4-100** has 0.5 nM (MTA+) and 80 nM (MTA absent) in Table 9a, versus 1,510 nM and >10,000 nM in Table 9b. The 156 nM record of the same ChEMBL molecule in US11479551 is **example 4-200, Table 9b**, supported by its compound-record alias. The discrepancy is not resolved by averaging: the table conditions differ, and distinct source examples are merged by the database molecule identity. Original patent chemical structure images were not independently decoded; the identity evidence is the ChEMBL record/molecule crosswalk.

## Reproduction and deliverables

The source snapshots and row-level outputs are in `data/curation/mirati/` (ignored by Git). Reproduction scripts are in `scripts/curation_mirati/`. Run from the project root:

```
python3 scripts/curation_mirati/download_sources.py --source-dir data/curation/mirati
python3 scripts/curation_mirati/curate_mirati.py --source-dir data/curation/mirati
python3 scripts/curation_mirati/map_mirati.py --source-dir data/curation/mirati --raw-dir data/raw
python3 scripts/curation_mirati/family_alias_candidates.py --source-dir data/curation/mirati
python3 scripts/curation_mirati/validate_mirati.py --source-dir data/curation/mirati
```

The download step preserves existing snapshots by default; `--refresh` deliberately changes inputs and may change hashes. The tracked source manifest is `manifests/mirati_patent_sources.json`; the compact summary is `reports/curation/mirati_summary.json`.

- `mirati_source_table_rows.csv`: 2,848 normalized source cells; original value, patent/date/family, example/occurrence, table/protocol, MTA condition, relation, value, source SHA256, HTML locator and line.
- `mirati_disclosure_comparison.csv`: 2,848 pairwise comparison rows, all identical.
- `mirati_activity_mapping.csv`: all 3,378 records, including exclusion/unresolved statuses.
- `mirati_mapping_candidates.csv`: direct candidate locators, including ambiguous alternatives.
- `mirati_family_alias_candidates.csv`: indirect candidates explicitly marked INFERRED.
- `mirati_compound_aliases.csv`: database example-name crosswalk and record provenance.
- `mirati_uncovered_numeric_source_cells.csv`: 12 cells not covered by direct candidates.
- `mirati_parse_summary.json`, `mirati_mapping_summary.json`, `download_manifest.json`: counts, source dates, file hashes and limitations.
- `mirati_validation.json`: engineering checks passed; no scientific/external validation claimed.

## Limits and gate implications

Verified source-table dates and protocol assignments repair the previous metadata-level uncertainty. They do **not** establish comparable measurements across protocols, independent temporal samples, direct biochemical target specificity, external pharmacology, or clinical validity. Both grant documents and their A1s belong to one family and cannot be treated as independent source groups. Preserve unresolved aliases, repeated source example labels, censored values, zeros, and missingness. The final benchmark excludes this Mirati family because these MTA-dependent and protocol-specific endpoints were not pooled with its admitted catalytic source classes.
