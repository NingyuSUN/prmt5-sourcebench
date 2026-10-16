# Data card — PRMT5 SourceBench

**920 parent structures with audited original-source biochemical activity classes.** This is a retrospective source-shift classification dataset; exact IC50 regression is not supported by the interval-only sources.

| Split | Parents | Positive / negative | Scaffolds | Families |
|---|---:|---:|---:|---:|
| Training | 604 | 164 / 440 | 285 | 3 |
| Validation | 70 | 39 / 31 | 44 | 1 |
| Final test | 246 | 69 / 177 | 89 | 1 |

The original ChEMBL 37 inventory contains 17,697 activities, 608 assays and 171 documents. Those are database records, not 17,697 usable independent compounds. The frozen cohort is the source-verified subset. An additional 2,973 compound records establish original example aliases. The [family/date ledger](../reports/curation/source_family_dates.csv) identifies every admitted family and its counts.

## Labels and source comparison

- Early Epizyme FlashPlate sources report A/B/C/D activity ranges; A<100 nM is positive. Original gaps/boundaries are preserved. Database numbers must not be treated as exact IC50 values.
- Pharmablock AlphaLISA reports numeric values. After excluding three old reference controls, 70 parents remain; IC50 ≤100 nM defines the binary target.
- SK ELISA reports +++=[1,100] nM, ++=(100,1000] nM, +=(1000,10000] nM. Numeric database values such as 50/550/5500 are surrogates, not measured point estimates. Four parents crossing the binary boundary are excluded; a parent with two negative bins retains both bins and remains negative.

All admitted assays concern PRMT5/MEP50-catalyzed histone H4 methylation, but SAM concentrations, enzyme/substrate conditions and readout differ. Early SAM 75 nM, Pharmablock 10 µM and SK approximately 1 µM must not be silently pooled onto a common potency scale. MTA-dependent Mirati sources were separately reviewed and are excluded from this cohort. See [original-table review](../reports/curation/LATER_SOURCE_REVIEW.md) and the [frozen protocol](PROTOCOL.md).

## Dates and exclusions

Training uses verified public-by bounds in 2017–2018, preceding the 2019 validation and 2021 test tables. Earlier family publications in 2014 were identified, but every exact structure–label pair was not reconstructed back to those earliest chemical figures. Priority dates and patent grant years do not automatically establish first disclosure. These are conservative table/disclosure bounds, not an exhaustive first-ever-publication claim.

607 training source rows collapse to 604 parents. The whole unresolved family 49958717 (1,347 source records) and 10 ambiguous earlier records are quarantined; 3 previously disclosed controls and 4 later threshold-conflict parents are excluded. The 145 included training source records carrying a potential-duplicate flag retain that flag and original evidence; explicit parent deduplication prevents counting them as independent compounds. No conflicting labels are averaged or majority-voted.

## Structural separation and remaining uncertainty

All split pairs have zero shared isomeric parent, nonisomeric identity, Murcko scaffold and patent family. Maximum training–test chiral Morgan similarity is 0.62963; 142/246 test compounds are below 0.5. Different scaffold strings alone do not prove full chemical independence. Primary model inputs omit chirality; the curated identity/source ledger retains it.

One test family cannot estimate future-family variability. Time, source, protocol and series change together. Original chemical figures are not all independently redrawn; compound-record/example correspondence is the documented provenance. No prospective, cellular, assay replication or clinical validation is claimed.

## Evidence and delivery

[Source manifests](../manifests/source_snapshot.json), [cohort gate](../reports/curation/cohort_gate.json), [offline reconstruction](../reports/curation/source_snapshot_replay.json) and [reproduction guide](REPRODUCIBILITY.md) bind exact inputs and outputs. Local source bundle:226 members, 38,080,957 bytes. The frozen cohort reproduces byte-for-byte; four auxiliary outputs have explicitly documented ordering-only differences. All raw source bytes and original values remain preserved. ChEMBL data are subject to their own attribution/share-alike terms; original patents and pretrained assets have separate provenance. The code license does not relicense external data.
