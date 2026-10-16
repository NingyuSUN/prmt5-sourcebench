# Epizyme PRMT5 原始活性表核查

对6个来源的1,964条IC50记录连接到原始Table 2。原表报告的是活性区间，不是精确连续数值；精确IC50回归不适用。采用100 nM分类阈值后，保守纳入三个家族的607行、604个母体结构、285个骨架，母体标签为164阳性、440阴性。

连接路径是活动记录record_id → ChEMBL来源化合物名称 → 本专利内编号 → 原表生化活性列。未直接重绘核验每张化学结构图。家族49958717在不同版本中存在编号变化及标签冲突，因此全族隔离，不能仅挑选数值一致的行来宣称身份匹配成功。

用于时间顺序的日期是直接核对来源版本的公开日期上界，不冒充首次公开日期。四份WO文本的2014-06-26公开日已核实，但逐化合物最早披露尚未完全证实。OCR仅作为本地辅助材料，不进入必需复现或标签。

工程核验检查可重现性、哈希和计数；不证明跨实验条件可比、外部药理有效或临床有效。550 nM代表值不能替代原始101–1000 nM区间。


Status: VERIFIED source-table and compound-record audit; NOT a training approval or scientific/clinical validation. Raw copies and row-level outputs are ignored under `data/curation/epizyme/`. Portable offline reproduction uses Python3 standard library only and does not use OCR/PDF libraries. Source URLs/hashes and verification time are recorded in `manifests/epizyme_patent_sources.json`. No training was run during this source review; completed modeling is reported separately in the final results.

## Canonical results

- `data/curation/epizyme/results/epizyme_activity_curation.csv`: all 1,964 ChEMBL IC50 activity rows across the six specified assays, including interval/censored observations that earlier `relation='='` filtering mishandled. Join is raw activity `record_id` -> ChEMBL compound record -> exact patent-specific compound alias -> indexed patent Table2 biochemical column. Original ChEMBL values/relations, candidate source identifiers, table-derived bins/bounds, exceptions, canonical parents/scaffolds and conservative publication bounds are preserved.
- `data/curation/epizyme/results/epizyme_curation_summary.json`: counts and exceptions; active threshold is **100nM**, not 1000nM. `A` is known <100nM; all B/C/D/E/F/G/*/** bins are known above100nM; unreported gap boundary values are not invented.
- `data/curation/epizyme/results/epizyme_protocols.json` and six `*_biochemical_protocol.txt`: full original assay methods and bin legends.
- `data/curation/epizyme/results/<patent>_table2.json`: exact parsed number plus biochemical/ICW/proliferation categories from machine-readable US source tables.
- `manifests/epizyme_patent_sources.json`: exact raw HTML/PDF/image SHA256 and final CSV SHA256.
- `data/curation/epizyme/optional_ocr/wo_ocr_candidate_rows.json`: 1,191 OCR candidate rows from original 2014 WO images. **Locator aid only; not manually verified, not used to assign earliest dates or labels.** OCR errors remain; these candidate rows are not label or date evidence.

## Verified central findings

No original source in this group reports exact continuous IC50 measurements. Table2 uses letter-coded intervals: A<100nM, B101–1000nM, C1001–10000nM, D10001–50000nM, E>50000nM, with US9777008 exceptions D10001–40000nM/E>40000nM. F>1000, G>5000, *>10000 and **>20000 occur where specified. ChEMBL stores bin boundary/midpoint representative numbers under `relation='='`; treating these as exact regression labels is not defensible.

The six patents form four Google-Patents family IDs:
- 50001259: US9777008B2, WO2014100695A1.
- 50029211: US9745291B2, WO2014100730A1.
- 50001260: US10150758B2, WO2014100716A1.
- 49958717: US10307413B2, US10391089B2 and US10980794B2, WO2014100719A2.
All four WO documents were published 2014-06-26 and contain structures plus Table2 activity images. Full row/chemical-identity equivalence to every later B2 is NOT established. In family49958717, numbering differs after276 and database aliases can erroneously connect different original-table labels. Do not infer first disclosure from family priority, grant year or same compound number alone.

`record_disclosed_by_date_bound` intentionally uses the exact indexed B2 source's public date as a conservative upper bound, **never a claim of first disclosure**. Earlier A1 same-number/all-three-category matches are separately shown as locators, without claiming structural identity. Every `earliest_exact_date_verified` is false. No date is derived from priority.

All six protocols explicitly describe full-length human PRMT5 and human MEP50 coexpressed in Sf9 insect cells; the assay target is human, not insect PRMT5. HistoneH4 residues1–15 substrate, 4nM enzyme,75nM tritiatedSAM,40nM peptide,30min preincubation at25C. Duration after adding SAM is not specified. Reported1%DMSO conflicts with stated1uL100%DMSO/51uL total; original text is preserved. No addedMTA is reported. Same written protocol does not establish identical reaction kinetics across sources.

## Conservative usable capacity, not validation

Excluding the entire unresolved family49958717 and ambiguous source aliases/missing compound278A leaves **607 source-label rows,604 parents,285 scaffolds** in the other three families; **164 positive parents and440 negative parents** at100nM. These indexed B2 documents were publicly available in2017/2018, before2019. Per-assay eligible rows/parents/scaffolds: CHEMBL5732693=305/304/170; CHEMBL5732555=272/271/97; CHEMBL5731839=30/29/18. Count parents once; row counts include duplicates. This is capacity only; training requires the endpoint-compatibility and split gates in the study protocol.

Ten source rows are quarantined for mixed-bin aliases or unmatched278A. Two source-vs-database inconsistencies are explicitly exposed: activity27074384 sourcecompound18 is C but ChEMBL550nM; activity28315127 exact US10980794 alias232 is B but ChEMBL100nM. The latter lies in unresolved family49958717 and is not in the conservative subset. At100nM, seven parents across that family have conflicting source-derived binary labels, affecting22rows; they are excluded from source-label eligibility. These may reflect alias/renumbering errors rather than real experimental discordance.

## Limits and reconstruction

- Chemical structures in original figures were not independently redrawn/validated. The mapping uses ChEMBL record_id and patent-specific source aliases; exact identifier ambiguity is retained.
- WO OCR is insufficient to certify exact earliest public dates. For the conservative pre2019 subset, explicit verified publication upper bounds suffice; retain exact-first-disclosure as unresolved.
- Do not use `eligible_source_label_only` as the final training filter: it does not exclude family49958717 globally and does not certify assay comparability or source independence. The final cohort explicitly excludes that family and applies the protocol split/endpoint gates.
- Reproduce from frozen original HTML and existing ChEMBL snapshots: `python3 scripts/curation_epizyme/reproduce.py --project-root .`. To independently verify deterministic outputs: `python3 scripts/curation_epizyme/reproduce.py --project-root . --output-dir data/reproduction/epizyme --verify-against data/curation/epizyme/results`. No network, OCR, package installation, or training is required.

## Exact conservative candidate filter

Select `family in {50001259,50029211,50001260}` AND `eligible_source_label_only == True` AND `record_disclosed_by_date_bound < 2019-01-01`. The eligibility flag requires no mapping/structure reason and no conflicting source binary labels for the canonical parent across the six assays. The family whitelist is mandatory: the flag alone admits unresolved49958717 rows. Deduplicate by `canonical_smiles` for parent counts, keeping earliest activity_id only as a deterministic representative; identical labels are checked before deduplication. Frozen output files are `data/curation/epizyme/results/epizyme_conservative_rows.csv` (607rows) and `epizyme_conservative_parents.csv` (604parents,164positive/440negative,285scaffolds).

### Why the entire family49958717 remains quarantined

Within the extended US10307413 table and shorter US10391089/US10980794 tables, common compound numbers after276 are not a trustworthy crosswalk. ChEMBL source aliases reuse some numbers across documents even where original categories differ. Seven parent structures carry conflicting100nM binary labels across22rows; independently validating every original chemical figure would be needed to determine whether these are mapping errors or real differences. Using only numerically matching rows would select on agreement and would not prove source identity. The entire family is excluded from this conservative set, including apparently agreeing records.

## Primary source references

Original Table2 and biochemical protocols: [US9777008B2](https://patents.google.com/patent/US9777008B2/en), [US9745291B2](https://patents.google.com/patent/US9745291B2/en), [US10150758B2](https://patents.google.com/patent/US10150758B2/en), [US10307413B2](https://patents.google.com/patent/US10307413B2/en), [US10391089B2](https://patents.google.com/patent/US10391089B2/en), [US10980794B2](https://patents.google.com/patent/US10980794B2/en). Four earliest family-level public texts: [WO2014100695A1](https://patents.google.com/patent/WO2014100695A1/en), [WO2014100730A1](https://patents.google.com/patent/WO2014100730A1/en), [WO2014100716A1](https://patents.google.com/patent/WO2014100716A1/en), [WO2014100719A2](https://patents.google.com/patent/WO2014100719A2/en). Their common2014-06-26 publication date is verified; compound-specific earliest date remains unresolved.

## Reproduction verification

VERIFIED: Two independent offline output directories reproduce all five principal output SHA256 values exactly;607rows/604parents/285scaffolds and164positive/440negative parents pass. Exact receipt: `data/curation/epizyme/results/reproduction_receipt.json`. The source manifest initially checked79entries without mismatch and now additionally tracks this receipt. No scientific/external or earliest-date validation is implied.
