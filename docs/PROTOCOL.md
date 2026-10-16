# PRMT5 source-shift classification protocol

Protocol: PRMT5-SSB-1, fixed on 2026-10-16 before model fitting. The machine-readable data freeze is `reports/curation/cohort_gate.json`; exact cohort bytes are bound by SHA256 there. The complete original audit remains available and is superseded only for the specific binary task below.

## Question and scope

Can structural representations learned from earlier publicly disclosed PRMT5 catalytic activity classes rank strong inhibitors in later, chemically distinct patent series? The binary target is source-reported IC50 at most 100 nM. Original bins are intervals, not precise measurements. This is retrospective transfer of source-defined labels across source, protocol, chemistry and time together. It is not a uniform potency-regression task, a prospective discovery study, or an isolated test of the effect of time.

All inputs are public. Reuse DrugPerturb chemistry standardization from commit `cfa256c201ef45ad54efdc78eb2cef529f211d51`, file SHA256 `0a059a1ec5528ff140ad5634a277f6fb96022acbfec4ab9aa8f3bb450698627b`. Preserve original activity records, structures, values, intervals, example identifiers, assay/family identifiers, date evidence and exclusions.

## Cohort and disclosure policy

| Role | Source families | Unique parents | Positive / negative | Murcko groups | Verified public date used |
|---|---|---:|---:|---:|---|
| Training | 50001259; 50029211; 50001260 | 604 | 164 / 440 | 285 | Directly indexed source publications: 2017-08-29, 2017-10-03, 2018-12-11 |
| Validation | 67846363, Pharmablock | 70 | 39 / 31 | 44 | Verified activity table: 2019-09-12 |
| Final test | 75338422, SK | 246 | 69 / 177 | 89 | Verified activity table: 2021-04-08 |

Admission cutoff for training source documents: `2019-09-01`. Three training families have WO publications on 2014-06-26, but row-by-row identity back to every early figure was not established. Use the directly verified later B2 publication as a conservative public-by upper bound; do not call all training labels first disclosed in 2014. Later table dates are verified disclosed dates, not an exhaustive worldwide first disclosure search. Patent priority dates are not disclosure dates.

Training: 607 source rows collapse to 604 isomeric canonical parents with consistent binary labels. Quarantine family 49958717 because the original chemical-example identity/renumbering is unresolved. Quarantine 10 ambiguous earlier records. Original-source bins override database surrogate numbers without destroying either. Later sources: exclude three previously disclosed reference controls; quarantine four SK parents whose original bins cross the binary threshold. Do not let post-cutoff labels decide training admission. Potential-duplicate database flags do not automatically discard a source-verified table measurement; preserve flags and deduplicate parents explicitly.

Across all three split pairs require no shared isomeric parent, nonisomeric identity, Murcko scaffold or expanded patent family. Parent/scaffold identity is checked using the pinned chemistry policy. Report nearest-neighbor fingerprint similarity rather than treating scaffold inequality as chemical independence. Do not average contradictory labels, impute missing measurements, turn zero into a small positive value, or convert assay scales.

## Assays and biological interpretation

All retained endpoints measure PRMT5/MEP50-catalyzed histone H4 methylation inhibition. Their scales differ: early FlashPlate uses SAM 75 nM; Pharmablock AlphaLISA uses SAM 10 micromolar; SK ELISA uses approximately 1 micromolar SAM and surface-coated substrate. Reaction timing, enzyme concentration and detection differ; early post-SAM reaction duration is unspecified. MTA-dependent Mirati experiments remain separate, including their different FlashPlate and HotSpot tables. Biochemical inhibition does not establish cellular activity, selectivity, exposure, toxicity, disease benefit or PRMT5 inhibitor binding mode.

## Representation and model budget

Primary common representation policy: canonical **nonisomeric** SMILES, while keeping isomeric parent identity in the cohort. This matches documented MoLFormer pretraining; it cannot distinguish stereoisomers. Count all within-split representation collisions, mixed-label collisions, unknown tokens and token lengths before fitting. Require no cross-split nonisomeric identity collision. Never silently truncate; any molecule above 202 tokens triggers protocol review before modeling, with preserved original cohort. Report representation loss and an empirical collision error bound. Primary fingerprint baselines also omit chirality for a matched information comparison; a prespecified chirality-aware logistic-regression sensitivity arm measures this limitation.

All model fitting uses the 604 training labels only. Validation labels select fixed configurations and operating thresholds; never refit on validation labels. Final-test labels are available to data curators for label checking, but not used by fitting/selection scripts. Freeze all configurations and prediction code before final evaluation. This is an openly curated retrospective split, not a blinded clinical test.

1. **Reference:** constant training prevalence, `164/604`; ranking reference AP equals each evaluation source's positive prevalence. Report expected random-selection hits and tie-aware ranking, not a lucky arbitrary sort.
2. **Similarity:** mean label of the five nearest training fingerprints, resolving equal-distance boundary ties by including all tied neighbors; Morgan radius 2, 2048 bits, no chirality. No tuning.
3. **Fingerprint logistic regression:** L2, `C={0.1,1,10}`, class weights fixed to balanced training counts, max iterations 5000. Select by validation average precision; ties prefer smaller C. Same settings for the chirality-aware sensitivity arm.
4. **Random forest:** 500 trees, max features `sqrt`, balanced class weights, `min_samples_leaf={2,5}`, seed 17. Select by validation AP; ties prefer larger leaf size.
5. **Frozen pretrained representation:** `ibm-research/MoLFormer-XL-both-10pct`, exact revision `361063d0ad524ef77cf39b08469f6be770dc550f`, mean-pooled 768-dimensional vector. Use reviewed local source and safetensors; record weight SHA256. Set `deterministic_eval=True`, call `eval()`, and require repeated/differently batched features to agree within absolute/relative tolerance `1e-5`. Training-only StandardScaler plus logistic head with the same three C values and selection rule as fingerprints.
6. **Partial fine-tuning:** adapt the last two transformer blocks and the same linear head, retaining frozen earlier layers; initialize from the selected frozen classifier, including its fixed training scaler. Two encoder learning rates `{1e-5,3e-5}`, head LR `1e-4`, AdamW weight decay `0.01`, effective batch 32, at most 20 epochs, validation AP early stopping with patience 4. Three seeds `{17,29,43}`; choose LR using mean validation AP, never select the best seed. Use deterministic encoder evaluation behavior during optimization to keep random features fixed; gradients still update the last two blocks. No architecture search. Save trainable parameter names/counts, weights, optimizer and random states, losses, validation trajectories, selected epochs and all seed predictions. The primary adapted prediction is the average probability across the three seeds at the selected LR; preserve each seed as sensitivity evidence.

The pretrained model is a modern masked-language model trained on public structure corpora; no target-label training was identified in publisher documentation. Exact pretraining membership is unverified. Unlabelled structure overlap must be disclosed separately from supervised target-label contamination. Known supervised exposure to held-out PRMT5 labels is a no-go. Do not claim a model available as of the 2019 cutoff.

## Evaluation fixed before fitting

Primary: average precision (AP), with source-positive prevalence beside it. AP asks whether positives appear early in a ranked list; it does not measure calibrated hit probability. Also report a top-10% screening budget: 7 of 70 validation molecules, 25 of 246 test molecules. At a tied boundary report expected hit count under random ordering within the tie, the attainable minimum/maximum, and a deterministic ID-ordered list explicitly labeled as arbitrary within ties. Do not silently present a fractional expectation as an observed integer hit count.

Secondary: ROC-AUC, balanced accuracy, sensitivity, specificity and integer TP/FN/FP/TN. Select each model's score threshold by maximum validation balanced accuracy; ties choose the higher threshold. Scores are not calibrated prospective hit probabilities. Save thresholds before test evaluation.

Compare all methods on the same test compounds. Use 2,000 paired scaffold-cluster bootstrap resamples (seed 20261016), retaining all members and multiplicity of each sampled test scaffold. Report percentile intervals for AP and paired AP differences versus the validation-selected strongest classical baseline; report invalid single-class resamples. These intervals describe uncertainty within this one observed patent family, not future-family variation. Raw metrics and denominators remain primary; small strata have descriptive results only. Report per-seed variation separately from bootstrap uncertainty.

Biological/error analysis is descriptive after the main result is fixed: source assay differences; training-neighbor similarity bands `<0.5`, `0.5–<0.7`, `>=0.7`; within-source activity bins; stereochemical representation collisions; predeclared physicochemical descriptors (MW, logP, H-bond donors/acceptors, TPSA). Inspect all available named reference controls as provenance/assay context, not hidden extra test labels or a new tuning set. Fingerprint attribution is a hypothesis about statistical association, not a molecular binding mechanism. No docking is required to create an unsupported mechanism narrative.

## Acceptance, failures and delivery

Scientific completion means faithful evaluation and interpretation, including a negative result. Improvement is not a required outcome. Engineering acceptance requires hash-verified input/output, frozen cohort equality, split isolation, finite features/loss, deterministic evaluation, correct metric edge cases, one-command scoring from sealed predictions and successful independent code/result review. These checks do not establish experimental validation.

Stop the affected modeling arm for unresolved label/identity leakage, corrupted source/weight bytes, unknown assays, nonfinite training, unsupported tokenizer inputs, or known supervised pretraining contamination. Investigate resource failure without silently changing data, model revision or evaluation rules. Keep failed runs and explain any protocol amendment before using final-test scores.

Deliver a data/source card, source-date/family ledger, model card, all selected and alternative validation runs, final predictions/metrics with uncertainty, a biology-facing report in English and Chinese, pinned runtime/requirements and small reproducible scoring fixtures. Raw source documents and model weights remain outside Git with acquisition manifests. The analysis reuses the pinned chemistry utilities. No wet experiment is included.
