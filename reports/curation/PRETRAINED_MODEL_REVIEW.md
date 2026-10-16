# Pretrained model provenance

The completed experiment uses `ibm-research/MoLFormer-XL-both-10pct` at revision `361063d0ad524ef77cf39b08469f6be770dc550f`. [All fitted configurations](../../docs/MODEL_CARD.md); [implementation](../../docs/ENCODER_IMPLEMENTATION.md).

The [publisher card](https://huggingface.co/ibm-research/MoLFormer-XL-both-10pct/blob/361063d0ad524ef77cf39b08469f6be770dc550f/README.md) describes masked-language pretraining on 10% each of ZINC15 and PubChem, canonical nonisomeric RDKit SMILES, and a 202-token limit. This is the released 10%+10% checkpoint, not the full 1.1-billion-molecule variant. [Paper](https://doi.org/10.1038/s42256-022-00580-7); [implementation](https://github.com/IBM/molformer).

The architecture has 12 layers, hidden width 768, 12 attention heads and random-feature attention. The experiment loads pinned safetensors weights and reviewed local custom-code files, preserves twelve pretrained random-feature buffers, and uses deterministic evaluation behavior during feature extraction and partial adaptation. All 920 structures passed input checks without truncation or unknown symbols. Runtime identities and assets are recorded in [preflight](../models/encoder_preflight.json) and [asset manifest](../../manifests/pretrained_model.json).

Source was inspected for direct network, subprocess and file-writing behavior before execution. This bounded inspection is not a comprehensive security audit. Exact source and weight hashes identify the inspected bytes; model assets are not redistributed.

No supervised PRMT5-label training was identified in inspected publisher documentation. Exact structure-corpus membership and comprehensive absence of label contamination remain unverified. Unlabelled structure overlap is distinct from supervised label leakage. A modern model supports retrospective comparison, not a claim that an equivalent predictor existed before the source publications.

The publisher declares Apache-2.0. Retain applicable license/notices when acquiring or redistributing assets. Actual parameter budgets, validation selection, all six adaptation trials and the negative held-out result are in the model card; no pretraining benchmark performance is transferred to this task by assumption.
