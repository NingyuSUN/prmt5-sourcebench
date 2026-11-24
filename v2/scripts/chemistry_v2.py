"""Deterministic chemical structure standardization and scaffold identity for PRMT5 v2.

Preserves stereochemistry in `canonical_smiles`, computes stereochemistry-stripped
`nonisomeric_smiles` to detect enantiomer/diastereomer collisions across splits,
and derives canonical Bemis-Murcko scaffolds (`scaffold_group`).
"""
from __future__ import annotations

import hashlib
from rdkit import Chem, RDLogger
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.Scaffolds import MurckoScaffold

RDLogger.DisableLog("rdApp.warning")
RDLogger.DisableLog("rdApp.error")

_LARGEST_FRAGMENT = rdMolStandardize.LargestFragmentChooser(preferOrganic=True)
_UNCHARGER = rdMolStandardize.Uncharger()
_TAUTOMER = rdMolStandardize.TautomerEnumerator()


def standardize_structure(smiles: str) -> dict[str, str]:
    """Standardize raw SMILES to parent isomeric SMILES, nonisomeric SMILES, and Murcko scaffold."""
    if not smiles or not isinstance(smiles, str):
        raise ValueError("Empty SMILES string")
    mol = Chem.MolFromSmiles(smiles.strip())
    if mol is None or mol.GetNumHeavyAtoms() == 0:
        raise ValueError(f"Invalid SMILES: {smiles!r}")

    clean = rdMolStandardize.Cleanup(mol)
    parent = _LARGEST_FRAGMENT.choose(clean)
    uncharged = _UNCHARGER.uncharge(parent)
    taut = _TAUTOMER.Canonicalize(uncharged)
    Chem.SanitizeMol(taut)

    iso_smiles = Chem.MolToSmiles(taut, canonical=True, isomericSmiles=True)
    # Re-parse canonical SMILES to guarantee idempotence
    iso_mol = Chem.MolFromSmiles(iso_smiles)
    if iso_mol is None:
        raise ValueError(f"Canonical SMILES failed re-parse: {iso_smiles!r}")
    iso_smiles = Chem.MolToSmiles(iso_mol, canonical=True, isomericSmiles=True)
    noniso_smiles = Chem.MolToSmiles(iso_mol, canonical=True, isomericSmiles=False)

    # Bemis-Murcko scaffold (nonisomeric ring/linker framework, matching v1 Murcko policy)
    achiral_mol = Chem.MolFromSmiles(noniso_smiles)
    scaff_mol = MurckoScaffold.GetScaffoldForMol(achiral_mol)
    scaff_smi = Chem.MolToSmiles(scaff_mol, canonical=True, isomericSmiles=False) if scaff_mol else ""
    scaffold_group = f"murcko:{scaff_smi}" if scaff_smi else f"acyclic:{noniso_smiles}"

    compound_id = hashlib.sha256(iso_smiles.encode("utf-8")).hexdigest()[:16]
    return {
        "compound_id": compound_id,
        "canonical_smiles": iso_smiles,
        "nonisomeric_smiles": noniso_smiles,
        "scaffold_group": scaffold_group,
        "has_defined_chirality": str("@" in iso_smiles),
    }
