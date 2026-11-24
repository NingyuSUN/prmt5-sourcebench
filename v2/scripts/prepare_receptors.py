#!/usr/bin/env python3
"""Prepare and superimpose PRMT5 co-crystal structures onto the 4X61 catalytic domain frame.

Strictly separates pre-2023 experimental crystal receptors (`receptors_pre2023/`, permitted
for Stage 2 featurization and docking) from post-2023 clinical co-crystals
(`diagnostic_post2023/`, strictly quarantined from featurization/training and used only
for post-hoc biophysical visualization).
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np
from Bio.PDB import MMCIFParser, Superimposer
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from rdkit import Chem

PRE2023_SPECS = [
    {"pdb_id": "4X61", "state": "State1_SAM_Substrate", "release_date": "2015-04-22", "resolution_A": 2.85, "cofactor": "SAM", "ligand_id": "3XV", "ligand_name": "EPZ015666", "vina_receptor": True},
    {"pdb_id": "5EML", "state": "State1_SAM_Substrate", "release_date": "2016-02-24", "resolution_A": 2.39, "cofactor": "SAM", "ligand_id": "5QK", "ligand_name": "Epizyme_Cmpd10", "vina_receptor": False},
    {"pdb_id": "7L1G", "state": "State1_SAM_Binary", "release_date": "2021-04-14", "resolution_A": 2.47, "cofactor": "SAM", "ligand_id": None, "ligand_name": None, "vina_receptor": False},
    {"pdb_id": "7S1S", "state": "State2_MTA_Cooperative", "release_date": "2022-01-26", "resolution_A": 2.62, "cofactor": "MTA", "ligand_id": "85K", "ligand_name": "MRTX1719", "vina_receptor": True},
    {"pdb_id": "8CSG", "state": "State2_MTA_Cooperative", "release_date": "2022-10-05", "resolution_A": 2.48, "cofactor": "MTA", "ligand_id": "PWL", "ligand_name": "Mirati_Frag1", "vina_receptor": False},
    {"pdb_id": "6CKC", "state": "State3_SAM_Nucleoside", "release_date": "2018-05-16", "resolution_A": 2.80, "cofactor": None, "ligand_id": "F5J", "ligand_name": "LLY-283", "vina_receptor": True},
    {"pdb_id": "7MX7", "state": "State3_SAM_Nucleoside", "release_date": "2021-11-17", "resolution_A": 2.49, "cofactor": None, "ligand_id": "ZR1", "ligand_name": "PF-06939999", "vina_receptor": False},
    {"pdb_id": "6RLL", "state": "State3_SAM_Nucleoside", "release_date": "2020-06-03", "resolution_A": 2.22, "cofactor": None, "ligand_id": "K8H", "ligand_name": "JNJ-44064146", "vina_receptor": False},
]

DIAGNOSTIC_POST2023_SPECS = [
    {"pdb_id": "8VEO", "state": "State2_MTA_Binary", "release_date": "2024-04-24", "resolution_A": 2.03, "cofactor": "MTA", "ligand_id": None, "ligand_name": None},
    {"pdb_id": "8VEY", "state": "State2_MTA_Cooperative", "release_date": "2024-04-24", "resolution_A": 2.44, "cofactor": "MTA", "ligand_id": "A1AAV", "ligand_name": "TNG908"},
    {"pdb_id": "9C10", "state": "State2_MTA_Cooperative", "release_date": "2024-10-02", "resolution_A": 2.85, "cofactor": "MTA", "ligand_id": "A1ATH", "ligand_name": "AMG_193"},
    {"pdb_id": "9N3R", "state": "State2_MTA_Cooperative", "release_date": "2025-03-12", "resolution_A": 2.47, "cofactor": "MTA", "ligand_id": "A1BV0", "ligand_name": "TNG462"},
    {"pdb_id": "9T43", "state": "State2_MTA_Cooperative", "release_date": "2026-06-10", "resolution_A": 2.18, "cofactor": "MTA", "ligand_id": "A1JTI", "ligand_name": "AZD3470"},
]

AROMATIC_RES_ATOMS = {
    "PHE": {"CG", "CD1", "CD2", "CE1", "CE2", "CZ"},
    "TYR": {"CG", "CD1", "CD2", "CE1", "CE2", "CZ"},
    "TRP": {"CG", "CD1", "CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"},
    "HIS": {"CG", "CD2", "CE1"},
    "SAM": {"C4", "C5", "C6", "C2", "C8"},
    "MTA": {"C4", "C5", "C6", "C2", "C8"},
}

ACCEPTOR_N = {("HIS", "ND1"), ("SAM", "N1"), ("SAM", "N3"), ("SAM", "N7"), ("MTA", "N1"), ("MTA", "N3"), ("MTA", "N7")}


def sha256_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _as_list(val: Any) -> list[str]:
    if val is None:
        return []
    return [val] if isinstance(val, str) else list(val)


def build_rdkit_mol_from_mmcif_residue(res, mmcif_dict: dict, comp_id: str) -> Chem.Mol:
    """Build an exact bond-order sanitized RDKit 3D Mol from Bio.PDB residue + _chem_comp_bond."""
    atom_map = {}
    rw = Chem.RWMol()
    coords = []

    # Map formal charges from _chem_comp_atom if present
    ca_cids = _as_list(mmcif_dict.get("_chem_comp_atom.comp_id"))
    ca_aids = _as_list(mmcif_dict.get("_chem_comp_atom.atom_id"))
    ca_chgs = _as_list(mmcif_dict.get("_chem_comp_atom.charge"))
    charge_lookup = {}
    for c, a, q in zip(ca_cids, ca_aids, ca_chgs):
        if c == comp_id:
            try:
                charge_lookup[a] = int(float(q))
            except ValueError:
                charge_lookup[a] = 0

    for atom in res:
        el = (atom.element or "").strip().capitalize()
        if el in ("H", "D", ""):
            continue
        if len(el) == 2:
            el = el[0].upper() + el[1].lower()
        rd_atom = Chem.Atom(el)
        aname = atom.name.strip()
        if aname in charge_lookup and charge_lookup[aname] != 0:
            rd_atom.SetFormalCharge(charge_lookup[aname])
        rd_atom.SetProp("pdb_atom_name", aname)
        idx = rw.AddAtom(rd_atom)
        atom_map[aname] = idx
        coords.append(atom.coord.astype(float))

    b_cids = _as_list(mmcif_dict.get("_chem_comp_bond.comp_id"))
    b_a1 = _as_list(mmcif_dict.get("_chem_comp_bond.atom_id_1"))
    b_a2 = _as_list(mmcif_dict.get("_chem_comp_bond.atom_id_2"))
    b_ord = _as_list(mmcif_dict.get("_chem_comp_bond.value_order"))
    order_map = {
        "SING": Chem.BondType.SINGLE,
        "DOUB": Chem.BondType.DOUBLE,
        "TRIP": Chem.BondType.TRIPLE,
        "AROM": Chem.BondType.AROMATIC,
    }
    for cid, a1, a2, bo in zip(b_cids, b_a1, b_a2, b_ord):
        if cid == comp_id and a1 in atom_map and a2 in atom_map:
            rw.AddBond(atom_map[a1], atom_map[a2], order_map.get(bo.upper(), Chem.BondType.SINGLE))

    mol = rw.GetMol()
    conf = Chem.Conformer(mol.GetNumAtoms())
    for i, pt in enumerate(coords):
        conf.SetAtomPosition(i, pt.tolist())
    mol.AddConformer(conf, assignId=True)
    try:
        Chem.SanitizeMol(mol)
    except Exception:
        Chem.SanitizeMol(mol, sanitizeOps=Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_KEKULIZE)
    return mol


def write_vina_receptor_pdbqt(model, pocket_center: np.ndarray, include_cofactor: str | None, out_path: Path, cutoff_A: float = 20.0) -> dict[str, Any]:
    """Write an AutoDock4-typed PDBQT receptor pocket including explicit polar H (HD) donors."""
    lines = []
    pocket_atoms_json = []
    serial = 1
    prev_c_coord = {}

    for ch in model:
        prev_c = None
        for r in ch:
            rn = r.resname.strip()
            is_prot = r.id[0] == " "
            is_cof = (include_cofactor is not None and rn == include_cofactor and ch.id == "A")
            if not (is_prot or is_cof):
                continue
            heavy_atoms = [a for a in r if (a.element or "").strip().upper() not in ("H", "D", "")]
            if not heavy_atoms:
                continue
            rcoords = np.array([a.coord for a in heavy_atoms], dtype=float)
            min_dist = float(np.linalg.norm(rcoords - pocket_center[None, :], axis=1).min())
            if min_dist > cutoff_A:
                if is_prot and "C" in r:
                    prev_c = r["C"].coord.astype(float)
                continue

            atom_by_name = {a.name.strip(): a.coord.astype(float) for a in heavy_atoms}
            for a in heavy_atoms:
                anm = a.name.strip()
                el = (a.element or "C").strip().upper()
                x, y, z = [float(v) for v in a.coord]
                if el == "C":
                    ad = "A" if anm in AROMATIC_RES_ATOMS.get(rn, set()) else "C"
                elif el == "N":
                    ad = "NA" if (rn, anm) in ACCEPTOR_N else "N"
                elif el == "O":
                    ad = "OA"
                elif el == "S":
                    ad = "SA"
                elif el == "P":
                    ad = "P"
                elif el in ("F", "FE", "ZN", "MG", "CL", "BR", "I"):
                    ad = el.capitalize()
                else:
                    ad = "C"

                resseq = int(r.id[1])
                lines.append(
                    f"ATOM  {serial:5d} {anm:<4s} {rn:>3s} {ch.id:1s}{resseq:4d}    "
                    f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00     0.000 {ad:<2s}"
                )
                pocket_atoms_json.append({
                    "chain": ch.id, "resname": rn, "resseq": resseq, "atom": anm,
                    "element": el, "ad4_type": ad, "xyz": [round(x, 4), round(y, 4), round(z, 4)],
                })
                serial += 1

                # Place geometric polar H (HD) on backbone amide N (non-PRO) so Vina scores H-bonds properly
                if is_prot and anm == "N" and rn != "PRO" and prev_c is not None and "CA" in atom_by_name:
                    v1 = prev_c - atom_by_name["N"]
                    v2 = atom_by_name["CA"] - atom_by_name["N"]
                    n1 = np.linalg.norm(v1)
                    n2 = np.linalg.norm(v2)
                    if n1 > 0.5 and n2 > 0.5:
                        bis = v1 / n1 + v2 / n2
                        bn = np.linalg.norm(bis)
                        if bn > 1e-3:
                            h_pt = atom_by_name["N"] - 1.01 * (bis / bn)
                            lines.append(
                                f"ATOM  {serial:5d} {'H':<4s} {rn:>3s} {ch.id:1s}{resseq:4d}    "
                                f"{h_pt[0]:8.3f}{h_pt[1]:8.3f}{h_pt[2]:8.3f}  1.00  0.00     0.000 {'HD':<2s}"
                            )
                            serial += 1

            # Place geometric polar H (HD) on key sidechain/cofactor hydroxyl & amine donors
            donor_pairs = [
                ("SER", "OG", "CB"), ("THR", "OG1", "CB"), ("TYR", "OH", "CZ"),
                ("LYS", "NZ", "CE"), ("CYS", "SG", "CB"),
                ("SAM", "N6", "C6"), ("SAM", "O2'", "C2'"), ("SAM", "O3'", "C3'"),
                ("MTA", "N6", "C6"), ("MTA", "O2'", "C2'"), ("MTA", "O3'", "C3'"),
            ]
            for drn, don, base in donor_pairs:
                if rn == drn and don in atom_by_name and base in atom_by_name:
                    vec = atom_by_name[don] - atom_by_name[base]
                    vn = np.linalg.norm(vec)
                    if vn > 0.5:
                        h_pt = atom_by_name[don] + 0.98 * (vec / vn)
                        resseq = int(r.id[1])
                        lines.append(
                            f"ATOM  {serial:5d} {'HD1':<4s} {rn:>3s} {ch.id:1s}{resseq:4d}    "
                            f"{h_pt[0]:8.3f}{h_pt[1]:8.3f}{h_pt[2]:8.3f}  1.00  0.00     0.000 {'HD':<2s}"
                        )
                        serial += 1

            if is_prot and "C" in atom_by_name:
                prev_c = atom_by_name["C"]

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return {"pdbqt_lines": len(lines), "heavy_atoms": len(pocket_atoms_json), "atoms": pocket_atoms_json}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    root = args.root.resolve()

    raw_dir = root / "data/structures/raw_mmcif"
    pre_dir = root / "data/structures/aligned/receptors_pre2023"
    diag_dir = root / "data/structures/aligned/diagnostic_post2023"
    pre_dir.mkdir(parents=True, exist_ok=True)
    diag_dir.mkdir(parents=True, exist_ok=True)

    cif_parser = MMCIFParser(QUIET=True)

    def load_entry(pdb_id: str):
        cif_path = raw_dir / f"pdb_0000{pdb_id.lower()}.cif.gz"
        if not cif_path.exists():
            raise FileNotFoundError(f"Missing PDB mmCIF: {cif_path}")
        with gzip.open(cif_path, "rt") as fh:
            mmcif_dict = MMCIF2Dict(fh)
        with gzip.open(cif_path, "rt") as fh:
            struct = cif_parser.get_structure(pdb_id, fh)
        return cif_path, mmcif_dict, struct[0]

    # Load 4X61 reference catalytic domain (residues 300-630 of chain A)
    ref_cif, ref_dict, ref_model = load_entry("4X61")
    ref_ca = {
        r.id[1]: r["CA"]
        for r in ref_model["A"]
        if r.id[0] == " " and "CA" in r and 300 <= r.id[1] <= 630
    }

    manifest_entries = []
    pocket_centers = {}

    all_specs = [(spec, False) for spec in PRE2023_SPECS] + [(spec, True) for spec in DIAGNOSTIC_POST2023_SPECS]

    for spec, is_diag in all_specs:
        pdb_id = spec["pdb_id"]
        out_base = diag_dir if is_diag else pre_dir
        cif_path, mmcif_dict, model = load_entry(pdb_id)

        common_res = sorted(set(ref_ca) & {r.id[1] for r in model["A"] if r.id[0] == " " and "CA" in r})
        sup = Superimposer()
        sup.set_atoms([ref_ca[i] for i in common_res], [model["A"][(" ", i, " ")]["CA"] for i in common_res])
        sup.apply(model.get_atoms())
        rmsd = float(sup.rms)

        rec_info: dict[str, Any] = {
            **spec,
            "quarantined_post2023_diagnostic_only": is_diag,
            "raw_mmcif_path": str(cif_path.relative_to(root)),
            "raw_mmcif_sha256": sha256_file(cif_path),
            "aligned_ca_count": len(common_res),
            "ca_rmsd_to_4X61_A": round(rmsd, 4),
        }

        # Extract cofactor SDF & JSON coordinates if present
        cof_id = spec["cofactor"]
        if cof_id:
            cof_res = next(r for r in model["A"] if r.resname.strip() == cof_id)
            cof_mol = build_rdkit_mol_from_mmcif_residue(cof_res, mmcif_dict, cof_id)
            cof_sdf = out_base / f"{pdb_id}_cofactor_{cof_id}.sdf"
            w = Chem.SDWriter(str(cof_sdf))
            cof_mol.SetProp("_Name", f"{pdb_id}_{cof_id}")
            w.write(cof_mol)
            w.close()
            cof_atoms = {
                a.name.strip(): [round(float(v), 4) for v in a.coord]
                for a in cof_res if (a.element or "").strip().upper() not in ("H", "D", "")
            }
            cof_json = out_base / f"{pdb_id}_cofactor_{cof_id}.json"
            cof_json.write_text(json.dumps({"pdb_id": pdb_id, "cofactor": cof_id, "atoms": cof_atoms}, indent=2) + "\n")
            rec_info["cofactor_sdf"] = str(cof_sdf.relative_to(root))
            rec_info["cofactor_sdf_sha256"] = sha256_file(cof_sdf)
            rec_info["cofactor_json"] = str(cof_json.relative_to(root))
            rec_info["cofactor_json_sha256"] = sha256_file(cof_json)

        # Extract crystallographic inhibitor ligand SDF if present
        lig_id = spec["ligand_id"]
        lig_center = None
        if lig_id:
            lig_res = next(r for r in model["A"] if r.resname.strip() == lig_id)
            lig_mol = build_rdkit_mol_from_mmcif_residue(lig_res, mmcif_dict, lig_id)
            lig_sdf = out_base / f"{pdb_id}_ligand_{lig_id}.sdf"
            w = Chem.SDWriter(str(lig_sdf))
            lig_mol.SetProp("_Name", f"{pdb_id}_{lig_id}_{spec['ligand_name']}")
            w.write(lig_mol)
            w.close()
            lig_coords = np.array([a.coord for a in lig_res if (a.element or "").strip().upper() not in ("H", "D", "")], dtype=float)
            lig_center = lig_coords.mean(axis=0)
            pocket_centers[pdb_id] = [round(float(v), 4) for v in lig_center]
            rec_info["ligand_sdf"] = str(lig_sdf.relative_to(root))
            rec_info["ligand_sdf_sha256"] = sha256_file(lig_sdf)
            rec_info["ligand_heavy_atoms"] = int(len(lig_coords))
            rec_info["pocket_center_xyz"] = pocket_centers[pdb_id]

        # Build Vina PDBQT + heavy-atom pocket JSON for pre-2023 Vina receptors
        if not is_diag and spec.get("vina_receptor") and lig_center is not None:
            pdbqt_path = out_base / f"{pdb_id}_receptor.pdbqt"
            pkt_data = write_vina_receptor_pdbqt(model, lig_center, cof_id, pdbqt_path, cutoff_A=20.0)
            pkt_json = out_base / f"{pdb_id}_pocket_atoms.json"
            pkt_json.write_text(json.dumps({
                "pdb_id": pdb_id,
                "state": spec["state"],
                "cofactor": cof_id,
                "pocket_center_xyz": pocket_centers[pdb_id],
                "heavy_atoms": pkt_data["heavy_atoms"],
                "atoms": pkt_data["atoms"],
            }, indent=2) + "\n")
            rec_info["receptor_pdbqt"] = str(pdbqt_path.relative_to(root))
            rec_info["receptor_pdbqt_sha256"] = sha256_file(pdbqt_path)
            rec_info["pocket_atoms_json"] = str(pkt_json.relative_to(root))
            rec_info["pocket_atoms_json_sha256"] = sha256_file(pkt_json)

        manifest_entries.append(rec_info)

    # Compute direct structural clash diagnostics across all crystal ligands vs 4X61(PRMT5•SAM) & 7S1S(PRMT5•MTA)
    sam_4x61 = json.loads((pre_dir / "4X61_cofactor_SAM.json").read_text())["atoms"]
    mta_7s1s = json.loads((pre_dir / "7S1S_cofactor_MTA.json").read_text())["atoms"]
    pkt_4x61 = json.loads((pre_dir / "4X61_pocket_atoms.json").read_text())["atoms"]
    pkt_7s1s = json.loads((pre_dir / "7S1S_pocket_atoms.json").read_text())["atoms"]

    sam_tail_names = ["N", "CA", "C", "O", "OXT", "CB", "CG", "SD", "CE"]
    sam_tail_pts = np.array([sam_4x61[k] for k in sam_tail_names], dtype=float)
    sam_ce_pt = np.array(sam_4x61["CE"], dtype=float)
    phe327_sam = np.array([a["xyz"] for a in pkt_4x61 if a["resname"] == "PHE" and a["resseq"] == 327], dtype=float)
    phe327_mta = np.array([a["xyz"] for a in pkt_7s1s if a["resname"] == "PHE" and a["resseq"] == 327], dtype=float)
    glu435_sam = np.array([a["xyz"] for a in pkt_4x61 if a["resname"] == "GLU" and a["resseq"] == 435], dtype=float)

    crystal_clash_table = []
    for entry in manifest_entries:
        if not entry.get("ligand_sdf"):
            continue
        suppl = Chem.SDMolSupplier(str(root / entry["ligand_sdf"]), removeHs=True)
        mol = suppl[0]
        pts = mol.GetConformer().GetPositions()
        d_sam_ce = float(np.linalg.norm(pts - sam_ce_pt[None, :], axis=1).min())
        d_sam_tail = float(np.linalg.norm(pts[:, None, :] - sam_tail_pts[None, :, :], axis=-1).min())
        d_phe327_sam = float(np.linalg.norm(pts[:, None, :] - phe327_sam[None, :, :], axis=-1).min())
        d_phe327_mta = float(np.linalg.norm(pts[:, None, :] - phe327_mta[None, :, :], axis=-1).min())
        d_glu435_sam = float(np.linalg.norm(pts[:, None, :] - glu435_sam[None, :, :], axis=-1).min())
        crystal_clash_table.append({
            "pdb_id": entry["pdb_id"],
            "ligand_id": entry["ligand_id"],
            "ligand_name": entry["ligand_name"],
            "state": entry["state"],
            "quarantined_post2023": entry["quarantined_post2023_diagnostic_only"],
            "min_dist_SAM_CE_A": round(d_sam_ce, 3),
            "min_dist_SAM_tail_all_A": round(d_sam_tail, 3),
            "min_dist_Glu435_SAM_state_A": round(d_glu435_sam, 3),
            "min_dist_Phe327_SAM_state_A": round(d_phe327_sam, 3),
            "min_dist_Phe327_MTA_state_A": round(d_phe327_mta, 3),
        })

    manifest_out = {
        "reference_frame": "PDB 4X61 chain A catalytic domain C-alpha residues 300-630 (331 atoms)",
        "featurization_cutoff_policy": "ONLY receptors deposited before 2023-01-01 (receptors_pre2023/) may be accessed by feature extraction or docking.",
        "receptors": manifest_entries,
        "crystallographic_steric_switch_audit": crystal_clash_table,
    }
    manifest_path = root / "manifests/pdb_receptors.json"
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(json.dumps(manifest_out, indent=2) + "\n", encoding="utf-8")
    print(f"Prepared {len(manifest_entries)} aligned PDB structures; manifest written to {manifest_path}")
    for row in crystal_clash_table:
        print(
            f"  {row['pdb_id']} ({row['ligand_name']:14s}, {row['state'][:6]}): "
            f"d(SAM:CE)={row['min_dist_SAM_CE_A']:.2f}A, d(Glu435_SAM)={row['min_dist_Glu435_SAM_state_A']:.2f}A, "
            f"d(Phe327_MTA)={row['min_dist_Phe327_MTA_state_A']:.2f}A"
        )


if __name__ == "__main__":
    main()
