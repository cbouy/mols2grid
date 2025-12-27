from typing import Any

from rdkit import Chem


def mol_to_smiles(mol: Chem.Mol | None) -> str | None:
    """Returns a SMILES from an RDKit molecule, or None if not an RDKit mol"""
    return Chem.MolToSmiles(mol) if mol else None


def mol_to_record(mol: Chem.Mol | None, mol_col: str = "mol") -> dict[str, Any]:
    """Function to create a dict of data from an RDKit molecule"""
    return {**mol.GetPropsAsDict(includePrivate=True), mol_col: mol} if mol else {}


def remove_coordinates(mol: Chem.Mol) -> Chem.Mol:
    """Removes the existing coordinates from the molecule. The molecule is
    modified inplace"""
    mol.RemoveAllConformers()
    return mol
