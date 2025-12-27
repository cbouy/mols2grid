import pytest
from rdkit import Chem
from rdkit.Chem.rdDepictor import Compute2DCoords

from mols2grid.chem import mol_to_record, mol_to_smiles, remove_coordinates


@pytest.mark.parametrize(
    ("smi", "exp"), [("CCO", "CCO"), ("blabla", None), (None, None)]
)
def test_mol_to_smiles(smi, exp):
    mol = Chem.MolFromSmiles(smi) if smi else smi
    assert mol_to_smiles(mol) == exp


def test_mol_to_record():
    mol = Chem.MolFromSmiles("CCO")
    props = {
        "NAME": "ethanol",
        "foo": 42,
        "_bar": 42.01,
        "__baz": 0,
    }
    for prop, value in props.items():
        if isinstance(value, int):
            mol.SetIntProp(prop, value)
        elif isinstance(value, float):
            mol.SetDoubleProp(prop, value)
        else:
            mol.SetProp(prop, value)
    new = mol_to_record(mol)
    assert "mol" in new
    new.pop("mol")
    assert new == props


def test_mol_to_record_none():
    new = mol_to_record(None)
    assert new == {}


def test_mol_to_record_overwrite_smiles():
    mol = Chem.MolFromSmiles("CCO")
    mol.SetProp("SMILES", "foo")
    new = mol_to_record(mol)
    assert new["SMILES"] == "foo"


def test_mol_to_record_custom_mol_col():
    mol = Chem.MolFromSmiles("CCO")
    new = mol_to_record(mol, mol_col="foo")
    assert new["foo"] is mol


def test_remove_coordinates():
    mol = Chem.MolFromSmiles("CCO")
    Compute2DCoords(mol)
    mol.GetConformer()
    new = remove_coordinates(mol)
    assert new is mol
    with pytest.raises(ValueError, match="Bad Conformer Id"):
        new.GetConformer()
