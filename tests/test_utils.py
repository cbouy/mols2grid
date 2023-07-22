import gzip
from tempfile import NamedTemporaryFile
from types import SimpleNamespace
from unittest.mock import Mock, patch

import pandas as pd
import pytest
from rdkit import Chem
from rdkit.Chem.rdDepictor import Compute2DCoords

from mols2grid import utils


def test_requires():
    @utils.requires("_not_a_module")
    def func():
        pass

    with pytest.raises(
        ModuleNotFoundError, match="The module '_not_a_module' is required"
    ):
        func()


@pytest.mark.parametrize(
    ["subset", "fmt", "style", "transform", "exp"],
    [
        (
            ["SMILES", "ID"],
            "<strong>{key}</strong>: {value}",
            {},
            {},
            '<strong>SMILES</strong>: <span class="copy-me">CCO</span><br>'
            '<strong>ID</strong>: <span class="copy-me">0</span>',
        ),
        (["ID"], "foo-{value}", {}, {}, 'foo-<span class="copy-me">0</span>'),
        (
            ["ID"],
            "{value}",
            {"ID": lambda x: "color: red"},
            {},
            '<span class="copy-me" style="color: red">0</span>',
        ),
        (
            ["Activity"],
            "{value}",
            {},
            {"Activity": lambda x: f"{x:.2f}"},
            '<span class="copy-me">42.01</span>',
        ),
        (
            ["Activity"],
            "{key}: {value}",
            {"Activity": lambda x: "color: red" if x > 40 else ""},
            {"Activity": lambda x: f"{x:.2f}"},
            'Activity: <span class="copy-me" style="color: red">42.01</span>',
        ),
    ],
)
def test_tooltip_formatter(subset, fmt, style, transform, exp):
    row = pd.Series(
        {
            "ID": 0,
            "SMILES": "CCO",
            "Activity": 42.012345,
        }
    )
    tooltip = utils.tooltip_formatter(row, subset, fmt, style, transform)
    assert tooltip == exp


@pytest.mark.parametrize(
    ["smi", "exp"], [("CCO", "CCO"), ("blabla", None), (None, None)]
)
def test_mol_to_smiles(smi, exp):
    if smi:
        mol = Chem.MolFromSmiles(smi)
    else:
        mol = smi
    assert utils.mol_to_smiles(mol) == exp


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
    new = utils.mol_to_record(mol)
    assert "mol" in new.keys()
    new.pop("mol")
    assert new == props


def test_mol_to_record_none():
    new = utils.mol_to_record(None)
    assert new == {}


def test_mol_to_record_overwrite_smiles():
    mol = Chem.MolFromSmiles("CCO")
    mol.SetProp("SMILES", "foo")
    new = utils.mol_to_record(mol)
    assert new["SMILES"] == "foo"


def test_mol_to_record_custom_mol_col():
    mol = Chem.MolFromSmiles("CCO")
    new = utils.mol_to_record(mol, mol_col="foo")
    assert new["foo"] is mol


@pytest.mark.parametrize("sdf_source", ["sdf_path", "sdf_file"])
def test_sdf_to_dataframe(sdf_source, request):
    sdf = request.getfixturevalue(sdf_source)
    df = utils.sdf_to_dataframe(sdf)
    exp = {
        "ID": 5,
        "NAME": "3-methylpentane",
        "SMILES": "CCC(C)CC",
        "SOL": -3.68,
        "SOL_classification": "(A) low",
        "_MolFileComments": "",
        "_MolFileInfo": "  SciTegic05121109362D",
        "_Name": "3-methylpentane",
    }
    new = df.iloc[0].drop(["mol"]).to_dict()
    assert new == exp


def test_sdf_to_dataframe_custom_mol_col(sdf_path):
    df = utils.sdf_to_dataframe(sdf_path, mol_col="foo")
    assert "mol" not in df.columns
    assert "foo" in df.columns


def test_sdf_to_df_gz(sdf_path, tmp_path):
    tmp_file = tmp_path / "mols.sdf.gz"
    with open(sdf_path, "rb") as fi, open(tmp_file, "wb") as tf:
        gz = gzip.compress(fi.read(), compresslevel=1)
        tf.write(gz)
        tf.flush()
        df = utils.sdf_to_dataframe(tmp_file).drop(columns=["mol"])
        ref = utils.sdf_to_dataframe(sdf_path).drop(columns=["mol"])
        assert (df == ref).values.all()


def test_remove_coordinates():
    mol = Chem.MolFromSmiles("CCO")
    Compute2DCoords(mol)
    mol.GetConformer()
    new = utils.remove_coordinates(mol)
    assert new is mol
    with pytest.raises(ValueError, match="Bad Conformer Id"):
        new.GetConformer()


@pytest.mark.parametrize(
    ["string", "expected"],
    [
        ("Mol", "Mol"),
        ("mol name", "mol-name"),
        ("mol  name", "mol-name"),
        ("mol-name", "mol-name"),
        ("mol- name", "mol--name"),
        ("mol\tname", "mol-name"),
        ("mol\nname", "mol-name"),
        ("mol \t\n name", "mol-name"),
    ],
)
def test_slugify(string, expected):
    assert utils.slugify(string) == expected


@pytest.mark.parametrize("value", [1, 2])
def test_callback_handler(value):
    callback = lambda x: x + 1
    mock = Mock(side_effect=callback)
    event = SimpleNamespace(new=str(value))
    utils.callback_handler(mock, event)
    mock.assert_called_once_with(value)


def test_is_running_within_streamlit():
    assert utils.is_running_within_streamlit() is False
    with patch(
        "mols2grid.utils._get_streamlit_script_run_ctx",
        create=True,
        new=lambda: object(),
    ):
        assert utils.is_running_within_streamlit() is True
    with patch(
        "mols2grid.utils._get_streamlit_script_run_ctx", create=True, new=lambda: None
    ):
        assert utils.is_running_within_streamlit() is False
