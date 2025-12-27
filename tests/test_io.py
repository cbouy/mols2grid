from pathlib import Path

import pandas as pd
import pytest

from mols2grid.io import read_mols_to_df, read_records
from mols2grid.utils import import_object

FIRST_MOL_DATA = {
    "SMILES": "CCC(C)CC",
    "_Name": "3-methylpentane",
    "_MolFileInfo": "SciTegic05121109362D",
    "_MolFileComments": "",
    "_MolFileChiralFlag": 0,
    "ID": 5,
    "NAME": "3-methylpentane",
    "SOL": -3.68,
    "SOL_classification": "(A) low",
}


@pytest.fixture
def csv_path(df: pd.DataFrame, tmp_path: Path):
    tmp_file = tmp_path / "mols.csv"
    df.drop("mol", axis=1).to_csv(tmp_file, index=False)
    return tmp_file


@pytest.fixture
def smi_path(df: pd.DataFrame, tmp_path: Path):
    tmp_file = tmp_path / "mols.smi"
    df.drop("mol", axis=1).to_csv(tmp_file, sep="\t", index=False)
    return tmp_file


@pytest.mark.parametrize("sdf_source", ["sdf_path", "sdf_file"])
def test_read_records(sdf_source, request):
    sdf = request.getfixturevalue(sdf_source)
    records = read_records(sdf)
    new = next(records)
    new["_MolFileInfo"] = new["_MolFileInfo"].strip()
    new.pop("mol")
    assert new == FIRST_MOL_DATA


def test_read_to_dataframe(sdf_path):
    df = read_mols_to_df(sdf_path)
    new = df.iloc[0].drop(["mol"]).to_dict()
    new["_MolFileInfo"] = new["_MolFileInfo"].strip()
    assert new == FIRST_MOL_DATA


def test_sdf_to_dataframe_custom_mol_col(sdf_path):
    df = read_mols_to_df(sdf_path, mol_col="foo")
    assert "mol" not in df.columns
    assert "foo" in df.columns


@pytest.mark.parametrize(
    ("ext", "compress_path", "kwargs"),
    [
        ("gz", "gzip.compress", {"compresslevel": 1}),
        ("xz", "lzma.compress", {}),
        ("bz2", "bz2.compress", {"compresslevel": 1}),
    ],
)
def test_compressed(ext, compress_path, kwargs, sdf_path, tmp_path):
    tmp_file = tmp_path / f"mols.sdf.{ext}"
    with open(tmp_file, "wb") as tf:
        compress = import_object(compress_path)
        dump = compress(sdf_path.read_bytes(), **kwargs)
        tf.write(dump)
        tf.flush()
        df = read_mols_to_df(tmp_file).drop(columns=["mol"])
        ref = read_mols_to_df(sdf_path).drop(columns=["mol"])
        assert (df == ref).to_numpy().all()


@pytest.mark.parametrize("source", ["csv_path", "smi_path"])
def test_custom_readers(request: pytest.FixtureRequest, source):
    path = request.getfixturevalue(source)
    records = read_records(path)
    new = next(records)
    new["_MolFileInfo"] = new["_MolFileInfo"].strip()
    new.pop("mol")
    assert new == {k: str(v) for k, v in FIRST_MOL_DATA.items()}


def test_smi_no_header(tmp_path: Path):
    df = pd.DataFrame(
        {"SMILES": ["CCC(C)CC", "CCO"], "ID": [1, 2], "NAME": ["foo", "bar"]}
    )
    tmp_file = tmp_path / "mols.smi"
    df.to_csv(tmp_file, sep=" ", index=False, header=False)
    records = read_records(tmp_file)
    new = next(records)
    new.pop("mol")
    assert new == {"SMILES": "CCC(C)CC", "TITLE": "1", "field_0": "foo"}


def test_invalid_smi_delimiter(tmp_path: Path):
    tmp_file = tmp_path / "mols.smi"
    with open(tmp_file, "w") as tf:
        tf.write("CC;1\nCC;2")
    records = read_records(tmp_file)
    with pytest.raises(ValueError, match="No valid columns found"):
        next(records)


def test_delimiter(tmp_path: Path):
    tmp_file = tmp_path / "mols.smi"
    with open(tmp_file, "w") as tf:
        tf.write("CC;1\nCC;2")
    records = read_records(tmp_file, delimiter=";")
    new = next(records)
    new.pop("mol")
    assert new == {"SMILES": "CC", "TITLE": "1"}


def test_smi_fieldnames(tmp_path: Path):
    tmp_file = tmp_path / "mols.smi"
    with open(tmp_file, "w") as tf:
        tf.write("CC 1\nCC 2")
    records = read_records(tmp_file, fieldnames=["SMILES", "id"])
    new = next(records)
    new.pop("mol")
    assert new == {"SMILES": "CC", "id": "1"}


def test_unknown_reader():
    with pytest.raises(ValueError, match="No reader found"):
        next(read_records("mols.bin"))
