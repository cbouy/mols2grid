from pathlib import Path

import pytest
from rdkit import RDConfig

from mols2grid import MolGrid, sdf_to_dataframe


@pytest.fixture(scope="module")
def sdf_path():
    return Path(RDConfig.RDDocsDir) / "Book" / "data" / "solubility.test.sdf"


@pytest.fixture(scope="module")
def sdf_file(sdf_path):
    return str(sdf_path)


@pytest.fixture(scope="module")
def smiles_records():
    return [{"SMILES": "C" * i, "ID": i} for i in range(1, 5)]


@pytest.fixture(scope="module")
def df(sdf_path):
    return sdf_to_dataframe(sdf_path).head(30)


@pytest.fixture(scope="module")
def small_df(df):
    return df.head(5)


@pytest.fixture(scope="module")
def grid_prerendered(df):
    return MolGrid(df, mol_col="mol", prerender=True)


@pytest.fixture(scope="module")
def grid(df):
    return MolGrid(df, mol_col="mol", size=(160, 120))


@pytest.fixture(scope="module")
def mols(small_df):
    return small_df["mol"]
