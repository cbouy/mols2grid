import pytest
from numpy.testing import assert_equal
from tempfile import NamedTemporaryFile
from pathlib import Path
from rdkit import RDConfig, Chem
from rdkit.Chem import Draw
from rdkit.Chem.rdDepictor import Compute2DCoords
from mols2grid import MolGrid, sdf_to_dataframe
from mols2grid.select import register

_len_df = 5
sdf = f"{RDConfig.RDDocsDir}/Book/data/solubility.test.sdf"

def _make_df():
    return sdf_to_dataframe(sdf).iloc[:_len_df]

default_df = _make_df()

def _make_grid(**kwargs):
    kwargs["mol_col"] = kwargs.pop("mol_col", "mol")
    return MolGrid(default_df, **kwargs)

@pytest.fixture(scope="module")
def df():
    return _make_df()

@pytest.fixture(scope="module")
def grid():
    return _make_grid()

@pytest.fixture(scope="module")
def mols():
    return default_df["mol"]

def test_no_input_specified():
    with pytest.raises(ValueError,
                       match="One of `smiles_col` or `mol_col` must be set"):
        _make_grid(smiles_col=None, mol_col=None)

def test_name_not_str():
    with pytest.raises(TypeError,
                       match="`name` must be a string. Currently of type int"):
        _make_grid(name=0)

def test_uses_svg(grid):
    assert grid.useSVG == True
    data = grid.dataframe.loc[0, "img"]
    assert "<svg " in data[:100]

def test_uses_png():
    grid = _make_grid(useSVG=False)
    assert grid.useSVG == False
    data = grid.dataframe.loc[0, "img"]
    assert data.startswith('<img src="data:image/png;base64')

def test_copy_df(df):
    grid = MolGrid(df, mol_col="mol")
    assert grid.dataframe is not df

def test_records_as_input():
    data = [{"SMILES": "C"*i, "ID": i} for i in range(1, 5)]
    MolGrid(data)

def test_mapping_deprecated():
    with pytest.warns(UserWarning, match="`mapping` is deprecated"):
        _make_grid(mapping={"SOL": "Solubility"})

def test_rename():
    grid = _make_grid(rename={"SOL": "Solubility"})
    assert "Solubility" in grid.dataframe.columns

def test_smi_only_makes_mol_temporary(df):
    df = df.copy()
    df["SMILES"] = df["mol"].apply(Chem.MolToSmiles)
    df = df.drop(columns=["mol"])
    grid = MolGrid(df)
    assert "mol" in grid._extra_columns
    assert "mol" not in grid.dataframe.columns

def test_keep_mol_col_if_input(grid):
    assert "mol" not in grid._extra_columns
    assert "mol" in grid.dataframe.columns

def test_remove_hydrogens():
    mol =  Chem.MolFromSmiles("C")
    mol = Chem.AddHs(mol)
    data = [{"mol": mol, "ID": 0}]
    grid = MolGrid(data, mol_col="mol", removeHs=True)
    new = grid.dataframe.loc[0, "mol"]
    assert new.GetNumAtoms() == 1

def test_remove_coordinates():
    mol = Chem.MolFromSmiles("C")
    mol = Chem.AddHs(mol)
    Compute2DCoords(mol)
    data = [{"mol": mol, "ID": 0}]
    grid = MolGrid(data, mol_col="mol", use_coords=False)
    new = grid.dataframe.loc[0, "mol"]
    with pytest.raises(ValueError, match="Bad Conformer Id"):
        new.GetConformer()

def test_generate_smi_from_mol_col(df):
    df = df.drop(columns=["SMILES"])
    grid = MolGrid(df, mol_col="mol")
    assert "SMILES" in grid.dataframe.columns

def test_index_when_none_molecules(df):
    df = df.copy()
    df.loc[1, "mol"] = None
    grid = MolGrid(df, mol_col="mol")
    assert len(grid.dataframe) == len(df) - 1
    assert_equal(grid.dataframe["mols2grid-id"].values, [0, 2, 3, 4])

def test_custom_moldrawoptions():
    opts = Draw.MolDrawOptions()
    opts.fixedBondLength = 42
    grid = _make_grid(MolDrawOptions=opts)
    assert grid.MolDrawOptions.fixedBondLength == 42

def test_moldrawoptions_as_kwargs():
    grid = _make_grid(fixedBondLength=42)
    assert grid.MolDrawOptions.fixedBondLength == 42

def test_update_current_grid_on_init():
    _make_grid(name="foo")
    assert register.current_selection == "foo"
    assert register.get_selection("foo") is register.get_selection()
    assert register.get_selection() == {}

def test_overwrite_warning():
    _make_grid()
    register._set_selection(0, "C")
    with pytest.warns(
        UserWarning,
        match=f"Overwriting non-empty 'default' grid selection: {{0: 'C'}}"
    ):
        _make_grid()
    assert register.get_selection() == {}

def test_from_mols(mols):
    grid = MolGrid.from_mols(mols)
    assert "mol" in grid.dataframe.columns

def test_from_mols_custom_mol_col(mols):
    grid = MolGrid.from_mols(mols, mol_col="rdmol")
    assert "rdmol" in grid.dataframe.columns
    assert "mol" not in grid.dataframe.columns

@pytest.mark.parametrize(["param", "attr", "value"], [
    ("useSVG", "useSVG", False),
    ("size", "img_size", (42, 42)),
    ("name", "_grid_id", "foobar"),
])
def test_from_mols_kwargs(mols, param, attr, value):
    grid = MolGrid.from_mols(mols, **{param: value})
    assert getattr(grid, attr) == value

def test_from_sdf():
    grid = MolGrid.from_sdf(sdf)
    assert "mol" in grid.dataframe.columns

def test_from_sdf_custom_mol_col():
    grid = MolGrid.from_sdf(sdf, mol_col="rdmol")
    assert "rdmol" in grid.dataframe.columns
    assert "mol" not in grid.dataframe.columns

@pytest.mark.parametrize(["param", "attr", "value"], [
    ("useSVG", "useSVG", False),
    ("size", "img_size", (42, 42)),
    ("name", "_grid_id", "foobar"),
])
def test_from_sdf_kwargs(param, attr, value):
    grid = MolGrid.from_sdf(sdf, **{param: value})
    assert getattr(grid, attr) == value

def test_template(grid):
    grid.template = "pages"
    assert grid._template is grid.template

def test_template_setter(grid):
    with pytest.raises(ValueError, match="template='foo' not supported"):
        grid.template = "foo"

def test_substruct_highlight():
    mol = Chem.MolFromSmiles("C")
    mol.__sssAtoms = [0]
    grid = MolGrid.from_mols([mol])
    svg = grid.dataframe.loc[0, "img"]
    assert "<ellipse" in svg

def test_mol_to_img_svg():
    mol = Chem.MolFromSmiles("C")
    grid = MolGrid.from_mols([mol])
    img = grid.mol_to_img(mol)
    assert "<svg " in img[:100]

def test_mol_to_img_png():
    mol = Chem.MolFromSmiles("C")
    grid = MolGrid.from_mols([mol])
    grid.useSVG = False
    grid._MolDraw2D = Draw.MolDraw2DCairo
    img = grid.mol_to_img(mol)
    assert img.startswith('<img src="data:image/png;base64')

def test_get_selection(df):
    grid = MolGrid(df, mol_col="mol")
    register._set_selection(0, "")
    new = grid.get_selection()
    assert_equal(new.values,
                 df.iloc[0:1].values)

def test_save(grid):
    with NamedTemporaryFile("w", suffix=".html") as f:
        grid.save(f.name)
        assert Path(f.name).is_file()

# TODO: test filters and renderers