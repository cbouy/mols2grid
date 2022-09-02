import pytest
from numpy.testing import assert_equal
from tempfile import NamedTemporaryFile
from pathlib import Path
from types import SimpleNamespace
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
    kwargs.setdefault("mol_col", "mol")
    return MolGrid(default_df, **kwargs)

@pytest.fixture(scope="module")
def df():
    return _make_df()

@pytest.fixture(scope="module")
def grid():
    return _make_grid(prerender=True)

@pytest.fixture(scope="module")
def grid_otf():
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
    grid = _make_grid(useSVG=False, prerender=True)
    assert grid.useSVG == False
    data = grid.dataframe.loc[0, "img"]
    assert data.startswith('<img src="data:image/png;base64')

def test_copy_df(df):
    grid = MolGrid(df, mol_col="mol")
    assert grid.dataframe is not df

def test_records_as_input():
    data = [{"SMILES": "C"*i, "ID": i} for i in range(1, 5)]
    MolGrid(data)

def test_rename():
    grid = _make_grid(rename={"SOL": "Solubility"})
    assert "Solubility" in grid.dataframe.columns

def test_smi_only_makes_mol_temporary(df):
    df = df.copy()
    df["SMILES"] = df["mol"].apply(Chem.MolToSmiles)
    df = df.drop(columns=["mol"])
    grid = MolGrid(df)
    assert "mol" not in grid._extra_columns
    assert "mol" not in grid.dataframe.columns

def test_keep_mol_col_if_input(grid):
    assert "mol" not in grid._extra_columns
    assert "mol" in grid.dataframe.columns

def test_remove_hydrogens():
    mol =  Chem.MolFromSmiles("C")
    mol = Chem.AddHs(mol)
    data = [{"mol": mol, "ID": 0}]
    grid = MolGrid(data, mol_col="mol", removeHs=True, prerender=True)
    new = grid.dataframe.loc[0, "mol"]
    assert new.GetNumAtoms() == 1

def test_remove_coordinates():
    mol = Chem.MolFromSmiles("C")
    mol = Chem.AddHs(mol)
    Compute2DCoords(mol)
    data = [{"mol": mol, "ID": 0}]
    grid = MolGrid(data, mol_col="mol", use_coords=False, prerender=True)
    new = grid.dataframe.loc[0, "mol"]
    with pytest.raises(ValueError, match="Bad Conformer Id"):
        new.GetConformer()

def test_generate_smi_from_mol_col(df):
    df = df.drop(columns=["SMILES"])
    grid = MolGrid(df, mol_col="mol", prerender=True)
    assert "SMILES" in grid.dataframe.columns

def test_index_when_none_molecules(df):
    df = df.copy()
    df.loc[1, "mol"] = None
    grid = MolGrid(df, mol_col="mol", prerender=True)
    assert len(grid.dataframe) == len(df) - 1
    assert_equal(grid.dataframe["mols2grid-id"].values, [0, 2, 3, 4])

def test_custom_moldrawoptions():
    opts = Draw.MolDrawOptions()
    opts.fixedBondLength = 42
    grid = _make_grid(MolDrawOptions=opts, prerender=True)
    assert grid.MolDrawOptions.fixedBondLength == 42

def test_moldrawoptions_as_kwargs():
    grid = _make_grid(fixedBondLength=42, prerender=True)
    assert grid.MolDrawOptions.fixedBondLength == 42

def test_update_current_grid_on_init():
    _make_grid(name="foo")
    assert register.current_selection == "foo"
    assert register.get_selection("foo") is register.get_selection()
    assert register.get_selection() == {}

def test_from_mols(mols):
    grid = MolGrid.from_mols(mols)
    assert "mol" in grid.dataframe.columns

def test_from_mols_custom_mol_col(mols):
    grid = MolGrid.from_mols(mols, mol_col="rdmol")
    assert "rdmol" in grid.dataframe.columns
    assert "mol" not in grid.dataframe.columns

@pytest.mark.parametrize(["param", "attr", "value"], [
    ("name", "_grid_id", "foobar"),
    ("coordGen", "prefer_coordGen", False),
    ("removeHs", "removeHs", True),
    ("useSVG", "useSVG", False),
    ("use_coords", "use_coords", False),
    ("size", "img_size", (42, 42)),
    ("prerender", "prerender", True),
    ("smiles_col", "smiles_col", "smi"),
    ("mol_col", "mol_col", "rdmol"),
])
def test_from_mols_kwargs(mols, param, attr, value):
    if param in ["useSVG"]:
        grid = MolGrid.from_mols(mols, prerender=True, **{param: value})
    else:
        grid = MolGrid.from_mols(mols, **{param: value})
    assert getattr(grid, attr) == value

def test_from_sdf():
    grid = MolGrid.from_sdf(sdf)
    assert "mol" in grid.dataframe.columns

def test_from_sdf_custom_mol_col():
    grid = MolGrid.from_sdf(sdf, mol_col="rdmol", prerender=True)
    assert "rdmol" in grid.dataframe.columns
    assert "mol" not in grid.dataframe.columns

@pytest.mark.parametrize(["param", "attr", "value"], [
    ("name", "_grid_id", "foobar"),
    ("coordGen", "prefer_coordGen", False),
    ("removeHs", "removeHs", True),
    ("useSVG", "useSVG", False),
    ("use_coords", "use_coords", False),
    ("size", "img_size", (42, 42)),
    ("prerender", "prerender", True),
    ("smiles_col", "smiles_col", "smi"),
    ("mol_col", "mol_col", "rdmol"),
])
def test_from_sdf_kwargs(param, attr, value):
    if param in ["useSVG"]:
        grid = MolGrid.from_sdf(sdf, prerender=True, **{param: value})
    else:
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
    grid = MolGrid.from_mols([mol], prerender=True)
    svg = grid.dataframe.loc[0, "img"]
    assert "<ellipse" in svg

def test_mol_to_img_svg():
    mol = Chem.MolFromSmiles("C")
    grid = MolGrid.from_mols([mol], prerender=True)
    img = grid.mol_to_img(mol)
    assert "<svg " in img[:100]

def test_mol_to_img_png():
    mol = Chem.MolFromSmiles("C")
    grid = MolGrid.from_mols([mol], prerender=True)
    grid.useSVG = False
    grid._MolDraw2D = Draw.MolDraw2DCairo
    img = grid.mol_to_img(mol)
    assert img.startswith('<img src="data:image/png;base64')

def test_get_selection(df):
    grid = MolGrid(df, mol_col="mol", name="grid")
    other = MolGrid(df, mol_col="mol", name="other")
    event = SimpleNamespace(new='{0: ""}')
    register.selection_updated("grid", event)
    assert register.current_selection == "grid"
    assert grid.get_selection().equals(df.head(1))
    assert other.get_selection().equals(df.head(0))  # empty dataframe
    register._clear()

def test_save(grid_otf):
    with NamedTemporaryFile("w", suffix=".html") as f:
        grid_otf.save(f.name)
        assert Path(f.name).is_file()

def test_render_pages(grid_otf):
    grid_otf.render(template="pages")

def test_render_table(grid):
    grid.render(template="table")

def test_render_wrong_template(grid_otf):
    with pytest.raises(ValueError, match="template='foo' not supported"):
        grid_otf.render(template="foo")

@pytest.mark.parametrize("kwargs", [
    dict(),
    dict(subset=["ID"]),
    dict(tooltip=["ID"]),
    dict(gap="5px"),
    dict(style={"ID": lambda x: "color: red" if x == 1 else ""}),
    dict(transform={"ID": lambda x: f"Id. #{x}"}),
])
def test_integration_table(grid, kwargs):
    grid.to_table(**kwargs)

def test_python_callback(grid_otf):
    html = grid_otf.to_pages(subset=["img"], callback=lambda data: None)
    assert "// trigger custom python callback" in html
    assert "// no kernel detected for callback" in html

def test_cache_selection():
    grid = _make_grid(name="cache")
    event = SimpleNamespace(new='{0: "CCO"}')
    register.selection_updated("cache", event)
    grid = _make_grid(name="cache", cache_selection=True)
    assert hasattr(grid, "_cached_selection")
    assert register.get_selection("cache") == grid._cached_selection

def test_cache_no_init():
    grid = _make_grid(name="new_cache", cache_selection=True)
    assert hasattr(grid, "_cached_selection")
    assert grid._cached_selection == {}
    assert "new_cache" in register.list_grids()

def test_no_cache_selection():
    grid = _make_grid(name="no_cache", cache_selection=False)
    assert grid._cached_selection == {}
    assert "no_cache" in register.list_grids()

def test_onthefly_render_png_error(df):
    with pytest.raises(ValueError,
                       match="On-the-fly rendering of PNG images not supported"):
        MolGrid(df, prerender=False, useSVG=False)

def test_substruct_highlight_prerender_error(grid):
    with pytest.raises(
        ValueError,
        match="Cannot highlight substructure search with prerendered images"
    ):
        grid.display(substruct_highlight=True)
    # disables substruct highlight when prerendering by default
    grid.display()

def test_use_coords_onthefly_error():
    with pytest.raises(ValueError,
                       match="Cannot use coordinates with on-the-fly rendering"):
        _make_grid(use_coords=True, prerender=False)

def test_sort_by_not_in_subset_or_tooltip(grid_otf):
    with pytest.raises(ValueError,
                       match="'_Name' is not an available field"):
        grid_otf.to_pages(subset=["ID", "img"], sort_by="_Name")

def test_subset_without_img_error(grid_otf):
    with pytest.raises(KeyError,
                       match="Please add the 'img' field in the `subset` parameter"):
        grid_otf.display(subset=["_Name"])

def test_table_no_prerender_error(grid_otf):
    with pytest.raises(ValueError,
        match="Please set `prerender=True` when using the 'table' template"
    ):
        grid_otf.display(template="table")

def test_replace_non_serializable_from_default_output():
    g =_make_grid()
    g.dataframe["non-serializable"] = g.dataframe[g.mol_col]
    html = g.render()
    assert "\\ud83e\\udd37\\u200d\\u2642\\ufe0f" in html  #🤷‍♂️
