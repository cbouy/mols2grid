from types import SimpleNamespace

import pytest
from numpy.testing import assert_equal
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.rdDepictor import Compute2DCoords

from mols2grid import MolGrid
from mols2grid.select import register


def get_grid(df, **kwargs):
    kwargs.setdefault("mol_col", "mol")
    return MolGrid(df, **kwargs)


def test_no_input_specified(small_df):
    with pytest.raises(
        ValueError, match="One of `smiles_col` or `mol_col` must be set"
    ):
        get_grid(small_df, smiles_col=None, mol_col=None)


def test_name_not_str(small_df):
    with pytest.raises(
        TypeError, match="`name` must be a string. Currently of type int"
    ):
        get_grid(small_df, name=0)


def test_uses_svg(grid_prerendered):
    assert grid_prerendered.useSVG == True
    data = grid_prerendered.dataframe.loc[0, "img"]
    assert "<svg " in data[:100]


def test_uses_png(small_df):
    grid = get_grid(small_df, useSVG=False, prerender=True)
    assert grid.useSVG == False
    data = grid.dataframe.loc[0, "img"]
    assert data.startswith('<img src="data:image/png;base64')


def test_copy_df(df):
    grid = MolGrid(df, mol_col="mol")
    assert grid.dataframe is not df


def test_records_as_input(smiles_records):
    MolGrid(smiles_records)


def test_rename(small_df):
    grid = get_grid(small_df, rename={"SOL": "Solubility"})
    assert "Solubility" in grid.dataframe.columns


def test_smi_only_makes_mol_temporary(small_df):
    df = small_df.copy()
    df["SMILES"] = df["mol"].apply(Chem.MolToSmiles)
    df = df.drop(columns=["mol"])
    grid = MolGrid(df)
    assert "mol" not in grid._extra_columns
    assert "mol" not in grid.dataframe.columns


def test_keep_mol_col_if_input(grid_prerendered):
    assert "mol" not in grid_prerendered._extra_columns
    assert "mol" in grid_prerendered.dataframe.columns


def test_remove_hydrogens():
    mol = Chem.MolFromSmiles("C")
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


def test_generate_smi_from_mol_col(small_df):
    df = small_df.drop(columns=["SMILES"])
    grid = MolGrid(df, mol_col="mol", prerender=True)
    assert "SMILES" in grid.dataframe.columns


def test_index_when_none_molecules(small_df):
    df = small_df.copy()
    df.loc[1, "mol"] = None
    grid = MolGrid(df, mol_col="mol", prerender=True)
    assert len(grid.dataframe) == len(df) - 1
    assert_equal(grid.dataframe["mols2grid-id"].values, [0, 2, 3, 4])


def test_custom_moldrawoptions(small_df):
    opts = Draw.MolDrawOptions()
    opts.fixedBondLength = 42
    grid = get_grid(small_df, MolDrawOptions=opts, prerender=True)
    assert grid.MolDrawOptions.fixedBondLength == 42


def test_moldrawoptions_as_kwargs(small_df):
    grid = get_grid(small_df, fixedBondLength=42, prerender=True)
    assert grid.MolDrawOptions.fixedBondLength == 42


def test_update_current_grid_on_init(small_df):
    get_grid(small_df, name="foo")
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


@pytest.mark.parametrize(
    ["param", "attr", "value"],
    [
        ("name", "_grid_id", "foobar"),
        ("coordGen", "prefer_coordGen", False),
        ("removeHs", "removeHs", True),
        ("useSVG", "useSVG", False),
        ("use_coords", "use_coords", False),
        ("size", "img_size", (42, 42)),
        ("prerender", "prerender", True),
        ("smiles_col", "smiles_col", "smi"),
        ("mol_col", "mol_col", "rdmol"),
    ],
)
def test_from_mols_kwargs(mols, param, attr, value):
    if param in ["useSVG"]:
        grid = MolGrid.from_mols(mols, prerender=True, **{param: value})
    else:
        grid = MolGrid.from_mols(mols, **{param: value})
    assert getattr(grid, attr) == value


@pytest.mark.parametrize("sdf_source", ["sdf_path", "sdf_file"])
def test_from_sdf(sdf_source, request):
    sdf = request.getfixturevalue(sdf_source)
    grid = MolGrid.from_sdf(sdf)
    assert "mol" in grid.dataframe.columns


def test_from_sdf_custom_mol_col(sdf_path):
    grid = MolGrid.from_sdf(sdf_path, mol_col="rdmol", prerender=True)
    assert "rdmol" in grid.dataframe.columns
    assert "mol" not in grid.dataframe.columns


@pytest.mark.parametrize(
    ["param", "attr", "value"],
    [
        ("name", "_grid_id", "foobar"),
        ("coordGen", "prefer_coordGen", False),
        ("removeHs", "removeHs", True),
        ("useSVG", "useSVG", False),
        ("use_coords", "use_coords", False),
        ("size", "img_size", (42, 42)),
        ("prerender", "prerender", True),
        ("smiles_col", "smiles_col", "smi"),
        ("mol_col", "mol_col", "rdmol"),
    ],
)
def test_from_sdf_kwargs(sdf_path, param, attr, value):
    if param in ["useSVG"]:
        grid = MolGrid.from_sdf(sdf_path, prerender=True, **{param: value})
    else:
        grid = MolGrid.from_sdf(sdf_path, **{param: value})
    assert getattr(grid, attr) == value


def test_template(grid_prerendered):
    grid_prerendered.template = "interactive"
    assert grid_prerendered._template is grid_prerendered.template


def test_template_setter(grid_prerendered):
    with pytest.raises(ValueError, match="template='foo' not supported"):
        grid_prerendered.template = "foo"


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


def test_get_selection(small_df):
    grid = MolGrid(small_df, mol_col="mol", name="grid")
    other = MolGrid(small_df, mol_col="mol", name="other")
    event = SimpleNamespace(new='{0: ""}')
    register.selection_updated("grid", event)
    assert register.current_selection == "grid"
    assert grid.get_selection().equals(small_df.head(1))
    assert other.get_selection().equals(small_df.head(0))  # empty dataframe
    register._clear()


def test_save(grid, tmp_path):
    output = tmp_path / "output.html"
    grid.save(output)
    assert output.is_file()


def test_render_interactive(grid):
    grid.render(template="interactive")


def test_render_static(grid_prerendered):
    grid_prerendered.render(template="static")


def test_render_wrong_template(grid):
    with pytest.raises(ValueError, match="template='foo' not supported"):
        grid.render(template="foo")


@pytest.mark.parametrize(
    "kwargs",
    [
        dict(),
        dict(subset=["ID"]),
        dict(tooltip=["ID"]),
        dict(gap=5),
        dict(style={"ID": lambda x: "color: red" if x == 1 else ""}),
        dict(transform={"ID": lambda x: f"Id. #{x}"}),
    ],
)
def test_integration_static(grid_prerendered, kwargs):
    grid_prerendered.to_static(**kwargs)


def test_python_callback(grid):
    html = grid.to_interactive(subset=["img"], callback=lambda data: None)
    assert "// Trigger custom python callback" in html
    assert "// No kernel detected for callback" in html


def test_cache_selection(small_df):
    grid = get_grid(small_df, name="cache")
    event = SimpleNamespace(new='{0: "CCO"}')
    register.selection_updated("cache", event)
    grid = get_grid(small_df, name="cache", cache_selection=True)
    assert hasattr(grid, "_cached_selection")
    assert register.get_selection("cache") == grid._cached_selection


def test_cache_no_init(small_df):
    grid = get_grid(small_df, name="new_cache", cache_selection=True)
    assert hasattr(grid, "_cached_selection")
    assert grid._cached_selection == {}
    assert "new_cache" in register.list_grids()


def test_no_cache_selection(small_df):
    grid = get_grid(small_df, name="no_cache", cache_selection=False)
    assert grid._cached_selection == {}
    assert "no_cache" in register.list_grids()


def test_onthefly_render_png_error(df):
    with pytest.raises(
        ValueError, match="On-the-fly rendering of PNG images not supported"
    ):
        MolGrid(df, prerender=False, useSVG=False)


def test_substruct_highlight_prerender_error(grid_prerendered):
    with pytest.raises(
        ValueError, match="Cannot highlight substructure search with prerendered images"
    ):
        grid_prerendered.display(substruct_highlight=True)
    # disables substruct highlight when prerendering by default
    grid_prerendered.display()


def test_use_coords_onthefly_error(small_df):
    with pytest.raises(
        ValueError, match="Cannot use coordinates with on-the-fly rendering"
    ):
        get_grid(small_df, use_coords=True, prerender=False)


def test_sort_by_not_in_subset_or_tooltip(grid):
    with pytest.raises(ValueError, match="'_Name' is not an available field"):
        grid.to_interactive(subset=["ID", "img"], sort_by="_Name")


def test_subset_without_img_no_error(grid):
    grid.display(subset=["_Name"])


def test_static_no_prerender_error(grid):
    with pytest.raises(
        ValueError, match="Please set `prerender=True` when using the 'static' template"
    ):
        grid.display(template="static")


def test_replace_non_serializable_from_default_output(small_df):
    grid = get_grid(small_df)
    grid.dataframe["non-serializable"] = grid.dataframe[grid.mol_col]
    html = grid.render()
    assert "\\ud83e\\udd37\\u200d\\u2642\\ufe0f" in html  # 🤷‍♂️
