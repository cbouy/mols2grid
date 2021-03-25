from functools import singledispatch
import inspect
from pandas import DataFrame, Series
from .molgrid import MolGrid


_SIGNATURE = {
    method: dict(inspect.signature(getattr(MolGrid, method)).parameters.items())
    for method in ["render", "to_pages", "to_table"]
}
for method in ["render", "to_pages", "to_table"]:
    _SIGNATURE[method].pop("self")
_SIGNATURE["render"].pop("kwargs")


def _prepare_kwargs(kwargs, kind):
    template = kwargs.pop("template", _SIGNATURE["render"]["template"].default)
    render_kwargs = {key: value
                     for key, value in kwargs.items()
                     if key in _SIGNATURE[f'to_{template}'].keys()}
    for key in render_kwargs.keys():
        kwargs.pop(key)
    if kind == "display":
        render_kwargs = {"width": kwargs.pop("width", "100%"),
                         "height": kwargs.pop("height", None),
                         **render_kwargs}
    elif kind == "save":
        render_kwargs = {"output": kwargs.pop("output"),
                         **render_kwargs}
    return template, kwargs, render_kwargs

@singledispatch
def display(arg, **kwargs):
    """Display molecules on an interactive grid
    
    Parameters
    ----------
    arg : pandas.DataFrame, SDF file or list of molecules
        The input containing your molecules
    smiles_col : str or None ("SMILES")
        If a pandas DataFrame is used, name of the column with SMILES
    mol_col : str or None (None)
        If a pandas DataFrame is used, name of the column with RDKit molecules
    useSVG : bool (True)
        Use SVG images or PNG
    coordGen : bool (True)
        Use the coordGen library instead of the RDKit one to depict the
        molecules in 2D
    use_coords : bool
        Use the coordinates of the molecules (only relevant when an SDF file, a
        list of molecules or a DataFrame of RDKit molecules were used as input)
    remove_Hs : bool
        Remove hydrogen atoms from the drawings
    size : tuple ((160, 120))
        Size of each image
    mapping : dict (None)
        Rename the properties in the final document
    template : str ("pages")
        One of {pages, table}
    width : str ("100%)
        Width of the frame displayed in the notebook
    height : int or None (None)
        Height of the frame displayed in the notebook. Use None for an
        automatic guess
    subset: list (None)
        Columns to be displayed in each cell of the grid. Each column's value
        will be displayed from top to bottom in the same order given here.
        Use "img" for the image of the molecule.
        Default: all columns (with "img" in first position)
    tooltip : list (None)
        Columns to be displayed as a tooltip when hovering/clicking on the
        image of a cell. Use `None` for no tooltip.
    tooltip_fmt : str ("<strong>{key}</strong>: {value}")
        Format string of each key/value pair in the tooltip
    tooltip_trigger : str ("click hover")
        Sequence of triggers for the tooltip: {click, hover, focus}
    tooltip_placement : str ("bottom")
        Position of the tooltip: {auto, top, bottom, left, right}
    n_cols : int (5)
        Number of columns per page
    n_rows` : int (3)
        Number of rows per page (only available for the "pages" template)
    border : str ("1px solid #cccccc")
        Styling of the border around each cell (CSS)
    gap : int or str (0)
        Size of the margin around each cell (CSS)
    fontsize : str ("12pt")
        Font size of the text displayed in each cell (CSS)
    fontfamily : str ("'DejaVu', sans-serif")
        Font used for the text in each cell (CSS)
    textalign : str ("center")
        Alignment of the text in each cell (CSS)
    hover_color : str ("#e7e7e7")
        Background color when hovering a cell (CSS)
    style : dict (None)
        CSS styling applied to each item in a cell. The dict must follow a
        `key: function` structure where the key must correspond to one of the
        columns in `subset` or `tooltip`. The function takes the item's value
        as input, and outputs a valid CSS styling. For example, if you want to
        color the text corresponding to the "Solubility" column in your
        dataframe:
        ```python
        style={"Solubility": lambda x: "color: red" if x < -5 else "color: black"}
        ```
    
    Notes
    -----
    You can also directly use RDKit's MolDrawOptions parameters as arguments.
    For more info, see https://www.rdkit.org/docs/source/rdkit.Chem.Draw.rdMolDraw2D.html#rdkit.Chem.Draw.rdMolDraw2D.MolDrawOptions
    """
    raise TypeError(f"No display method registered for type {type(arg)!r}")

@display.register(DataFrame)
@display.register(dict)
def _(df, **kwargs):
    template, kwargs, render_kwargs = _prepare_kwargs(kwargs, "display")
    return MolGrid(df, **kwargs).display(template=template, **render_kwargs)

@display.register(str)
def _(sdf, **kwargs):
    template, kwargs, render_kwargs = _prepare_kwargs(kwargs, "display")
    return MolGrid.from_sdf(sdf, **kwargs).display(template=template,
                                                   **render_kwargs)

@display.register(Series)
@display.register(list)
def _(mols, **kwargs):
    template, kwargs, render_kwargs = _prepare_kwargs(kwargs, "display")
    return MolGrid.from_mols(mols, **kwargs).display(template=template,
                                                     **render_kwargs)

@singledispatch
def save(arg, **kwargs):
    """Generate an interactive grid of molecules and save it
    
    Parameters
    ----------
    output : str
        Name and path of the output document
    
    See `mols2grid.display` for the full list of arguments
    """
    raise TypeError(f"No display method registered for type {type(arg)!r}")

@save.register(DataFrame)
def _(df, **kwargs):
    template, kwargs, render_kwargs = _prepare_kwargs(kwargs, "save")
    return MolGrid(df, **kwargs).save(template=template, **render_kwargs)

@save.register(str)
def _(sdf, **kwargs):
    template, kwargs, render_kwargs = _prepare_kwargs(kwargs, "save")
    return MolGrid.from_sdf(sdf, **kwargs).save(template=template,
                                                **render_kwargs)

@save.register(Series)
@save.register(list)
def _(mols, **kwargs):
    template, kwargs, render_kwargs = _prepare_kwargs(kwargs, "save")
    return MolGrid.from_mols(mols, **kwargs).save(template=template,
                                                  **render_kwargs)