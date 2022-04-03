from functools import singledispatch
import inspect
from pandas import DataFrame, Series
from .molgrid import MolGrid


_SIGNATURE = {
    method: dict(inspect.signature(getattr(MolGrid, method)).parameters.items())
    for method in ["render", "to_pages", "to_table", "display"]
}
for method in ["render", "to_pages", "to_table", "display"]:
    _SIGNATURE[method].pop("self")
    if method in ["render", "display"]:
        _SIGNATURE[method].pop("kwargs")


def _prepare_kwargs(kwargs, kind):
    """Separate kwargs for the init and render methods of MolGrid"""
    template = kwargs.pop("template", _SIGNATURE["render"]["template"].default)
    render_kwargs = {param: kwargs.pop(param, sig.default)
                     for param, sig in _SIGNATURE[f'to_{template}'].items()}
    if kind == "display":
        render_kwargs.update({param: kwargs.pop(param, sig.default)
                              for param, sig in _SIGNATURE["display"].items()})
    return template, kwargs, render_kwargs

@singledispatch
def display(arg, **kwargs):
    """Display molecules on an interactive grid

    Parameters
    ----------
    arg : pandas.DataFrame, SDF file or list of molecules
        The input containing your molecules
    smiles_col : str or None, default="SMILES"
        If a pandas DataFrame is used, name of the column with SMILES
    mol_col : str or None, default=None
        If a pandas DataFrame is used, name of the column with RDKit molecules
    useSVG : bool, default=True
        Use SVG images or PNG
    coordGen : bool, default=True
        Use the coordGen library instead of the RDKit one to depict the
        molecules in 2D
    use_coords : bool, default=False
        Use the coordinates of the molecules (only relevant when an SDF file, a
        list of molecules or a DataFrame of RDKit molecules were used as input)
    removeHs : bool, default=False
        Remove hydrogen atoms from the drawings
    size : tuple, default=(160, 120)
        Size of each image
    prerender : bool, default=False
        Prerender images for the entire dataset, or generate them on-the-fly
        when needed
    substruct_highlight : bool, default=True
        Highlight substructure when using the SMARTS search. Only available
        when ``prerender=False``
    single_highlight : bool, default=False
        Highlight only the first match of the substructure query
    MolDrawOptions : rdkit.Chem.Draw.rdMolDraw2D.MolDrawOptions or None, default=None
        Drawing options. Useful for making highly customized drawings
    rename : dict or None, default=None
        Rename the properties in the final document
    name : str, default="default"
        Name of the grid. Used when retrieving selections from multiple grids
        at the same time
    template : str, default="pages"
        Either ``"pages"`` or ``"table"``
    width : str, default="100%
        Width of the frame displayed in the notebook
    height : int or None, default=None
        Height of the frame displayed in the notebook. Use ``None`` for an
        automatic guess
    subset: list or None, default=None
        Columns to be displayed in each cell of the grid. Each column's value
        will be displayed from top to bottom in the same order given here.
        Use ``"img"`` for the image of the molecule, and ``"mols2grid-id"`` for
        the molecule's index in your input file.
    tooltip : list, None or False, default=None
        Columns to be displayed as a tooltip when hovering/clicking on the
        image of a cell. Use ``False`` for no tooltip.
    tooltip_fmt : str, default="<strong>{key}</strong>: {value}"
        Format string of each key/value pair in the tooltip
    tooltip_trigger : str, default="click hover"
        Sequence of triggers for the tooltip: ``click``, ``hover`` or ``focus``
    tooltip_placement : str, default="bottom"
        Position of the tooltip: ``auto``, ``top``, ``bottom``, ``left`` or
        ``right``
    n_cols : int, default=5
        Number of columns per page
    n_rows` : int, default=3
        Only available for the "pages" template. Number of rows per page
    border : str, default="1px solid #cccccc"
        Styling of the border around each cell (CSS)
    gap : int, default=0
        Size of the margin around each cell in px
    fontsize : str, default="12pt"
        Font size of the text displayed in each cell (CSS)
    fontfamily : str, default="'DejaVu', sans-serif"
        Font used for the text in each cell (CSS)
    textalign : str, default="center"
        Alignment of the text in each cell (CSS)
    hover_color : str, default="#e7e7e7"
        Background color when hovering a cell (CSS)
    style : dict or None, default=None
        CSS styling applied to each item in a cell. The dict must follow a
        ``key: function`` structure where the key must correspond to one of the
        columns in ``subset`` or ``tooltip``. The function takes the item's
        value as input, and outputs a valid CSS styling. For example, if you
        want to color the text corresponding to the "Solubility" column in your
        dataframe::

            style={"Solubility": lambda x: "color: red" if x < -5 else ""}

        You can also style a whole cell using the ``__all__`` key, the
        corresponding function then has access to all values for each cell::

            style={"__all__": lambda x: "color: red" if x["Solubility"] < -5 else ""}

    selection : bool, default=True
        Only available for the "pages" template. Enables the selection of
        molecules and displays a checkbox at the top of each cell. To access
        your selection (index and SMILES), use :func:`mols2grid.get_selection()`
        or the export options in the bottom checkbox dropdown menu.
    cache_selection : bool, default=False
        Restores the selection from a previous grid with the same name
    transform : dict or None, default=None
        Functions applied to specific items in all cells. The dict must follow
        a ``key: function`` structure where the key must correspond to one of
        the columns in ``subset`` or ``tooltip``. The function takes the item's
        value as input and transforms it, for example::

            transform={"Solubility": lambda x: f"{x:.2f}",
                       "Melting point": lambda x: f"MP: {5/9*(x-32):.1f}°C"}

        These transformations only affect columns in ``subset`` and
        ``tooltip``, and do not interfere with ``style``.
    custom_css : str or None, default=None
        Only available for the "pages" template. Custom CSS properties applied
        to the content of the HTML document.
    custom_header : str or None, default=None
        Custom libraries to be loaded in the header of the document
    callback : str or None, default=None
        Only available for the "pages" template. JavaScript or Python
        callback to be executed when clicking on an image. A dictionnary
        containing the data for the full cell is directly available as
        ``data`` in JS. For Python, the callback function must have
        ``data`` as the first argument to the function. All the values in
        the ``data`` dict are parsed as strings, except "mols2grid-id"
        which is always an integer. Note that fields containing spaces in
        their name will be replaced by hyphens, i.e. "mol weight" becomes
        available as ``data["mol-weight"]``.
    sort_by : str or None, default=None
        Sort the grid according to the following field (which must be present
        in ``subset`` or ``tooltip``).

    Returns
    -------
    view : IPython.core.display.HTML

    Notes
    -----
    You can also directly use RDKit's :class:`~rdkit.Chem.Draw.rdMolDraw2D.MolDrawOptions`
    parameters as arguments.
    Additionally, ``atomColourPalette`` is available to customize the atom
    palette if you're not prerendering image (``prerender=False``).

    If ``subset=None, tooltip=None``, the index and image will be directly
    displayed on the grid while the remaining fields will be in the tooltip.
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
    arg : pandas.DataFrame, SDF file or list of molecules
        The input containing your molecules
    output : str
        Name and path of the output document

    Notes
    -----
    See :func:`display` for the full list of arguments
    """
    raise TypeError(f"No save method registered for type {type(arg)!r}")

@save.register(DataFrame)
def _(df, **kwargs):
    template, kwargs, render_kwargs = _prepare_kwargs(kwargs, "save")
    output = kwargs.pop("output")
    return MolGrid(df, **kwargs).save(output, template=template,
                                      **render_kwargs)

@save.register(str)
def _(sdf, **kwargs):
    template, kwargs, render_kwargs = _prepare_kwargs(kwargs, "save")
    output = kwargs.pop("output")
    return MolGrid.from_sdf(sdf, **kwargs).save(output, template=template,
                                                **render_kwargs)

@save.register(Series)
@save.register(list)
def _(mols, **kwargs):
    template, kwargs, render_kwargs = _prepare_kwargs(kwargs, "save")
    output = kwargs.pop("output")
    return MolGrid.from_mols(mols, **kwargs).save(output, template=template,
                                                  **render_kwargs)