import inspect
from functools import singledispatch
from pathlib import Path

from pandas import DataFrame, Series

from .molgrid import MolGrid

_SIGNATURE = {
    method: dict(inspect.signature(getattr(MolGrid, method)).parameters.items())
    for method in ["render", "to_interactive", "to_static", "display"]
}
for method in ["render", "to_interactive", "to_static", "display"]:
    _SIGNATURE[method].pop("self")
    if method in ["render", "display"]:
        _SIGNATURE[method].pop("kwargs")


def _prepare_kwargs(kwargs, kind):
    """Separate kwargs for the init and render methods of MolGrid"""
    template = kwargs.pop("template", _SIGNATURE["render"]["template"].default)
    render_kwargs = {
        param: kwargs.pop(param, sig.default)
        for param, sig in _SIGNATURE[f"to_{template}"].items()
    }
    if kind == "display":
        render_kwargs.update(
            {
                param: kwargs.pop(param, sig.default)
                for param, sig in _SIGNATURE["display"].items()
            }
        )
    return template, kwargs, render_kwargs


@singledispatch
def display(arg, **kwargs):
    """Display molecules on an interactive grid.

    Parameters: Data
    ----------------
    arg : pandas.DataFrame, SDF file or list of molecules
        The input containing your molecules.
    smiles_col : str or None, default="SMILES"
        If a pandas dataframe is used, name of the column with SMILES.
    mol_col : str or None, default=None
        If a pandas dataframe is used, name of the column with RDKit molecules.
        If available, coordinates and atom/bonds annotations from this will be
        used for depiction.

    Parameters: Display
    -------------------
    template : str, default="interactive"
        Either ``"interactive"`` or ``"static"``. See ``render()`` for more details.
    size : tuple, default=(130, 90)
        The size of the drawing canvas. The cell minimum width is set to the
        width of the image, so if the cell padding is increased, the image will
        be displayed smaller.
    useSVG : bool, default=True
        Use SVG images instead of PNG.
    prerender : bool, default=False
        Prerender images for the entire dataset, or generate them on-the-fly.
        Prerendering is slow and memory-hungry, but required when ``template="static"``
        or ``useSVG=False``.
    subset: list or None, default=None
        Columns to be displayed in each cell of the grid. Each column's
        value will be displayed from top to bottom in the order provided.
        The ``"img"`` and ``"mols2grid-id"`` columns are displayed by default,
        however you can still add the ``"img"`` column if you wish to change
        the display order.
    tooltip : list, None or False, default=None
        Columns to be displayed inside the tooltip. When no subset is set,
        all columns will be listed in the tooltip by default. Use ``False``
        to hide the tooltip.
    tooltip_fmt : str, default="<strong>{key}</strong>: {value}"
        Format string of each key/value pair in the tooltip.
    tooltip_trigger : str, default="focus"
        Only available for the "static" template.
        Sequence of triggers for the tooltip: ``click``, ``hover`` or ``focus``
    tooltip_placement : str, default="auto"
        Position of the tooltip: ``auto``, ``top``, ``bottom``, ``left`` or
        ``right``
    transform : dict or None, default=None
        Functions applied to specific items in all cells. The dict must follow
        a ``key: function`` structure where the key must correspond to one of
        the columns in ``subset`` or ``tooltip``. The function takes the item's
        value as input and transforms it, for example::

            transform={
                "Solubility": lambda x: f"{x:.2f}",
                "Melting point": lambda x: f"MP: {5/9*(x-32):.1f}°C"
            }

        These transformations only affect columns in ``subset`` and
        ``tooltip``, and do not interfere with ``style``.
    sort_by : str or None, default=None
        Sort the grid according to the following field (which must be
        present in ``subset`` or ``tooltip``).
    truncate: bool, default=True/False
        Whether to truncate the text in each cell if it's too long.
        Defaults to ``True`` for interactive grids, ``False`` for static grid.
    n_items_per_page, default=24
        Only available for the "interactive" template.
        Number of items to display per page. A multiple of 12 is recommended
        for optimal display.
    n_cols : int, default=5
        Only available for the "static" template.
        Number of columns in the table.
    selection : bool, default=True
        Only available for the "interactive" template.
        Enables the selection of molecules and displays a checkbox at the
        top of each cell. In the context of a Jupyter Notebook, this gives
        you access to your selection (index and SMILES) through
        :func:`mols2grid.get_selection()` or :meth:`MolGrid.get_selection()`.
        In all cases, you can export your selection by clicking on the triple-dot menu.
    cache_selection : bool, default=False
        Only available for the "interactive" template.
        Restores the selection from a previous grid with the same name.
    use_iframe : bool, default=False
        Whether to use an iframe to display the grid. When the grid is displayed
        inside a Jupyter Notebook or JupyterLab, this will default to ``True``.
    iframe_width : str, default="100%
        Width of the iframe
    iframe_height : int or None, default=None
        Height of the frame. When set to ``None``, the height is set dynamically
        based on the content.

    Parameters: Mols
    ----------------
    removeHs : bool, default=False
        Remove hydrogen atoms from the drawings.
    use_coords : bool, default=False
        Use the coordinates of the molecules (only relevant when an SDF file, a
        list of molecules or a DataFrame of RDKit molecules were used as input.)
    coordGen : bool, default=True
        Use the CoordGen library instead of the RDKit one to depict the
        molecules in 2D.
    MolDrawOptions : rdkit.Chem.Draw.rdMolDraw2D.MolDrawOptions or None, default=None
        Drawing options. Useful for making highly customized drawings.
    substruct_highlight : bool or None, default=None
        Highlight substructure when using the SMARTS search. Active by default
        when ``prerender=False``.
    single_highlight : bool, default=False
        Highlight only the first match of the substructure query.

    Parameters: CSS
    ---------------
    border : str, default="1px solid #cccccc"
        Styling of the border around each cell.
    gap : int, default=0
        Size in pixels of the gap between cells.
    pad : int, default=10
        Size in pixels of the cell padding.
    fontsize : str, default="12px"
        Font size of the text displayed in each cell.
    fontfamily : str, default="'DejaVu', sans-serif"
        Font used for the text in each cell.
    textalign : str, default="center"
        Alignment of the text in each cell.
    background_color : str, default="white"
        Only available for the "interactive" template.
        Background color of a cell.
    hover_color : str, default="rgba(0,0,0,0.05)"
        Only available for the "interactive" template.
        Background color when hovering a cell
    custom_css : str or None, default=None
        Custom CSS properties applied to the generated HTML. Please note that
        the CSS will apply to the entire page if no iframe is used (see
        ``use_iframe`` for more details).
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

    Parameters: Customization
    -------------------------
    name : str, default="default"
        Name of the grid. Used when retrieving selections from multiple grids
        at the same time
    rename : dict or None, default=None
        Rename the properties in the final document.
    custom_header : str or None, default=None
        Custom libraries to be loaded in the header of the document.
    callback : str, callable or None, default=None
        Only available for the "interactive" template.
        JavaScript or Python callback to be executed when clicking on an image.
        A dictionnary containing the data for the full cell is directly available
        as ``data`` in JS. For Python, the callback function must have ``data``
        as the first argument to the function. All the values in the ``data`` dict
        are parsed as strings, except "mols2grid-id" which is always an integer.
        Note that fields containing spaces in their name will be replaced by
        hyphens, i.e. "mol weight" becomes available as ``data["mol-weight"]``.

    Returns
    -------
    view : IPython.core.display.HTML

    Notes
    -----
    You can also directly use RDKit's :class:`~rdkit.Chem.Draw.rdMolDraw2D.MolDrawOptions`
    parameters as arguments.
    Additionally, ``atomColourPalette`` is available to customize the atom
    palette if you're not prerendering image (``prerender=False``).

    .. versionadded:: 0.1.0
        Added ``sort_by``, ``custom_css``, ``custom_header`` and ``callback``
        arguments. Added the ability to style an entire cell with
        ``style={"__all__": <function>}``.

    .. versionadded:: 0.2.0
        Added ``substruct_highlight`` argument.

    .. versionchanged:: 0.2.2
        If both ``subset`` and ``tooltip`` are ``None``, the index and image
        will be directly displayed on the grid while the remaining fields will
        be in the tooltip.

    .. versionchanged:: 1.0.0
        ``callback`` can now be a *lambda* function. If ``prerender=True``,
        substructure highlighting will be automatically disabled if it wasn't
        explicitely set to ``True`` instead of raising an error.

    """
    raise TypeError(f"No display method registered for type {type(arg)!r}")


@display.register(DataFrame)
@display.register(dict)
def _(df, **kwargs):
    template, kwargs, render_kwargs = _prepare_kwargs(kwargs, "display")
    return MolGrid(df, **kwargs).display(template=template, **render_kwargs)


@display.register(str)
@display.register(Path)
def _(sdf, **kwargs):
    template, kwargs, render_kwargs = _prepare_kwargs(kwargs, "display")
    return MolGrid.from_sdf(sdf, **kwargs).display(template=template, **render_kwargs)


@display.register(Series)
@display.register(list)
@display.register(tuple)
def _(mols, **kwargs):
    template, kwargs, render_kwargs = _prepare_kwargs(kwargs, "display")
    return MolGrid.from_mols(mols, **kwargs).display(template=template, **render_kwargs)


@singledispatch
def save(arg, **kwargs):
    """Generate an interactive grid of molecules and save it.

    Parameters
    ----------
    arg : pandas.DataFrame, SDF file or list of molecules
        The input containing your molecules.
    output : str
        Name and path of the output document.

    Notes
    -----
    See :func:`display` for the full list of arguments.
    """
    raise TypeError(f"No save method registered for type {type(arg)!r}")


@save.register(DataFrame)
def _(df, **kwargs):
    template, kwargs, render_kwargs = _prepare_kwargs(kwargs, "save")
    output = kwargs.pop("output")
    return MolGrid(df, **kwargs).save(output, template=template, **render_kwargs)


@save.register(str)
@save.register(Path)
def _(sdf, **kwargs):
    template, kwargs, render_kwargs = _prepare_kwargs(kwargs, "save")
    output = kwargs.pop("output")
    return MolGrid.from_sdf(sdf, **kwargs).save(
        output, template=template, **render_kwargs
    )


@save.register(Series)
@save.register(list)
@save.register(tuple)
def _(mols, **kwargs):
    template, kwargs, render_kwargs = _prepare_kwargs(kwargs, "save")
    output = kwargs.pop("output")
    return MolGrid.from_mols(mols, **kwargs).save(
        output, template=template, **render_kwargs
    )
