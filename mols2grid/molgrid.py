import json
import warnings
from base64 import b64encode
from functools import partial
from html import escape

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw

from .callbacks import _JSCallback
from .select import register
from .utils import (
    callback_handler,
    env,
    is_running_within_streamlit,
    mol_to_record,
    mol_to_smiles,
    remove_coordinates,
    requires,
    sdf_to_dataframe,
    slugify,
    tooltip_formatter,
)
from .widget import MolGridWidget

try:
    from IPython.display import HTML, Javascript, display
except ModuleNotFoundError:
    pass
else:
    warnings.filterwarnings("ignore", "Consider using IPython.display.IFrame instead.")


# Detect if mols2grid is running inside a Jupyter Notebook/Lab.
# If it is, we wrap the HTML in an iframe.
try:
    get_ipython()  # This is callable only in Jupyter Notebook.
    is_jupyter = True
except NameError:
    is_jupyter = False


class MolGrid:
    """Class that handles drawing molecules, rendering the HTML document and
    saving or displaying it in a Jupyter Notebook.

    Parameters: Data
    ----------------
    df : pandas.DataFrame, dict or list, required
        Dataframe containing a SMILES or mol column, or dictionary containing a
        list of SMILES, or list of dictionnaries containing a SMILES field.
    smiles_col : str or None, default="SMILES"
        Name of the SMILES column in the dataframe, if available.
    mol_col : str or None, default=None
        Name of an RDKit molecule column. If available, coordinates and
        atom/bonds annotations from this will be used for depiction.

    Parameters: Display
    -------------------
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
    cache_selection : bool, default=False
        Restores the selection from a previous grid with the same name.

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

    Parameters: Customization
    -------------------------
    name : str, default="default"
        Name of the grid. Used when retrieving selections from multiple grids
        at the same time.
    rename : dict or None, default=None
        Rename the properties in the final document.
    kwargs : object
        :class:`~rdkit.Chem.Draw.rdMolDraw2D.MolDrawOptions` attributes, and
        the additional ``atomColourPalette``.

    Notes
    -----
    On-the-fly rendering of images does not read the atom colour palette
    from the :class:`~rdkit.Chem.Draw.rdMolDraw2D.MolDrawOptions` parameter.
    If this is needed, use the following::

        MolGrid(df, atomColourPalette={1: (.8, 0, 1)})

    .. versionchanged:: 0.1.0
        Added ``rename`` parameter to replace ``mapping``.

    .. versionadded:: 0.2.0
        Added ``prerender`` and ``cache_selection`` parameters.

    .. versionchanged:: 0.2.0
        Images are now generated on-the-fly. ``use_coords`` is now ``False`` by
        default to avoid a systematic error when using ``MolGrid.from_sdf``.
    """

    def __init__(
        self,
        df,
        smiles_col="SMILES",
        mol_col=None,
        #
        size=(130, 90),
        useSVG=True,
        prerender=False,
        cache_selection=False,
        #
        removeHs=False,
        use_coords=False,
        coordGen=True,
        MolDrawOptions=None,
        #
        name="default",
        rename=None,
        **kwargs,
    ):
        if not (smiles_col or mol_col):
            raise ValueError("One of `smiles_col` or `mol_col` must be set")
        if not isinstance(name, str):
            raise TypeError(
                f"`name` must be a string. Currently of type {type(name).__name__}"
            )
        if not prerender:
            if not useSVG:
                raise ValueError("On-the-fly rendering of PNG images not supported")
            if use_coords and mol_col:
                raise ValueError("Cannot use coordinates with on-the-fly rendering")
        self.prefer_coordGen = coordGen
        self.removeHs = removeHs
        self.useSVG = useSVG
        self.use_coords = use_coords
        self.img_size = size
        self.prerender = prerender
        self.smiles_col = smiles_col
        self.mol_col = mol_col
        if isinstance(df, pd.DataFrame):
            dataframe = df.copy()
        else:
            # List of dicts or other input formats for dataframes.
            dataframe = pd.DataFrame(df)
        if rename:
            dataframe.rename(columns=rename, inplace=True)
        self._extra_columns = ["img", "mols2grid-id"]
        # Add index.
        dataframe["mols2grid-id"] = list(range(len(dataframe)))
        # Generate drawing options.
        if prerender:
            Draw.rdDepictor.SetPreferCoordGen(coordGen)
            opts = MolDrawOptions or Draw.MolDrawOptions()
            for key, value in kwargs.items():
                setattr(opts, key, value)
            self.MolDrawOptions = opts
            self._MolDraw2D = Draw.MolDraw2DSVG if useSVG else Draw.MolDraw2DCairo
        else:
            opts = {}
            if MolDrawOptions:
                for key in dir(MolDrawOptions):
                    value = getattr(MolDrawOptions, key)
                    if not (
                        key.startswith("_")
                        or callable(value)
                        or value.__class__.__module__ != "builtins"
                    ):
                        opts[key] = value
            opts.update(kwargs)
            opts.update({"width": self.img_size[0], "height": self.img_size[1]})
            self.json_draw_opts = json.dumps(opts)
        # Prepare smiles and images.
        self._prepare_dataframe(dataframe)
        self.dataframe = dataframe
        # Register instance.
        self._grid_id = name
        if cache_selection:
            try:
                self._cached_selection = register.get_selection(name)
            except KeyError:
                self._cached_selection = {}
                register._init_grid(name)
            else:
                register._update_current_grid(name)
        else:
            self._cached_selection = {}
            register._init_grid(name)
        # Create widget.
        widget = MolGridWidget(grid_id=name, selection=str(self._cached_selection))
        selection_handler = partial(register.selection_updated, name)
        widget.observe(selection_handler, names=["selection"])
        # Register widget JS-side.
        display(widget)
        self.widget = widget

    @classmethod
    def from_mols(cls, mols, **kwargs):
        """Set up the dataframe used by mols2grid directly from a list of RDKit
        molecules.

        Parameters
        ----------
        mols : list
            List of RDKit molecules
        kwargs : object
            Other arguments passed on initialization
        """
        mol_col = kwargs.pop("mol_col", "mol")
        df = pd.DataFrame([mol_to_record(mol, mol_col=mol_col) for mol in mols])
        return cls(df, mol_col=mol_col, **kwargs)

    @classmethod
    def from_sdf(cls, sdf_file, **kwargs):
        """Set up the dataframe used by mols2grid directly from an SDFile.

        Parameters
        ----------
        sdf_file : str, pathlib.Path
            Path to the SDF file (.sdf or .sdf.gz)
        kwargs : object
            Other arguments passed on initialization


        .. versionchanged:: 0.2.0
            Added support for `.sdf.gz` files
        """
        mol_col = kwargs.pop("mol_col", "mol")
        df = sdf_to_dataframe(sdf_file, mol_col=mol_col)
        return cls(df, mol_col=mol_col, **kwargs)

    @property
    def template(self):
        """Kind of grid displayed, one of:

        * interactive
        * static
        """
        return self._template

    @template.setter
    def template(self, value):
        if value not in ["interactive", "static"]:
            raise ValueError(
                f"template={value!r} not supported. "
                "Use either 'interactive' or 'static'."
            )
        self._template = value

    def draw_mol(self, mol):
        """Draw a molecule."""
        d2d = self._MolDraw2D(*self.img_size)
        d2d.SetDrawOptions(self.MolDrawOptions)
        hl_atoms = getattr(mol, "__sssAtoms", [])
        d2d.DrawMolecule(mol, highlightAtoms=hl_atoms)
        d2d.FinishDrawing()
        return d2d.GetDrawingText()

    def mol_to_img(self, mol):
        """Convert an RDKit mol to an inline PNG image containing a drawing of the
        molecule."""
        img = self.draw_mol(mol)
        if self.useSVG:
            return img
        data = b64encode(img).decode()
        return f'<img src="data:image/png;base64,{data}">'

    def _prepare_dataframe(self, dataframe):
        """Prepares the dataframe with SMILES and images depending on user input. The
        dataframe is modified inplace."""
        if self.prerender:
            if self.mol_col:
                keep_mols = True
            else:
                # Make temporary mol column if not present.
                self.mol_col = "mol"
                keep_mols = False
                dataframe[self.mol_col] = dataframe[self.smiles_col].apply(
                    Chem.MolFromSmiles
                )
            # Drop empty mols.
            dataframe.dropna(axis=0, subset=[self.mol_col], inplace=True)
            # Modify mol according to user pref.
            if not self.use_coords:
                dataframe[self.mol_col] = dataframe[self.mol_col].apply(
                    remove_coordinates
                )
            if self.removeHs:
                dataframe[self.mol_col] = dataframe[self.mol_col].apply(Chem.RemoveHs)
            # Render.
            dataframe["img"] = dataframe[self.mol_col].apply(self.mol_to_img)
            # Cleanup.
            if not keep_mols:
                dataframe.drop(columns=self.mol_col, inplace=True)
                self.mol_col = None
        else:
            dataframe["img"] = None
        # Generate smiles col if not present or needs to be updated.
        if self.mol_col and self.smiles_col not in dataframe.columns:
            dataframe[self.smiles_col] = dataframe[self.mol_col].apply(mol_to_smiles)

    def render(self, template="interactive", **kwargs):
        """Returns the HTML document corresponding to the "interactive" or "static"
        template. See :meth:`to_interactive` and :meth:`to_static` for the full list
        of arguments.

        Parameters
        ----------
        template : str
            What kind of grid to draw:

            * interactive
                An interactive grid that layouts the original set of
                molecules on several pages, allowing for selecting molecules and
                filtering them using text or substructure queries.
            * static
                A simple table with all molecules displayed at once, similarly
                to RDKit's :func:`~rdkit.Chem.Draw.rdMolDraw2D.MolsToGridImage`.
                This template is mainly used for printing on paper or in a PDF
                file. Most of the interactive actions aren't available.
        """
        self.template = template
        return getattr(self, f"to_{self.template}")(**kwargs)

    def to_interactive(
        self,
        # Display
        subset=None,
        tooltip=None,
        tooltip_fmt="<strong>{key}</strong>: {value}",
        tooltip_placement="auto",
        transform=None,
        sort_by=None,
        use_iframe=False,
        truncate=True,
        n_items_per_page=24,
        selection=True,
        # Mols
        substruct_highlight=None,
        single_highlight=False,
        # CSS
        border="1px solid #cccccc",
        gap=0,
        pad=10,
        fontsize="12px",
        fontfamily="'DejaVu', sans-serif",
        textalign="center",
        background_color="white",
        hover_color="rgba(0,0,0,0.05)",
        custom_css=None,
        style=None,
        # Customization
        custom_header=None,
        callback=None,
        **kwargs,
    ):
        """Returns the HTML document for the "interactive" template.

        Parameters: Display
        -------------------
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
        tooltip_placement : str, default="auto"
            Position of the tooltip: ``auto``, ``top``, ``bottom``, ``left`` or
            ``right``.
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
        use_iframe : bool, default=False
            Whether to use an iframe to display the grid. When the grid is displayed
            inside a Jupyter Notebook or JupyterLab, this will default to ``True``.
        truncate: bool, default=True/False
            Whether to truncate the text in each cell if it's too long.
            Defaults to ``True`` for interactive grids, ``False`` for static grid.
        n_items_per_page, default=24
            Number of items to display per page. A multiple of 12 is recommended
            for optimal display.
        selection : bool, default=True
            Enables the selection of molecules and displays a checkbox at the
            top of each cell. In the context of a Jupyter Notebook, this gives
            you access to your selection (index and SMILES) through
            :func:`mols2grid.get_selection()` or :meth:`MolGrid.get_selection()`.
            In all cases, you can export your selection by clicking on the triple-dot menu.

        Parameters: Mols
        ----------------
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
            Background color of a cell.
        hover_color : str, default="rgba(0,0,0,0.05)"
            Background color when hovering a cell.
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
            custom_header : str or None, default=None
                Custom libraries to be loaded in the header of the document.
            callback : str, callable or None, default=None
                JavaScript or Python callback to be executed when clicking on an image.
                A dictionnary containing the data for the full cell is directly available
                as ``data`` in JS. For Python, the callback function must have ``data``
                as the first argument to the function. All the values in the ``data`` dict
                are parsed as strings, except "mols2grid-id" which is always an integer.
                Note that fields containing spaces in their name will be replaced by
                hyphens, i.e. "mol weight" becomes available as ``data["mol-weight"]``.

        Returns
        -------
        html_document : str

        Notes
        -----
        If ``subset=None, tooltip=None``, the index and image will be directly
        displayed on the grid while the remaining fields will be in the
        tooltip.

        The cell width is defined by the size[0] parameter.

        .. versionadded:: 0.1.0
            Added ``sort_by``, ``custom_css``, ``custom_header`` and
            ``callback`` arguments.
            Added the ability to style an entire cell with
            ``style={"__all__": <function>}``.

        .. versionadded:: 0.2.0
            Added ``substruct_highlight`` argument.

        .. versionchanged:: 0.2.2
            If both ``subset`` and ``tooltip`` are ``None``, the index and
            image will be directly displayed on the grid while the remaining
            fields will be in the tooltip.

        .. versionchanged:: 1.0.0
            ``callback`` can now be a *lambda* function. If ``prerender=True``,
            substructure highlighting will be automatically disabled if it
            wasn't explicitely set to ``True`` instead of raising an error.

        """
        if substruct_highlight is None:
            substruct_highlight = not self.prerender
        if substruct_highlight and self.prerender:
            raise ValueError(
                "Cannot highlight substructure search with prerendered images"
            )
        if self.mol_col:
            df = self.dataframe.drop(columns=self.mol_col).copy()
        else:
            df = self.dataframe.copy()
        smiles = self.smiles_col
        content = []  # Gets filled with the HTML content of each cell.
        column_map = {}

        if subset is None:
            if tooltip is None:
                subset = ["mols2grid-id", "img"]
                tooltip = [x for x in df.columns.tolist() if x not in subset]
            else:
                # When no subset is defined, all columns are displayed.
                subset = df.columns.tolist()
        else:
            # work on a copy
            subset = subset[:]

        if "mols2grid-id" not in subset:
            subset.insert(0, "mols2grid-id")
        if "img" not in subset:
            subset.insert(0, "img")

        # Always make sure the image comes first.
        # subset = [subset.pop(subset.index("img"))] + subset
        #
        # This was removed at Cedric's request, so you can choose
        # to have certain properties displayed above the image.

        # Define fields that are searchable and sortable.
        search_cols = [f"data-{col}" for col in subset if col != "img"]
        if tooltip:
            search_cols.extend([f"data-{col}" for col in tooltip])
            for col in tooltip:
                if col not in subset:
                    s = f'<div class="data data-{slugify(col)}" style="display: none;"></div>'
                    content.append(s)
                    column_map[col] = f"data-{col}"
        else:
            tooltip = []
        sort_cols = search_cols[:]
        sort_cols = ["data-mols2grid-id"] + sort_cols

        # Get unique list but keep order.
        sort_cols = list(dict.fromkeys(sort_cols))
        if style is None:
            style = {}
        if transform is None:
            transform = {}
        value_names = list(set(subset + [smiles] + tooltip))
        value_names = [f"data-{col}" for col in value_names]

        # Force id, SMILES, and tooltip values to be present in the data.
        final_columns = subset[:]
        final_columns.extend(["mols2grid-id", smiles])
        if tooltip:
            final_columns.extend(tooltip)
        final_columns = list(set(final_columns))

        # Make a copy of id shown explicitly.
        id_name = "mols2grid-id-display"
        df[id_name] = df["mols2grid-id"]
        value_names.append(f"data-{id_name}")
        final_columns.append(id_name)
        subset = [id_name if x == "mols2grid-id" else x for x in subset]
        id_display_html = f'<div class="data-{id_name}"></div>'

        # Organize data.
        temp = []
        for col in subset:
            if col == "mols2grid-id-display":
                s = ""  # Avoid an empty div to be created for the display id.
            elif col == "img" and tooltip:
                s = f'<a class="data data-{col}"></a>'
            else:
                if style.get(col):
                    s = (
                        f'<div class="data data-{slugify(col)} copy-me style-{slugify(col)}" '
                        'style=""></div>'
                    )
                else:
                    s = f'<div class="data data-{slugify(col)} copy-me"></div>'
            temp.append(s)
            column_map[col] = f"data-{col}"
        content = temp + content

        # Add but hide SMILES div if not present.
        if smiles not in (subset + tooltip):
            s = f'<div class="data data-{slugify(smiles)} copy-me" style="display: none;"></div>'
            content.append(s)
            column_map[smiles] = f"data-{smiles}"

        # Set mapping for list.js.
        if "__all__" in style.keys():
            whole_cell_style = True
            x = "[{data: ['mols2grid-id', 'cellstyle']}, "
        else:
            whole_cell_style = False
            x = "[{data: ['mols2grid-id']}, "
        value_names = [slugify(c) for c in value_names]
        value_names = x + str(value_names)[1:]

        # Apply CSS styles.
        for col, func in style.items():
            if col == "__all__":
                name = "cellstyle"
                df[name] = df.apply(func, axis=1)
            else:
                name = f"style-{slugify(col)}"
                df[name] = df[col].apply(func)
            final_columns.append(name)
            value_names = value_names[:-1] + f", {{ attr: 'style', name: {name!r} }}]"

        # Create tooltip.
        if tooltip:
            df["m2g-tooltip"] = df.apply(
                tooltip_formatter, axis=1, args=(tooltip, tooltip_fmt, style, transform)
            )
            final_columns += ["m2g-tooltip"]
            value_names = (
                value_names[:-1] + ", {attr: 'data-content', name: 'm2g-tooltip'}]"
            )
            info_btn_html = '<div class="m2g-info">i</div>'
        else:
            info_btn_html = ""

        # Apply custom user function.
        for col, func in transform.items():
            df[col] = df[col].apply(func)

        # Add checkboxes.
        if selection:
            if self._cached_selection:
                df["cached_checkbox"] = False
                df.loc[
                    df["mols2grid-id"].isin(self._cached_selection.keys()),
                    "cached_checkbox",
                ] = True
                final_columns += ["cached_checkbox"]
                value_names = (
                    value_names[:-1] + ", {attr: 'checked', name: 'cached_checkbox'}]"
                )
            checkbox_html = (
                '<input type="checkbox" tabindex="-1" '
                'class="position-relative float-left cached_checkbox">'
            )
        else:
            checkbox_html = ""

        # Add callback button.
        if callback:
            callback_btn_html = '<div class="m2g-callback"></div>'
        else:
            callback_btn_html = ""

        # Generate cell HTML.
        item = (
            '<div class="m2g-cell" data-mols2grid-id="0" tabindex="0">'
            '<div class="m2g-cb-wrap">{checkbox_html}<div class="m2g-cb"></div>{id_display_html}</div>'
            '<div class="m2g-cell-actions">{info_btn_html}{callback_btn_html}</div>'
            "{content}"
            "{tooltip_html}"
            "</div>"
        )

        item = item.format(
            checkbox_html=checkbox_html,
            id_display_html=id_display_html,
            info_btn_html=info_btn_html,
            callback_btn_html=callback_btn_html,
            content="".join(content),
            tooltip_html='<div class="m2g-tooltip" data-toggle="popover" data-content="."></div>'
            if tooltip
            else "",
        )

        # Callback
        if isinstance(callback, _JSCallback):
            if custom_header and callback.library_src:
                custom_header = callback.library_src + custom_header
            else:
                custom_header = callback.library_src
            callback = callback.code
        if callable(callback):
            callback_type = "python"
            cb_handler = partial(callback_handler, callback)
            self.widget.observe(cb_handler, names=["callback_kwargs"])
        else:
            callback_type = "js"

        # Sort
        if sort_by and sort_by != "mols2grid-id":
            if sort_by in (subset + tooltip):
                sort_by = f"data-{slugify(sort_by)}"
            else:
                raise ValueError(
                    f"{sort_by!r} is not an available field in " "`subset` or `tooltip`"
                )
        else:
            sort_by = "mols2grid-id"

        # Slugify remaining vars.
        column_map = {k: slugify(v) for k, v in column_map.items()}
        sort_cols = [slugify(c) for c in sort_cols]
        search_cols = [slugify(c) for c in search_cols]
        smiles = slugify(smiles)
        df = df[final_columns].rename(columns=column_map).sort_values(sort_by)

        template = env.get_template("interactive.html")
        template_kwargs = dict(
            tooltip=tooltip,
            tooltip_placement=repr(tooltip_placement),
            n_items_per_page=n_items_per_page,
            selection=selection,
            truncate=truncate,
            sort_by=sort_by,
            use_iframe=use_iframe,
            #
            border=border,
            gap=gap,
            gap_px="-1px -1px 0 0" if gap == 0 else f"{gap}px",
            pad=pad,
            fontsize=fontsize,
            fontfamily=fontfamily,
            textalign=textalign,
            background_color=background_color,
            hover_color=hover_color,
            #
            iframe_padding=18,
            cell_width=self.img_size[0],
            image_width=self.img_size[0],
            image_height=self.img_size[1],
            #
            item=item,
            item_repr=repr(item),
            value_names=value_names,
            search_cols=search_cols,
            data=json.dumps(
                df.to_dict("records"), indent=None, default=lambda x: "🤷‍♂️"
            ),
            cached_selection=(
                [
                    list(self._cached_selection.keys()),
                    list(self._cached_selection.values()),
                ]
                if self._cached_selection
                else False
            ),
            smiles_col=smiles,
            sort_cols=sort_cols,
            grid_id=self._grid_id,
            whole_cell_style=whole_cell_style,
            custom_css=custom_css or "",
            custom_header=custom_header or "",
            callback=callback,
            callback_type=callback_type,
            removeHs=self.removeHs,
            prefer_coordGen=self.prefer_coordGen,
            onthefly=not self.prerender,
            substruct_highlight=substruct_highlight,
            json_draw_opts=getattr(self, "json_draw_opts", ""),
            single_highlight=single_highlight,
        )
        return template.render(**template_kwargs)

    def get_selection(self):
        """Retrieve the dataframe subset corresponding to your selection.

        Returns
        -------
        pandas.DataFrame
        """
        sel = list(register.get_selection(self._grid_id).keys())
        return self.dataframe.loc[self.dataframe["mols2grid-id"].isin(sel)].drop(
            columns=self._extra_columns
        )

    def filter(self, mask):
        """Filters the grid using a mask (boolean array).

        Parameters
        ----------
        mask : list, pandas.Series or numpy.ndarray
            Boolean array: ``True`` when the item should be displayed,
            ``False`` if it should be filtered out.
        """
        if isinstance(mask, (pd.Series, np.ndarray)):
            mask = mask.tolist()
        if is_running_within_streamlit():
            filtering_script = env.get_template("js/filter.js").render(
                grid_id=self._grid_id, mask=json.dumps(mask)
            )
            return Javascript(filtering_script)
        else:
            self.widget.filter_mask = mask

    def filter_by_index(self, indices):
        """Filters the grid using the dataframe's index."""

        # convert index to mask
        mask = self.dataframe.index.isin(indices)
        return self.filter(mask)

    def to_static(
        self,
        # Display
        subset=None,
        tooltip=None,
        tooltip_fmt="<strong>{key}</strong>: {value}",
        tooltip_trigger="focus",
        tooltip_placement="auto",
        transform=None,
        sort_by=None,
        use_iframe=False,
        truncate=False,
        n_cols=5,
        # CSS Styling
        border="1px solid #cccccc",
        gap=0,
        pad=10,
        fontsize="12px",
        fontfamily="'DejaVu', sans-serif",
        textalign="center",
        custom_css=None,
        style=None,
        # Customization
        custom_header=None,
        **kwargs,
    ):
        """Returns the HTML document for the "static" template

        Parameters: Display
        -------------------
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
            Sequence of triggers for the tooltip: ``click``, ``hover`` or ``focus``
        tooltip_placement : str, default="auto"
            Position of the tooltip: ``auto``, ``top``, ``bottom``, ``left`` or
            ``right``.
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
        use_iframe : bool, default=False
            Whether to use an iframe to display the grid. When the grid is displayed
            inside a Jupyter Notebook or JupyterLab, this will default to ``True``.
        truncate: bool, default=False
            Whether to truncate the text in each cell if it's too long.
        n_cols : int, default=5
            Number of columns in the table.

        Parameters: CSS
        --------------
        border : str, default="1px solid #cccccc"
            Styling of the border around each cell.
        gap : int, default=0
            Size in pixels of the gap between cells.
        pad: int, default=10
            Size in pixels of the cell padding.
        fontsize : str, default="12pt"
            Font size of the text displayed in each cell.
        fontfamily : str, default="'DejaVu', sans-serif"
            Font used for the text in each cell.
        textalign : str, default="center"
            Alignment of the text in each cell.
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
        -----------------------
        custom_header : str or None
            Custom libraries to be loaded in the header of the document

        Returns
        -------
        html_document : str

        Notes
        -----
        If ``subset=None, tooltip=None``, the index and image will be directly
        displayed on the grid while the remaining fields will be in the
        tooltip.

        .. versionadded:: 0.1.0
            Added the ability to style an entire cell with
            ``style={"__all__": <function>}``

        .. versionadded:: 0.2.2
            Added ``sort_by``, ``custom_css``, ``custom_header`` arguments.

        .. versionchanged:: 0.2.2
            If both ``subset`` and ``tooltip`` are ``None``, the index and
            image will be directly displayed on the grid while the remaining
            fields will be in the tooltip.
        """
        if not self.prerender:
            raise ValueError(
                "Please set `prerender=True` when using the 'static' template"
            )
        tr = []
        data = []
        sort_by = sort_by or "mols2grid-id"
        df = self.dataframe.sort_values(sort_by).reset_index(drop=True)

        if subset is None:
            if tooltip is None:
                subset = ["mols2grid-id", "img"]
                tooltip = [x for x in df.columns.tolist() if x not in subset]
            else:
                # When no subset is defined, all columns are displayed.
                subset = df.columns.tolist()
        else:
            # work on a copy
            subset = subset[:]

        if "mols2grid-id" not in subset:
            subset.insert(0, "mols2grid-id")
        if "img" not in subset:
            subset.insert(0, "img")

        # Always make surer the image comes first.
        subset = [subset.pop(subset.index("img"))] + subset

        if style is None:
            style = {}
        if transform is None:
            transform = {}
        if not tooltip:
            tooltip = []

        for i, row in df.iterrows():
            ncell = i + 1
            nrow, ncol = divmod(i, n_cols)
            popover = tooltip_formatter(row, tooltip, tooltip_fmt, style, transform)
            td = [
                f'<td class="col-{ncol} m2g-tooltip" tabindex="0" data-toggle="popover" data-content="{escape(popover)}">'
            ]
            if "__all__" in style.keys():
                s = style["__all__"](row)
                div = [f'<div class="m2g-cell-{i}" style="{s}">']
            else:
                div = [f'<div class="m2g-cell-{i}">']
            for col in subset:
                v = row[col]
                if col == "img" and tooltip:
                    item = f'<div class="data data-img">' f"{v}</div>"
                else:
                    func = style.get(col)
                    slug_col = slugify(col)
                    if func:
                        item = f'<div class="data copy-me data-{slug_col}" style="{func(v)}">'
                    else:
                        item = f'<div class="data copy-me data-{slug_col}">'
                    func = transform.get(col)
                    v = func(v) if func else v
                    item += f"{v}</div>"
                div.append(item)
            div.append("</div>")
            td.append("\n".join(div))
            td.append("</td>")
            tr.append("\n".join(td))

            if (ncell % n_cols == 0) or (ncell == len(df)):
                cell = [f'<tr class="row-{nrow}">']
                cell.append("\n".join(tr))
                cell.append("</tr>")
                data.append("\n".join(cell))
                tr = []

        template = env.get_template("static.html")
        template_kwargs = dict(
            tooltip=tooltip,
            tooltip_trigger=repr(tooltip_trigger),
            tooltip_placement=repr(tooltip_placement),
            use_iframe=use_iframe,
            truncate=truncate,
            #
            border=border,
            gap=gap,
            pad=pad,
            textalign=textalign,
            fontsize=fontsize,
            fontfamily=fontfamily,
            #
            iframe_padding=18,
            cell_width=self.img_size[0],
            #
            custom_css=custom_css or "",
            custom_header=custom_header or "",
            data="\n".join(data),
        )
        return template.render(**template_kwargs)

    @requires("IPython.display")
    def display(
        self,
        use_iframe=False,
        iframe_width="100%",
        iframe_height=None,
        iframe_allow="clipboard-write",
        iframe_sandbox="allow-scripts allow-same-origin allow-downloads allow-popups allow-modals",
        **kwargs,
    ):
        """Render and display the grid in a Jupyter notebook.

        Returns
        -------
        view : IPython.core.display.HTML
        """
        use_iframe = is_jupyter or use_iframe
        doc = self.render(**kwargs, use_iframe=use_iframe)
        if use_iframe:
            # Render HTML in iframe.
            iframe = env.get_template("html/iframe.html").render(
                width=iframe_width,
                height=iframe_height,
                allow=iframe_allow,
                sandbox=iframe_sandbox,
                doc=escape(doc),
            )
            return HTML(iframe)
        else:
            # Render HTML regularly.
            return HTML(doc)

    def save(self, output, **kwargs):
        """Render and save the grid in an HTML document."""
        with open(output, "w", encoding="utf-8") as f:
            f.write(self.render(**kwargs))
