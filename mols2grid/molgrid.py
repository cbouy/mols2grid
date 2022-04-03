import warnings
from base64 import b64encode
from html import escape
import json
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Draw
from .utils import (env,
                    requires,
                    tooltip_formatter,
                    mol_to_record,
                    mol_to_smiles,
                    sdf_to_dataframe,
                    remove_coordinates,
                    slugify)
from .select import register
try:
    from IPython.display import HTML, Javascript
except ModuleNotFoundError:
    pass
else:
    warnings.filterwarnings("ignore",
                            "Consider using IPython.display.IFrame instead")


class MolGrid:
    """Class that handles drawing molecules, rendering the HTML document and
    saving or displaying it in a notebook

    Parameters
    ----------
    df : pandas.DataFrame, dict or list
        Dataframe containing a SMILES or mol column, or dictionary containing a
        list of SMILES, or list of dictionnaries containing a SMILES field
    smiles_col : str or None
        Name of the SMILES column in the dataframe, if available
    mol_col : str or None
        Name of an RDKit molecule column. If available, coordinates and
        atom/bonds annotations from this will be used for depiction
    removeHs : bool
        Remove hydrogen atoms from the drawings
    use_coords : bool
        Use the existing coordinates of the molecule
    coordGen : bool
        Sets whether or not the CoordGen library should be preferred to the
        RDKit depiction library
    useSVG : bool
        Use SVG instead of PNG
    size : tuple
        The size of the drawing canvas
    MolDrawOptions : rdkit.Chem.Draw.rdMolDraw2D.MolDrawOptions or None
        Drawing options. Useful for making highly customized drawings
    rename : dict or None
        Rename the properties/fields stored in the molecule
    name : str
        Name of the grid. Used when retrieving selections from multiple grids
        at the same time
    cache_selection : bool
        Restores the selection from a previous grid with the same name
    prerender : bool
        Prerender images for the entire dataset, or generate them on-the-fly
        when needed
    kwargs : object
        :class:`~rdkit.Chem.Draw.rdMolDraw2D.MolDrawOptions` attributes, and
        the additional ``atomColourPalette``

    Notes
    -----
    On-the-fly rendering of images does not read the atom colour palette
    from the :class:`~rdkit.Chem.Draw.rdMolDraw2D.MolDrawOptions` parameter.
    If this is needed, use the following::

        MolGrid(df, atomColourPalette={1: (.8, 0, 1)})

    .. versionchanged:: 0.1.0
        Added ``rename`` parameter to replace ``mapping``

    .. versionadded:: 0.2.0
        Added ``prerender`` and ``cache_selection`` parameters

    .. versionchanged:: 0.2.0
        Images are now generated on-the-fly. ``use_coords`` is now ``False`` by
        default to avoid a systematic error when using ``MolGrid.from_sdf``
    """

    def __init__(self, df, smiles_col="SMILES", mol_col=None, removeHs=False,
        use_coords=False, coordGen=True, useSVG=True, size=(160, 120),
        MolDrawOptions=None, rename=None, name="default", prerender=False,
        cache_selection=False,
        **kwargs):

        if not (smiles_col or mol_col):
            raise ValueError("One of `smiles_col` or `mol_col` must be set")
        if not isinstance(name, str):
            raise TypeError(
                f"`name` must be a string. Currently of type {type(name).__name__}")
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
            # list of dicts or other input formats for dataframes
            dataframe = pd.DataFrame(df)
        mapping = kwargs.pop("mapping", None)
        if mapping:
            warnings.warn(
                "`mapping` is deprecated and will be removed soon. Consider "
                "using `rename` in the future."
            )
        rename = rename or mapping
        if rename:
            dataframe.rename(columns=rename, inplace=True)
        self._extra_columns = ["img", "mols2grid-id"]
        # add index
        dataframe["mols2grid-id"] = list(range(len(dataframe)))
        # generate drawing options
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
                    if not (key.startswith("_")
                            or callable(value)
                            or value.__class__.__module__ != "builtins"):
                        opts[key] = value
            opts.update(kwargs)
            opts.update({"width": self.img_size[0],
                         "height": self.img_size[1]})
            self.json_draw_opts = json.dumps(opts)
        # prepare smiles and images
        self._prepare_dataframe(dataframe)
        self.dataframe = dataframe
        # register instance
        self._grid_id = name
        if cache_selection:
            try:
                self._cached_selection = register.get_selection(name)
            except KeyError:
                self._cached_selection = None
                register._init_grid(name)
            else:
                register._update_current_grid(name)
        else:
            self._cached_selection = None
            register._init_grid(name)

    @classmethod
    def from_mols(cls, mols, **kwargs):
        """Set up the dataframe used by mols2grid directly from a list of RDKit
        molecules

        Parameters
        ----------
        mols : list
            List of RDKit molecules
        kwargs : object
            Other arguments passed on initialization
        """
        mol_col = kwargs.pop("mol_col", "mol")
        df = pd.DataFrame([mol_to_record(mol, mol_col=mol_col)
                           for mol in mols])
        return cls(df, mol_col=mol_col, **kwargs)

    @classmethod
    def from_sdf(cls, sdf_file, **kwargs):
        """Set up the dataframe used by mols2grid directly from an SDFile

        Parameters
        ----------
        sdf_file : str
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

            * pages
            * table
        """
        return self._template

    @template.setter
    def template(self, value):
        if value not in ["pages", "table"]:
            raise ValueError(f"template={value!r} not supported. "
                             "Use one of 'pages' or 'table'")
        self._template = value

    def draw_mol(self, mol):
        """Draw a molecule"""
        d2d = self._MolDraw2D(*self.img_size)
        d2d.SetDrawOptions(self.MolDrawOptions)
        hl_atoms = getattr(mol, "__sssAtoms", [])
        d2d.DrawMolecule(mol, highlightAtoms=hl_atoms)
        d2d.FinishDrawing()
        return d2d.GetDrawingText()

    def mol_to_img(self, mol):
        """Convert an RDKit mol to an HTML image containing a drawing of the
        molecule"""
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
                # make temporary mol column if not present
                self.mol_col = "mol"
                keep_mols = False
                dataframe[self.mol_col] = (
                    dataframe[self.smiles_col].apply(Chem.MolFromSmiles))
            # drop empty mols
            dataframe.dropna(axis=0, subset=[self.mol_col], inplace=True)
            # modify mol according to user pref
            if not self.use_coords:
                dataframe[self.mol_col] = (
                    dataframe[self.mol_col].apply(remove_coordinates))
            if self.removeHs:
                dataframe[self.mol_col] = dataframe[self.mol_col].apply(Chem.RemoveHs)
            # render
            dataframe["img"] = dataframe[self.mol_col].apply(self.mol_to_img)
            # cleanup
            if not keep_mols:
                dataframe.drop(columns=self.mol_col, inplace=True)
                self.mol_col = None
        else:
            dataframe["img"] = None
        # generate smiles col if not present or needs to be updated
        if self.mol_col and self.smiles_col not in dataframe.columns:
            dataframe[self.smiles_col] = dataframe[self.mol_col].apply(mol_to_smiles)

    def render(self, template="pages", **kwargs):
        """Returns the HTML document corresponding to the "pages" or "table"
        template. See :meth:`to_pages` and :meth:`to_table` for the full list
        of arguments

        Parameters
        ----------
        template : str
            Kind of grid to draw:

            * table
                A very simple table where all molecules are displayed on the
                document, similarly to RDKit's :func:`~rdkit.Chem.Draw.rdMolDraw2D.MolsToGridImage`.
                This template is mainly used for printing on paper or in a PDF
                file. Most of the interactive actions aren't available.
            * pages
                A more interactive version that layouts the original set of
                molecules on several pages, allows for selecting molecules and
                filtering them using text or substructure queries.
        """
        self.template = template
        return getattr(self, f"to_{self.template}")(**kwargs)

    def to_pages(self, subset=None, tooltip=None, n_cols=5, n_rows=3,
                 border="1px solid #cccccc", gap=0,
                 fontsize="12pt", fontfamily="'DejaVu', sans-serif",
                 textalign="center", tooltip_fmt="<strong>{key}</strong>: {value}",
                 tooltip_trigger="click hover", tooltip_placement="bottom",
                 hover_color="#e7e7e7", style=None, selection=True, transform=None,
                 custom_css=None, custom_header=None, callback=None, sort_by=None,
                 substruct_highlight=True, single_highlight=False):
        """Returns the HTML document for the "pages" template

        Parameters
        ----------
        subset : list or None
            Columns to be displayed in each cell of the grid. Each column's
            value will be displayed from top to bottom in the same order given
            here. Use ``"img"`` for the image of the molecule, and
            ``"mols2grid-id"`` for the molecule's index in your input file.
        tooltip : list or None
            Columns to be displayed as a tooltip when hovering/clicking on the
            image of a cell.
        tooltip_fmt : str
            Format string of each key/value pair in the tooltip
        tooltip_trigger : str
            Sequence of triggers for the tooltip: ``click``, ``hover`` or
            ``focus``
        tooltip_placement : str
            Position of the tooltip: ``auto``, ``top``, ``bottom``, ``left``
            or ``right``
        n_cols : int
            Number of columns per page
        n_rows : int
            Number of rows per page
        border : str
            Styling of the border around each cell (CSS)
        gap : int
            Size of the margin around each cell in px
        fontsize : str
            Font size of the text displayed in each cell (CSS)
        fontfamily : str
            Font used for the text in each cell (CSS)
        textalign : str
            Alignment of the text in each cell (CSS)
        hover_color : str
            Background color when hovering a cell (CSS)
        style : dict or None
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

        selection : bool
            Enables the selection of molecules and displays a checkbox at the
            top of each cell. In the context of a Jupyter notebook, this gives
            you access to your selection (index and SMILES) through :func:`mols2grid.get_selection()`
            or :meth:`MolGrid.get_selection()`. In all cases, you can export your
            selection by clicking on the ☑ icon.
        transform : dict or None
            Functions applied to specific items in all cells. The dict must follow
            a ``key: function`` structure where the key must correspond to one of
            the columns in ``subset`` or ``tooltip``. The function takes the item's
            value as input and transforms it, for example::

                transform={"Solubility": lambda x: f"{x:.2f}",
                           "Melting point": lambda x: f"MP: {5/9*(x-32):.1f}°C"}

            These transformations only affect columns in ``subset`` and
            ``tooltip``, and do not interfere with ``style``.
        custom_css : str or None
            Custom CSS properties applied to the content of the HTML document
        custom_header : str or None
            Custom libraries to be loaded in the header of the document
        callback : str or callable
            Only available for the "pages" template. JavaScript or Python
            callback to be executed when clicking on an image. A dictionnary
            containing the data for the full cell is directly available as
            ``data`` in JS. For Python, the callback function must have
            ``data`` as the first argument to the function. All the values in
            the ``data`` dict are parsed as strings, except "mols2grid-id"
            which is always an integer. Note that fields containing spaces in
            their name will be replaced by hyphens, i.e. "mol weight" becomes
            available as ``data["mol-weight"]``.
        sort_by : str or None
            Sort the grid according to the following field (which must be
            present in ``subset`` or ``tooltip``).
        substruct_highlight : bool
            Highlight substructure when using the SMARTS search. Only available
            when ``prerender=False``
        single_highlight : bool
            Highlight only the first match of the substructure query

        Returns
        -------
        html_document : str

        Notes
        -----
        If ``subset=None, tooltip=None``, the index and image will be directly
        displayed on the grid while the remaining fields will be in the
        tooltip.

        .. versionadded:: 0.1.0
            Added ``sort_by``, ``custom_css``, ``custom_header`` and
            ``callback`` arguments.
            Added the ability to style an entire cell with
            ``style={"__all__": <function>}``.

        .. versionadded:: 0.2.0
            Added ``substruct_highlight`` argument

        .. versionchanged:: 0.2.2
            If both ``subset`` and ``tooltip`` are ``None``, the index and
            image will be directly displayed on the grid while the remaining
            fields will be in the tooltip.
        """
        if substruct_highlight and self.prerender:
            raise ValueError(
                "Cannot highlight substructure search with prerendered images")
        if self.mol_col:
            df = self.dataframe.drop(columns=self.mol_col).copy()
        else:
            df = self.dataframe.copy()
        cell_width = self.img_size[0]
        smiles = self.smiles_col
        content = []
        column_map = {}
        width = n_cols * (cell_width + 2 * (gap + 2))

        if subset is None:
            if tooltip is None:
                subset = ["mols2grid-id", "img"]
                tooltip = [x for x in df.columns.tolist() if x not in subset]
            else:
                subset = df.columns.tolist()
                subset = [subset.pop(subset.index("img"))] + subset

        if "img" not in subset:
            raise KeyError("Please add the 'img' field in the `subset` parameter")

        # define fields that are searchable and sortable
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
        sort_cols = ["mols2grid-id"] + sort_cols
        # get unique list but keep order
        sort_cols = list(dict.fromkeys(sort_cols))
        if style is None:
            style = {}
        if transform is None:
            transform = {}
        value_names = list(set(subset + [smiles] + tooltip))
        value_names = [f"data-{col}" for col in value_names]

        # force id, SMILES, and tooltip values to be present in the data
        final_columns = subset[:]
        final_columns.extend(["mols2grid-id", smiles])
        if tooltip:
            final_columns.extend(tooltip)
        final_columns = list(set(final_columns))

        # make a copy if id shown explicitely
        if "mols2grid-id" in subset:
            id_name = "mols2grid-id-copy"
            df[id_name] = df["mols2grid-id"]
            value_names.append(f"data-{id_name}")
            final_columns.append(id_name)
            subset = [id_name if x == "mols2grid-id" else x for x in subset]
        # organize data
        temp = []
        for col in subset:
            if col == "img" and tooltip:
                s = (f'<a tabindex="0" class="data data-{col} mols2grid-tooltip" '
                     'data-toggle="popover" data-content="."></a>')  
            else:
                if style.get(col):
                    s = (f'<div class="data data-{slugify(col)} style-{slugify(col)}" '
                         'style=""></div>')
                else:
                    s = f'<div class="data data-{slugify(col)}"></div>'
            temp.append(s)
            column_map[col] = f"data-{col}"
        content = temp + content
        # add but hide SMILES div if not present
        if smiles not in (subset + tooltip):
            s = f'<div class="data data-{slugify(smiles)}" style="display: none;"></div>'
            content.append(s)
            column_map[smiles] = f"data-{smiles}"
        # set mapping for list.js
        if "__all__" in style.keys():
            whole_cell_style = True
            x = "[{data: ['mols2grid-id', 'cellstyle']}, "
        else:
            whole_cell_style = False
            x = "[{data: ['mols2grid-id']}, "
        value_names = [slugify(c) for c in value_names]
        value_names = x + str(value_names)[1:]

        # apply CSS styles
        for col, func in style.items():
            if col == "__all__":
                name = "cellstyle"
                df[name] = df.apply(func, axis=1)
            else:
                name = f"style-{slugify(col)}"
                df[name] = df[col].apply(func)
            final_columns.append(name)
            value_names = value_names[:-1] + f", {{ attr: 'style', name: {name!r} }}]"

        if tooltip:
            df["mols2grid-tooltip"] = df.apply(tooltip_formatter, axis=1,
                                               args=(tooltip, tooltip_fmt, style,
                                                     transform))
            final_columns += ["mols2grid-tooltip"]
            value_names = (value_names[:-1] +
                           ", {attr: 'data-content', name: 'mols2grid-tooltip'}]")

        # apply custom user function
        for col, func in transform.items():
            df[col] = df[col].apply(func)

        if selection:
            if self._cached_selection:
                df["cached_checkbox"] = False
                df.loc[df["mols2grid-id"].isin(self._cached_selection.keys()), "cached_checkbox"] = True
                final_columns += ["cached_checkbox"]
                value_names = (value_names[:-1] +
                           ", {attr: 'checked', name: 'cached_checkbox'}]")
            checkbox = ('<input type="checkbox" '
                        'class="position-relative float-left cached_checkbox">')
        else:
            checkbox = ""
        if whole_cell_style:
            item = ('<div class="cell" data-mols2grid-id="0" '
                    'data-cellstyle="0">{checkbox}{content}</div>')
        else:
            item = ('<div class="cell" data-mols2grid-id="0">'
                    '{checkbox}{content}</div>')
        item = item.format(checkbox=checkbox, content="".join(content))

        # callback
        if callable(callback):
            if callback.__name__ == "<lambda>":
                raise TypeError(
                    "Lambda functions are not supported as callbacks. Please "
                    "use a regular function instead.")
            # automatically register callback in Google Colab
            try:
                from google import colab
            except:
                pass
            else:
                colab.output.register_callback(callback.__name__, callback)
            callback_type = "python"
            callback = callback.__name__
        else:
            callback_type = "js"

        if sort_by and sort_by != "mols2grid-id":
            if sort_by in (subset + tooltip):
                sort_by = f"data-{slugify(sort_by)}"
            else:
                raise ValueError(f"{sort_by!r} is not an available field in "
                                 "`subset` or `tooltip`")
        else:
            sort_by = "mols2grid-id"
        
        # slugify remaining vars
        column_map = {k: slugify(v) for k, v in column_map.items()}
        sort_cols = [slugify(c) for c in sort_cols]
        search_cols = [slugify(c) for c in search_cols]
        smiles = slugify(smiles)
        df = df[final_columns].rename(columns=column_map).sort_values(sort_by)

        template = env.get_template('pages.html')
        template_kwargs = dict(
            width = width,
            height = self.img_size[1],
            border = border,
            textalign = textalign,
            cell_width = cell_width,
            fontfamily = fontfamily,
            fontsize = fontsize,
            gap = f"{gap}px",
            hover_color = hover_color,
            item = item,
            item_repr = repr(item),
            value_names = value_names,
            tooltip = tooltip,
            tooltip_trigger = repr(tooltip_trigger),
            tooltip_placement = repr(tooltip_placement),
            n_items_per_page = n_rows * n_cols,
            search_cols = search_cols,
            data = json.dumps(df.to_dict("records")),
            selection = selection,
            cached_selection = ([list(self._cached_selection.keys()),
                                 list(self._cached_selection.values())]
                                if self._cached_selection else False),
            smiles_col = smiles,
            sort_cols = sort_cols,
            grid_id = self._grid_id,
            whole_cell_style = whole_cell_style,
            custom_css = custom_css or "",
            custom_header = custom_header or "",
            callback = callback,
            callback_type = callback_type,
            sort_by = sort_by,
            removeHs = self.removeHs,
            prefer_coordGen = self.prefer_coordGen,
            onthefly = not self.prerender,
            substruct_highlight = substruct_highlight,
            json_draw_opts = getattr(self, "json_draw_opts", ""),
            single_highlight = single_highlight,
        )
        return template.render(**template_kwargs)

    def get_selection(self):
        """Retrieve the dataframe subset corresponding to your selection

        Returns
        -------
        pandas.DataFrame
        """
        sel = list(register.get_selection().keys())
        return (self.dataframe.loc[self.dataframe["mols2grid-id"].isin(sel)]
                              .drop(columns=self._extra_columns))

    def filter(self, mask):
        """Filters the grid using a mask (boolean array)

        Parameters
        ----------
        mask : list, pandas.Series or numpy.ndarray
            Boolean array: `True` when the item should be displayed, `False` if it should
            be filtered out. 
        """
        # convert mask to mols2grid-id
        ids = self.dataframe.loc[mask]["mols2grid-id"]
        return self._filter_by_id(ids)

    def filter_by_index(self, indices):
        """Filters the grid using the dataframe's index"""
        # convert index to mols2grid-id
        ids = self.dataframe.loc[self.dataframe.index.isin(indices)]["mols2grid-id"]
        return self._filter_by_id(ids)

    def _filter_by_id(self, ids):
        """Filters the grid using the values in the ``mols2grid-id`` column"""
        if isinstance(ids, (pd.Series, np.ndarray)):
            ids = ids.to_list()
        code = env.get_template('js/filter.js').render(
            grid_id = self._grid_id,
            ids = ids)
        return Javascript(code)

    def to_table(self, subset=None, tooltip=None, n_cols=5,
                 cell_width=160, border="1px solid #cccccc", gap=0,
                 fontsize="12pt", fontfamily="'DejaVu', sans-serif",
                 textalign="center", tooltip_fmt="<strong>{key}</strong>: {value}",
                 tooltip_trigger="click hover", tooltip_placement="bottom",
                 hover_color="#e7e7e7", style=None, transform=None,
                 custom_css=None, custom_header=None, sort_by=None):
        """Returns the HTML document for the "table" template

        Parameters
        ----------
        subset : list or None
            Columns to be displayed in each cell of the grid. Each column's
            value will be displayed from top to bottom in the same order given
            here. Use ``"img"`` for the image of the molecule, and
            ``"mols2grid-id"`` for the molecule's index in your input file.
        tooltip : list or None
            Columns to be displayed as a tooltip when hovering/clicking on the
            image of a cell.
        tooltip_fmt : str
            Format string of each key/value pair in the tooltip
        tooltip_trigger : str
            Sequence of triggers for the tooltip: ``click``, ``hover`` or
            ``focus``
        tooltip_placement : str
            Position of the tooltip: ``auto``, ``top``, ``bottom``, ``left``
            or ``right``
        n_cols : int
            Number of columns in the table
        border : str
            Styling of the border around each cell (CSS)
        gap : int or str
            Size of the margin around each cell (CSS)
        fontsize : str
            Font size of the text displayed in each cell (CSS)
        fontfamily : str
            Font used for the text in each cell (CSS)
        textalign : str
            Alignment of the text in each cell (CSS)
        hover_color : str
            Background color when hovering a cell (CSS)
        style : dict or None
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

        transform : dict or None
            Functions applied to specific items in all cells. The dict must follow
            a ``key: function`` structure where the key must correspond to one of
            the columns in ``subset`` or ``tooltip``. The function takes the item's
            value as input and transforms it, for example::

                transform={"Solubility": lambda x: f"{x:.2f}",
                           "Melting point": lambda x: f"MP: {5/9*(x-32):.1f}°C"}

            These transformations only affect columns in ``subset`` and
            ``tooltip``, and do not interfere with ``style``.
        custom_css : str or None
            Custom CSS properties applied to the content of the HTML document
        custom_header : str or None
            Custom libraries to be loaded in the header of the document
        sort_by : str or None
            Sort the table according to the following field
        
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
                "Please set `prerender=True` when using the 'table' template")
        tr = []
        data = []
        sort_by = sort_by or "mols2grid-id"
        df = self.dataframe.sort_values(sort_by).reset_index(drop=True)
        cell_width = self.img_size[0]

        if subset is None:
            if tooltip is None:
                subset = ["mols2grid-id", "img"]
                tooltip = [x for x in df.columns.tolist() if x not in subset]
            else:
                subset = df.columns.tolist()
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
            td = [f'<td class="col-{ncol}">']
            if "__all__" in style.keys():
                s = style["__all__"](row)
                div = [f'<div class="cell-{i}" style="{s}">']
            else:
                div = [f'<div class="cell-{i}">']
            for col in subset:
                v = row[col]
                if col == "img" and tooltip:
                    popover = tooltip_formatter(row, tooltip, tooltip_fmt, style,
                                                transform)
                    item = (f'<div class="data data-img mols2grid-tooltip" '
                            f'data-toggle="popover" data-content="{escape(popover)}">'
                            f'{v}</div>')
                else:
                    func = style.get(col)
                    slug_col = slugify(col)
                    if func:
                        item = f'<div class="data data-{slug_col}" style="{func(v)}">'
                    else:
                        item = f'<div class="data data-{slug_col}">'
                    func = transform.get(col)
                    v = func(v) if func else v
                    item += f'{v}</div>'
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

        template = env.get_template('table.html')
        template_kwargs = dict(
            border = border,
            textalign = textalign,
            cell_width = cell_width,
            fontfamily = fontfamily,
            fontsize = fontsize,
            gap = gap,
            hover_color = hover_color,
            tooltip = tooltip,
            tooltip_trigger = repr(tooltip_trigger),
            tooltip_placement = repr(tooltip_placement),
            custom_css = custom_css or "",
            custom_header = custom_header or "",
            data = "\n".join(data),
        )
        return template.render(**template_kwargs)

    @requires("IPython.display")
    def display(self, width="100%", height=None, iframe_allow="clipboard-write",
                **kwargs):
        """Render and display the grid in a Jupyter notebook
        
        Returns
        -------
        view : IPython.core.display.HTML
        """
        doc = self.render(**kwargs)
        iframe = (env.get_template("html/iframe.html")
                     .render(width=width, height=height, padding=18,
                             allow=iframe_allow, doc=escape(doc)))
        return HTML(iframe)

    def save(self, output, **kwargs):
        """Render and save the grid in an HTML document"""
        with open(output, "w", encoding="utf-8") as f:
            f.write(self.render(**kwargs))
