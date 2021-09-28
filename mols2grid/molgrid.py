import warnings
from base64 import b64encode
from html import escape
import json
from pathlib import Path
import pandas as pd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Draw
from jinja2 import Environment, FileSystemLoader
from .utils import (requires,
                    tooltip_formatter,
                    mol_to_record,
                    mol_to_smiles,
                    sdf_to_dataframe,
                    remove_coordinates)
from .select import register
try:
    from IPython.display import HTML, Javascript
except ModuleNotFoundError:
    pass
else:
    warnings.filterwarnings("ignore",
                            "Consider using IPython.display.IFrame instead")


env = Environment(loader=FileSystemLoader(Path(__file__).parent / 'templates'),
                  autoescape=False)


class MolGrid:
    """Class that handles drawing molecules, rendering the HTML document and
    saving or displaying it in a notebook
    """

    def __init__(self, df, smiles_col="SMILES", mol_col=None, removeHs=False,
        use_coords=True, coordGen=True, useSVG=True, size=(160, 120),
        MolDrawOptions=None, rename=None, name="default", **kwargs):
        """
        Parameters
        ----------
        df : pandas.DataFrame or dict or list
            Dataframe containing a SMILES or mol column, or dictionary
            containing a list of SMILES, or list of dictionnaries containing a
            SMILES field
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
        MolDrawOptions : rdkit.Chem.Draw.MolDrawOptions or None
            Drawing options. Useful for making highly customized drawings
        rename : dict or None
            Rename the properties/fields stored in the molecule
        name : str
            Name of the grid. Used when retrieving selections from multiple
            grids at the same time
        kwargs : object
            MolDrawOptions attributes

        Notes
        -----
        The list of supported MolDrawOptions attributes are available in
        https://www.rdkit.org/docs/source/rdkit.Chem.Draw.rdMolDraw2D.html#rdkit.Chem.Draw.rdMolDraw2D.MolDrawOptions

        ..versionchanged: 0.0.7
            Added `rename` argument to replace `mapping`
        """
        if not (smiles_col or mol_col):
            raise ValueError("One of `smiles_col` or `mol_col` must be set")
        if not isinstance(name, str):
            raise TypeError(
                f"`name` must be a string. Currently of type {type(name).__name__}")
        Draw.rdDepictor.SetPreferCoordGen(coordGen)
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
        # generate temporary RDKit molecules
        if smiles_col and not mol_col:
            mol_col = "mol"
            keep_mols = False
            self._extra_columns.append(mol_col)
            dataframe[mol_col] = dataframe[smiles_col].apply(Chem.MolFromSmiles)
        else:
            keep_mols = True
        # remove hydrogens
        if removeHs:
            dataframe[mol_col] = dataframe[mol_col].apply(Chem.RemoveHs)
        if not use_coords:
            dataframe[mol_col] = dataframe[mol_col].apply(remove_coordinates)
        # generate smiles col
        if mol_col and (smiles_col not in dataframe.columns):
            dataframe[smiles_col] = dataframe[mol_col].apply(mol_to_smiles)
        # add index
        dataframe["mols2grid-id"] = list(range(len(dataframe)))
        # drop None
        dataframe.dropna(axis=0, subset=[mol_col], inplace=True)
        # generate drawings
        self.useSVG = useSVG
        opts = MolDrawOptions or Draw.MolDrawOptions()
        for key, value in kwargs.items():
            setattr(opts, key, value)
        self.MolDrawOptions = opts
        self._MolDraw2D = Draw.MolDraw2DSVG if useSVG else Draw.MolDraw2DCairo
        self.img_size = size
        dataframe["img"] = dataframe[mol_col].apply(self.mol_to_img)
        if keep_mols:
            self.dataframe = dataframe
        else:
            self.dataframe = dataframe.drop(columns=mol_col)
            mol_col = None
        self.smiles_col = smiles_col
        self.mol_col = mol_col
        # register instance
        self._grid_id = name
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
            Path to the SDF file
        kwargs : object
            Other arguments passed on initialization
        """
        mol_col = kwargs.pop("mol_col", "mol")
        df = sdf_to_dataframe(sdf_file, mol_col=mol_col)
        return cls(df, mol_col=mol_col, **kwargs)

    @property
    def template(self):
        """Kind of grid displayed, one of:
            - pages
            - table
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

    def render(self, template="pages", **kwargs):
        """Returns the HTML document corresponding to the "pages" or "table"
        template. See `to_pages` and `to_table` for the list of arguments
        
        Parameters
        ----------
        template : str
            Kind of grid to draw:
                * "table" is a very simple table where all molecules are
                  displayed on the document, the main usecase is printing to
                  PDF or on paper.
                * "pages" is a more interactive version that splits the
                  original data into several pages.
        """
        self.template = template
        return getattr(self, f"to_{self.template}")(**kwargs)
    
    def to_pages(self, subset=None, tooltip=None,
                 cell_width=160, n_cols=5, n_rows=3,
                 border="1px solid #cccccc", gap=0,
                 fontsize="12pt", fontfamily="'DejaVu', sans-serif",
                 textalign="center", tooltip_fmt="<strong>{key}</strong>: {value}",
                 tooltip_trigger="click hover", tooltip_placement="bottom",
                 hover_color="#e7e7e7", style=None, selection=True, transform=None,
                 custom_css=None, callback=None):
        """Returns the HTML document for the "pages" template
        
        Parameters
        ----------
        subset : list or None
            Columns to be displayed in each cell of the grid. Each
            column's value will be displayed from top to bottom in the same
            order given here. Use `"img"` for the image of the molecule.
            Default: all columns (with "img" in first position)
        tooltip : list or None
            Columns to be displayed as a tooltip when hovering/clicking on the
            image of a cell. Use `None` for no tooltip.
        tooltip_fmt : str
            Format string of each key/value pair in the tooltip
        tooltip_trigger : str
            Sequence of triggers for the tooltip: (click, hover, focus)
        tooltip_placement : str
            Position of the tooltip: auto, top, bottom, left, right
        n_cols : int
            Number of columns per page
        n_rows : int
            Number of rows per page
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
            CSS styling applied to specific items in all cells. The dict must follow a
            `key: function` structure where the key must correspond to one of
            the columns in `subset` or `tooltip`. The function takes the item's value as
            input, and outputs a valid CSS styling, for example
            `style={"Solubility": lambda x: "color: red" if x < -5 else "color: black"}`
            if you want to color the text corresponding to the "Solubility"
            column in your dataframe. You can also style a whole cell using the `__all__`
            key, the corresponding function then has access to all values for each cell:
            `style={"__all__": lambda x: "color: red" if x["Solubility"] < -5 else 
            "color: black"}`
        selection : bool
            Enables the selection of molecules and displays a checkbox at the top of each
            cell. This is only usefull in the context of a Jupyter notebook, which gives
            you access to your selection (index and SMILES) through `mols2grid.selection`
        transform : dict or None
            Functions applied to specific items in all cells. The dict must follow a
            `key: function` structure where the key must correspond to one of the columns
            in `subset` or `tooltip`. The function takes the item's value as input and 
            transforms it, for example:
            `transform={"Solubility": lambda x: f"{x:.2f}",
                        "Melting point": lambda x: f"MP: {5/9*(x-32):.1f}°C"}`
            will round the solubility to 2 decimals, and display the melting point in
            Celsius instead of Fahrenheit with a single digit precision and some text
            before (MP) and after (°C) the value. These transformations only affect
            columns in `subset` and `tooltip`, and do not interfere with `style`.
        custom_css : str
            Custom CSS properties applied to the content of the HTML document.
        callback : str
            JavaScript callback to be executed when clicking on an image. A dictionnary
            containing the data for the full cell is available as `data`. All the values
            are parsed as strings, except "mols2grid-id" which is always an integer.
        """
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
            subset = df.columns.tolist()
            subset = [subset.pop(subset.index("img"))] + subset
        # define fields that are searchable and sortable
        search_cols = [f"data-{col}" for col in subset if col != "img"]
        if tooltip:
            search_cols.append("mols2grid-tooltip")
            sort_cols = search_cols[:-1]
            sort_cols.extend([f"data-{col}" for col in tooltip])
            for col in tooltip:
                if col not in subset:
                    s = f'<div class="data data-{col}" style="display: none;"></div>'
                    content.append(s)
                    column_map[col] = f"data-{col}"
        else:
            sort_cols = search_cols[:]
        sort_cols = ["mols2grid-id"] + sort_cols
        # get unique list but keep order
        sort_cols = list(dict.fromkeys(sort_cols))
        if style is None:
            style = {}
        if transform is None:
            transform = {}
        if tooltip:
            value_names = list(set(subset + [smiles] + tooltip))
        else:
            value_names = list(set(subset + [smiles]))
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
        for col in subset:
            if col == "img" and tooltip:
                s = (f'<a tabindex="0" class="data data-{col} mols2grid-tooltip" '
                     'data-toggle="popover" data-content="foo"></a>')  
            else:
                if style.get(col):
                    s = f'<div class="data data-{col} style-{col}" style=""></div>'
                else:
                    s = f'<div class="data data-{col}"></div>'
            content.append(s)
            column_map[col] = f"data-{col}"
        # add but hide SMILES div if not present
        if smiles not in subset:
            s = f'<div class="data data-{smiles}" style="display: none;"></div>'
            content.append(s)
            column_map[smiles] = f"data-{smiles}"
        # set mapping for list.js
        if "__all__" in style.keys():
            whole_cell_style = True
            x = "[{data: ['mols2grid-id', 'cellstyle']}, "
        else:
            whole_cell_style = False
            x = "[{data: ['mols2grid-id']}, "
        value_names = x + str(value_names)[1:]

        # apply CSS styles
        for col, func in style.items():
            if col == "__all__":
                name = "cellstyle"
                df[name] = df.apply(func, axis=1)
            else:
                name = f"style-{col}"
                df[name] = df[col].apply(func)
            final_columns.append(name)
            value_names = value_names[:-1] + f", {{ attr: 'style', name: {name!r} }}]"

        if tooltip:
            df["mols2grid-tooltip"] = df.apply(tooltip_formatter, axis=1,
                                               args=(tooltip, tooltip_fmt, style,
                                                     transform))
            final_columns = final_columns + ["mols2grid-tooltip"]
            value_names = (value_names[:-1] +
                           ", {attr: 'data-content', name: 'mols2grid-tooltip'}]")

        # apply custom user function
        for col, func in transform.items():
            df[col] = df[col].apply(func)

        if selection:
            checkbox = '<input type="checkbox" class="position-relative float-left">'
        else:
            checkbox = ""
        if whole_cell_style:
            item = ('<div class="cell" data-mols2grid-id="0" '
                    'data-cellstyle="0">{checkbox}{content}</div>')
        else:
            item = ('<div class="cell" data-mols2grid-id="0">'
                    '{checkbox}{content}</div>')
        item = item.format(checkbox=checkbox, content="".join(content))

        df = df[final_columns].rename(columns=column_map)

        template = env.get_template('pages.html')
        template_kwargs = dict(
            width = width,
            border = border,
            textalign = textalign,
            cell_width = cell_width,
            fontfamily = fontfamily,
            fontsize = fontsize,
            gap = gap,
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
            smiles_col = smiles,
            sort_cols = sort_cols,
            grid_id = self._grid_id,
            whole_cell_style = whole_cell_style,
            custom_css = custom_css if custom_css else "",
            callback = callback,
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
        mask : list, pd.Series, np.ndarray
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
        """Filters the grid using the values in the `mols2grid-id` column"""
        if isinstance(ids, (pd.Series, np.ndarray)):
            ids = ids.to_list()
        code = env.get_template('js/filter.js').render(
            grid_id = self._grid_id,
            ids = ids)
        return Javascript(code)
    
    def to_table(self, subset=None, tooltip=None, n_cols=6,
                 cell_width=160, border="1px solid #cccccc", gap=0,
                 fontsize="12pt", fontfamily="'DejaVu', sans-serif",
                 textalign="center", tooltip_fmt="<strong>{key}</strong>: {value}",
                 tooltip_trigger="click hover", tooltip_placement="bottom",
                 hover_color="#e7e7e7", style=None, transform=None):
        """Returns the HTML document for the "table" template
        
        Parameters
        ----------
        subset : list or None
            Columns to be displayed in each cell of the grid. Each
            column's value will be displayed from top to bottom in the same
            order given here. Use `"img"` for the image of the molecule.
            Default: all columns (with "img" in first position)
        tooltip : list or None
            Columns to be displayed as a tooltip when hovering/clicking on the
            image of a cell. Use `None` for no tooltip.
        tooltip_fmt : str
            Format string of each key/value pair in the tooltip
        tooltip_trigger : str
            Sequence of triggers for the tooltip: (click, hover, focus)
        tooltip_placement : str
            Position of the tooltip: auto, top, bottom, left, right
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
            CSS styling applied to specific items in all cells. The dict must follow a
            `key: function` structure where the key must correspond to one of
            the columns in `subset` or `tooltip`. The function takes the item's value as
            input, and outputs a valid CSS styling, for example
            `style={"Solubility": lambda x: "color: red" if x < -5 else "color: black"}`
            if you want to color the text corresponding to the "Solubility"
            column in your dataframe
        transform : dict or None
            Functions applied to specific items in all cells. The dict must follow a
            `key: function` structure where the key must correspond to one of the columns
            in `subset`. The function takes the item's value as input and transforms it,
            for example:
            `transform={"Solubility": lambda x: f"{x:.2f}",
                    "Melting point": lambda x: f"MP: {5/9*(x-32):.1f}°C"}`
            will round the solubility to 2 decimals, and display the melting point in
            Celsius instead of Fahrenheit with a single digit precision and some text
            before (MP) and after (°C) the value. These transformations only affect
            columns in `subset` and `tooltip`, and are applied independantly from `style`
        """
        tr = []
        data = []
        df = self.dataframe
        cell_width = self.img_size[0]

        if subset is None:
            subset = df.columns.tolist()
            subset = [subset.pop(subset.index("img"))] + subset
        if style is None:
            style = {}
        if transform is None:
            transform = {}

        for i, row in df.iterrows():
            ncell = i + 1
            nrow, ncol = divmod(i, n_cols)
            td = [f'<td class="col-{ncol}>"']
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
                    func = transform.get(col)
                    v = func(v) if func else v
                    item = (f'<div class="data data-{col} mols2grid-tooltip" data-toggle="popover" '
                            f'data-content="{escape(popover)}">{v}</div>')
                else:
                    func = style.get(col)
                    if func:
                        item = f'<div class="data data-{col}" style="{func(v)}">'
                    else:
                        item = f'<div class="data data-{col}">'
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
            data = "\n".join(data),
        )
        return template.render(**template_kwargs)

    @requires("IPython.display")
    def display(self, width="100%", height=None, **kwargs):
        """Render and display the grid in a Jupyter notebook"""
        code = self.render(**kwargs)
        if height:
            iframe = (
                '<iframe class="mols2grid-iframe" frameborder="0" '
                'width={width} height="{height}" srcdoc="{code}"></iframe>')
        else:
            # automatically resize iframe
            iframe = (
                '<iframe class="mols2grid-iframe" frameborder="0" '
                'onload="javascript:(function(o){{o.style.height=o.contentWindow.document.body.scrollHeight+18+\'px\';}}(this));" '
                'width={width} height="10px" srcdoc="{code}"></iframe>'
            )
        return HTML(iframe.format(width=width, height=height,
                                  code=escape(code)))

    def save(self, output, **kwargs):
        """Render and save the grid in an HTML document"""
        with open(output, "w") as f:
            f.write(self.render(**kwargs))