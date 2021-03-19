import warnings
from base64 import b64encode
from html import escape
from pathlib import Path
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
from jinja2 import Environment, FileSystemLoader
from .utils import requires, tooltip_formatter
try:
    from IPython.display import HTML
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
    def __init__(self, df, smiles_col="SMILES", coordGen=True, useSVG=True,
        **kwargs):
        """
        Parameters
        ----------
        df : pandas.DataFrame
            Dataframe containing a SMILES column and some other information 
            about each molecule
        smiles_col : str
            Name of the SMILES column in the dataframe
        coordGen : bool
            Sets whether or not the CoordGen library should be preferred to the
            RDKit depiction library
        useSVG : bool
            Use SVG instead of PNG
        kwargs : object
            Arguments passed to the drawing function
        """
        Draw.rdDepictor.SetPreferCoordGen(coordGen)
        self.useSVG = useSVG
        dataframe = df.copy()
        dataframe["img"] = dataframe[smiles_col].apply(self.smi_to_img,
                                                       **kwargs)
        self.dataframe = dataframe

    @classmethod
    def from_mols(cls, mols, mapping=None, **kwargs):
        """Set up the dataframe used by mols2grid directly from a list of RDKit
        molecules.
        
        Parameters
        ----------
        mols : list
            List of RDKit molecules
        mapping : dict or None
            Rename the properties stored in the molecule
        kwargs : object
            Other arguments passed on initialization
        """
        df = pd.DataFrame([{"SMILES": Chem.MolToSmiles(mol),
                            **mol.GetPropsAsDict()}
                           for mol in mols if mol])
        if mapping:
            df.rename(columns=mapping, inplace=True)
        return cls(df, **kwargs)

    @classmethod
    def from_sdf(cls, sdf_file, mapping=None, **kwargs):
        """Set up the dataframe used by mols2grid directly from an SDFile
        
        Parameters
        ----------
        sdf_file : str
            Path to the SDF file
        mapping : dict or None
            Rename fields of the SDFile
        kwargs : object
            Other arguments passed on initialization
        """
        df = pd.DataFrame([{"SMILES": Chem.MolToSmiles(mol),
                            **mol.GetPropsAsDict()}
                           for mol in Chem.SDMolSupplier(sdf_file) if mol])
        if mapping:
            df.rename(columns=mapping, inplace=True)
        return cls(df, **kwargs)

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

    def draw_mol(self, mol, size=(160, 120), **kwargs):
        """Draw a molecule

        Parameters
        ----------
        mol : rdkit.Chem.rdchem.Mol
            The molecule to draw
        size : tuple
            The size of the drawing canvas
        **kwargs : object
            Attributes of the rdkit.Chem.Draw.rdMolDraw2D.MolDrawOptions class
            like `fixedBondLength=35, bondLineWidth=2`
        
        Notes
        -----
        The list of supported MolDrawOptions attributes are available in
        https://www.rdkit.org/docs/source/rdkit.Chem.Draw.rdMolDraw2D.html#rdkit.Chem.Draw.rdMolDraw2D.MolDrawOptions
        """
        if self.useSVG:
            d2d = Draw.MolDraw2DSVG(*size)
        else:
            d2d = Draw.MolDraw2DCairo(*size)
        opts = Draw.MolDrawOptions()
        for key, value in kwargs.items():
            setattr(opts, key, value)
        d2d.SetDrawOptions(opts)
        d2d.DrawMolecule(mol)
        d2d.FinishDrawing()
        return d2d.GetDrawingText()
    
    def mol_to_img(self, mol, **kwargs):
        """Convert an RDKit mol to an HTML img tag containing a drawing of the
        molecule"""
        img = self.draw_mol(mol, **kwargs)
        if self.useSVG:
            return img
        data = b64encode(img).decode()
        return f'<img src="data:image/png;base64,{data}" alt="molecule">'
    
    def smi_to_img(self, smi, **kwargs):
        """Convert a SMILES string to an HTML img tag containing a drawing of
        the molecule"""
        mol = Chem.MolFromSmiles(smi)
        return self.mol_to_img(mol, **kwargs)
    
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
                 hover_color="#e7e7e7", style=None):
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
        cell_width : int
            Max width of each cell, in pixels
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
            CSS styling applied to each item in a cell. The dict must follow a
            `key: function` structure where the key must correspond to one of
            the columns in `subset` or `tooltip`. The function takes the item's value as
            input, and outputs a valid CSS styling, for example
            `style={"Solubility": lambda x: "color: red" if x < -5 else "color: black"}`
            if you want to color the text corresponding to the "Solubility"
            column in your dataframe.
        """
        df = self.dataframe.copy()
        if subset is None:
            subset = df.columns.tolist()
            subset = [subset.pop(subset.index("img"))] + subset
        if style is None:
            style = {}
        columns = [f"data-{col}" for col in subset]
        width = n_cols * (cell_width + 2 * (gap + 2))
        content = []
        final_columns = ["mols2grid-id"] + subset
        column_map = {}
        value_names = "[{data: ['mols2grid-id']}, " + str(columns)[1:]
        for col in subset:
            if col == "img" and tooltip:
                s = (f'<div class="data data-{col} mols2grid-tooltip" '
                     'data-toggle="popover" data-content="foo"></div>')  
            else:
                if style.get(col):
                    s = f'<div class="data data-{col} style-{col}" style=""></div>'
                else:
                    s = f'<div class="data data-{col}"></div>'
            content.append(s)
            column_map[col] = f"data-{col}"
            
        if tooltip:
            df["mols2grid-tooltip"] = df.apply(tooltip_formatter, axis=1,
                                             args=(tooltip, tooltip_fmt, style))
            final_columns = final_columns + ["mols2grid-tooltip"]
            value_names = (value_names[:-1] +
                           ", {attr: 'data-content', name: 'mols2grid-tooltip'}]")
            
        for col, func in style.items():
            name = f"style-{col}"
            df[name] = df[col].apply(func)
            final_columns.append(name)
            value_names = value_names[:-1] + f", {{ attr: 'style', name: {name!r} }}]"
        
        item = f'<div class="cell" data-mols2grid-id="0">{"".join(content)}</div>'
        df["mols2grid-id"] = [str(i) for i in range(len(df))]
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
            data = df.to_dict("records"),
        )
        return template.render(**template_kwargs)
    
    def to_table(self, subset=None, tooltip=None, n_cols=6,
                 cell_width=160, border="1px solid #cccccc", gap=0,
                 fontsize="12pt", fontfamily="'DejaVu', sans-serif",
                 textalign="center", tooltip_fmt="<strong>{key}</strong>: {value}",
                 tooltip_trigger="click hover", tooltip_placement="bottom",
                 hover_color="#e7e7e7", style=None):
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
        cell_width : int
            Max width of each cell, in pixels
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
            `key: function` structure where the key must correspond to one of
            the columns in `subset` or `tooltip`. The function takes the item's value as
            input, and outputs a valid CSS styling, for example
            `style={"Solubility": lambda x: "color: red" if x < -5 else "color: black"}`
            if you want to color the text corresponding to the "Solubility"
            column in your dataframe.
        """
        tr = []
        data = []
        df = self.dataframe

        if subset is None:
            subset = df.columns.tolist()
            subset = [subset.pop(subset.index("img"))] + subset
        if style is None:
            style = {}

        for i, row in df.iterrows():
            ncell = i + 1
            nrow, ncol = divmod(i, n_cols)
            td = [f'<td class="col-{ncol}>"']
            div = [f'<div class="cell-{i}">']
            for col in subset:
                v = row[col]
                if col == "img" and tooltip:
                    popover = tooltip_formatter(row, tooltip, tooltip_fmt, style)
                    item = (f'<div class="data data-{col} mols2grid-tooltip" data-toggle="popover" '
                            f'data-content="{escape(popover)}">{v}</div>')
                else:
                    func = style.get(col)
                    if func:
                        item = f'<div class="data data-{col}" style="{func(v)}">{v}</div>'
                    else:
                        item = f'<div class="data data-{col}">{v}</div>'
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
    def display(self, width="100%", height=600, **kwargs):
        """Render and display the grid in a Jupyter notebook"""
        code = self.render(**kwargs)
        iframe = ('<iframe class="mols2grid-iframe" width={width} '
                  'height="{height}" frameborder="0" srcdoc="{code}"></iframe>')
        return HTML(iframe.format(width=width, height=height,
                                  code=escape(code)))

    def save(self, output, **kwargs):
        """Render and save the grid in an HTML document"""
        with open(output, "w") as f:
            f.write(self.render(**kwargs))