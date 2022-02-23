# ![mols2grid logo](https://user-images.githubusercontent.com/27850535/154588465-43dc5d82-ee2d-4178-a2f3-e06000bc87c9.png) mols2grid

[![Pypi version](https://img.shields.io/pypi/v/mols2grid.svg)](https://pypi.python.org/pypi/mols2grid)
![Conda version](https://img.shields.io/conda/vn/conda-forge/mols2grid)

[![Tests status](https://github.com/cbouy/mols2grid/workflows/CI/badge.svg)](https://github.com/cbouy/mols2grid/actions/workflows/ci.yml)
[![Code coverage](https://codecov.io/gh/cbouy/mols2grid/branch/master/graph/badge.svg?token=QDI1XQSDUC)](https://codecov.io/gh/cbouy/mols2grid)
[![Build status](https://github.com/cbouy/mols2grid/workflows/build/badge.svg)](https://github.com/cbouy/mols2grid/actions/workflows/build.yml)

[![Powered by RDKit](https://img.shields.io/static/v1?label=Powered%20by&message=RDKit&color=3838ff&style=flat&logo=data:image/x-icon;base64,AAABAAEAEBAQAAAAAABoAwAAFgAAACgAAAAQAAAAIAAAAAEAGAAAAAAAAAMAABILAAASCwAAAAAAAAAAAADc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nz/FBT/FBT/FBT/FBT/FBT/FBTc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nz/FBT/PBT/PBT/PBT/PBT/PBT/PBT/FBTc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nz/FBT/PBT/ZGT/ZGT/ZGT/ZGT/ZGT/ZGT/PBT/FBTc3Nzc3Nzc3Nzc3Nzc3Nz/FBT/PBT/ZGT/ZGT/ZGT/ZGT/ZGT/ZGT/ZGT/ZGT/PBT/FBTc3Nzc3Nzc3Nz/FBT/PBT/ZGT/ZGT/ZGT/jIz/jIz/jIz/jIz/ZGT/ZGT/ZGT/PBT/FBTc3Nzc3Nz/FBT/PBT/ZGT/ZGT/jIz/jIz/jIz/jIz/jIz/jIz/ZGT/ZGT/PBT/FBTc3Nzc3Nz/FBT/PBT/ZGT/ZGT/jIz/jIz/tLT/tLT/jIz/jIz/ZGT/ZGT/PBT/FBTc3Nzc3Nz/FBT/PBT/ZGT/ZGT/jIz/jIz/tLT/tLT/jIz/jIz/ZGT/ZGT/PBT/FBTc3Nzc3Nz/FBT/PBT/ZGT/ZGT/jIz/jIz/jIz/jIz/jIz/jIz/ZGT/ZGT/PBT/FBTc3Nzc3Nz/FBT/PBT/ZGT/ZGT/ZGT/jIz/jIz/jIz/jIz/ZGT/ZGT/ZGT/PBT/FBTc3Nzc3Nzc3Nz/FBT/PBT/ZGT/ZGT/ZGT/ZGT/ZGT/ZGT/ZGT/ZGT/PBT/FBTc3Nzc3Nzc3Nzc3Nzc3Nz/FBT/PBT/ZGT/ZGT/ZGT/ZGT/ZGT/ZGT/PBT/FBTc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nz/FBT/PBT/PBT/PBT/PBT/PBT/PBT/FBTc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nz/FBT/FBT/FBT/FBT/FBT/FBTc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nzc3Nz/////+B////AP///gB///wAP//4AB//+AAf//gAH//4AB//+AAf//gAH//8AD///gB///8A////gf////////)](https://www.rdkit.org/)
[![Knime Hub](https://img.shields.io/static/v1?label=Available%20on&message=KNIME&color=ffd500&style=flat&logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABAAAAAQCAMAAAAoLQ9TAAAABGdBTUEAALGPC/xhBQAAAAFzUkdCAK7OHOkAAABdUExURUxpcf/VAP/VAP/VAP/VAP/VAP/VAP/VAP/VAP/VAP/VAP/VAP/VAP/VAP/VAP/VAP/VAP/VAP/VAP/VAP/VAP/VAP/VAP/VAP/VAP/VAP/VAP/VAP/VAP/VAP/VAPfHMOMAAAAfdFJOUwCyq6CEYtAEApEspncGDpjxVlAYgDSdiEBHbMrCHtmwXwB/AAAAT3pUWHRSYXcgcHJvZmlsZSB0eXBlIGlwdGMAAHic48osKEnmUgADIwsuYwsTIxNLkxQDEyBEgDTDZAMjs1Qgy9jUyMTMxBzEB8uASKBKLgAolQ7jKiJtHwAAAIxJREFUGNNdjFkSgyAQBYdtYADZVNxz/2NGjSlj+q9fvWqAD1rDk1Ke3nJqH4NnpH7d4iCFvV1XVJ3r7u6URPZiHb8eeFJ25sjDNahlKRDUkq7u5njd32ZC3A433g2+h3bKCuUx9FHOecyV/CzXfi/KSJG9EjJB0lEAS9UxxriINMiOLJim0SfNiYF/3szTBp6mEP9HAAAAAElFTkSuQmCC)](https://hub.knime.com/cbouy/spaces/Private/latest/Interactive%20Grid%20of%20Molecules~OZIyk4YLBNvXq-xU)
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/cbouy/mols2grid/blob/master/demo.ipynb)

mols2grid is an interactive chemical viewer for 2D structures of small molecules, based on RDKit.

![Demo showing mols2grid's integration in a Jupyter notebook](https://raw.githubusercontent.com/cbouy/mols2grid/master/demo.png)

➡️ Try the [demo notebook](https://colab.research.google.com/github/cbouy/mols2grid/blob/master/demo.ipynb) on Google Colab, or the more [advanced notebook](https://colab.research.google.com/github/rdkit/UGM_2021/blob/main/Notebooks/Bouysset_mols2grid.ipynb) from the RDKit UGM 2021.

## 🐍 Installation
---

mols2grid was developped for Python 3.7+ and requires rdkit (>=2020.03.1), pandas and jinja2 as dependencies.  
The easiest way to install it is from conda:
```shell
conda install -c conda-forge mols2grid
```

Alternatively, you can also use pip:
```shell
pip install mols2grid
```

It is fully compatible with Jupyter Notebook and Google Colab and can run on Streamlit.  
It also works with Visual Studio notebooks and Jupyterlab, except for accessing your selection with `mols2grid.get_selection()` and copying it to the clipboard, but you can still export it as a CSV or SMILES file.

<img alt="knime logo" align="left" style="padding:6px" src="https://www.knime.com/sites/default/files/favicons/favicon-32x32.png"/>
<p>You can also use mols2grid directly in <a href="https://www.knime.com/">KNIME</a>, by looking for the `Interactive Grid of Molecules` component on the Knime HUB.<br/>
Make sure you have setup <a href="https://docs.knime.com/latest/python_installation_guide">Knime's Python integration</a> for the node to work.</p>

## 📜 Usage
---

```python
import mols2grid

mols2grid.display("path/to/molecules.sdf",
                  # RDKit's MolDrawOptions parameters
                  fixedBondLength=25,
                  # rename fields for the output document
                  rename={"SOL": "Solubility",
                          "SOL_classification": "Class",
                          "NAME": "Name"},
                  # set what's displayed on the grid
                  subset=["ID", "img", "Solubility"],
                  # set what's displayed on the tooltips
                  tooltip=["Name", "smiles", "Class", "Solubility"],
                  # style for the grid labels and tooltips
                  style={"Solubility": lambda x: "color: red" if x < -3 else "color: black"},
                  # change the precision and format (or other transformations)
                  transform={"Solubility": lambda x: f"{x:+.2f}"})
```

#### Input parameters

You can setup the grid from various inputs:
* a pandas **DataFrame** (with a column of SMILES or RDKit molecules, controlled by the `smiles_col` and `mol_col` parameters),
* a list of **RDKit molecules** (with properties accessible through the `mol.GetPropsAsDict()` method),
* or an **SDF file** (`.sdf` or `.sdf.gz`)

You can also rename each field of your input with the `rename` parameter. Please note that 3 fields are automatically added regardless of your input: `mols2grid-id`, `SMILES` and `img`. If a "SMILES" field already exists, it will not be overwritten.

#### Parameters for the drawing of each molecule

* `useSVG=True`: use SVG images or PNG
* `coordGen=True`: use the coordGen library instead of the RDKit one to depict the molecules in 2D
* `removeHs=False`: remove explicit hydrogen atoms from the drawings
* `size=(160, 120)`: size of each image
* `use_coords=False`: use the coordinates of the input molecules if available
* `MolDrawOptions=None`: RDKit's [MolDrawOptions](https://www.rdkit.org/docs/source/rdkit.Chem.Draw.rdMolDraw2D.html#rdkit.Chem.Draw.rdMolDraw2D.MolDrawOptions) class. Useful for making highly customized drawings. You can also leave this to `None`, and directly use the attributes of this class as parameters like `addStereoAnnotation=True`. Additionally, `atomColourPalette` is made available to modify the color palette if you are not prerendering images.
* `prerender=False`: Prerender images for the entire dataset, or generate them on-the-fly when needed
* `cache_selection=False` : Restores the selection from a previous grid with the same name.

#### Parameters for the grid
  
You can control the general look of the document through the `template` argument:
* `template="pages"` (default) which is displayed above. It integrates nicely with Jupyter notebooks and has a search bar
* `template="table"`, which displays the full list of molecules (no pages). Useful if you ever need to print the full list of molecules on paper (or print to PDF)

Both templates can be configured with the same parameters (a lot of which are [CSS](https://www.w3schools.com/cssref/) declarations).
For the `pages` template, the following parameters are available:

* `subset=None`: list or None  
    Columns to be displayed in each cell of the grid. Each column's value will be displayed from top to bottom in the same order given here. Use `"img"` for the image of the molecule. Default: all columns (with "img" in first position)
* `tooltip=None`: list or None  
    Columns to be displayed as a tooltip when hovering/clicking on the image of a cell. Use `None` for no tooltip.
* `tooltip_fmt="<strong>{key}</strong>: {value}"`: str  
    Format string of each key/value pair in the tooltip
* `tooltip_trigger="click hover"`: str  
    Sequence of triggers for the tooltip: (click, hover, focus)
* `tooltip_placement="bottom"`: str  
    Position of the tooltip: auto, top, bottom, left, right
* `n_cols=5`: int  
    Number of columns per page
* `n_rows=3` : int  
    Number of rows per page
* `border="1px solid #cccccc"`: str  
    Styling of the border around each cell (CSS)
* `gap=0`: int
    Size of the margin around each cell in px
* `fontsize="12pt"`: str  
    Font size of the text displayed in each cell (CSS)
* `fontfamily"'DejaVu', sans-serif"`: str  
    Font used for the text in each cell (CSS)
* `textalign="center"`: str  
    Alignment of the text in each cell (CSS)
* `hover_color="#e7e7e7"`: str  
    Background color when hovering a cell (CSS)
* `style=None`: dict or None  
    CSS styling applied to each item in a cell. The dict must follow a `key: function` structure where the key must correspond to one of the columns in `subset` or `tooltip`. The function takes the item's value as input, and outputs a valid CSS styling. For example, if you want to color the text corresponding to the "Solubility"
    column in your dataframe:
    ```python
    style={"Solubility": lambda x: "color: red" if x < -3 else ""}
    ```
    You can also style a whole cell using `__all__` as a key, the corresponding function then has access to all values for each cell:
    ```python
    style={"__all__": lambda x: "background-color: yellow" if x["Solubility"] < -5 else ""}
    ```
* `transform=None`: dict or None  
    Functions applied to specific items in all cells. The dict must follow a `key: function` structure where the key must correspond to one of the columns in `subset` or `tooltip`. The function takes the item's value as input and transforms it. For example, to round the "Solubility" to 2 decimals, and display the "Melting point" in Celsius instead of Fahrenheit with a single digit precision and some text before ("MP") and after ("°C") the value:
    ```python
    transform={"Solubility": lambda x: f"{x:.2f}",
               "Melting point": lambda x: f"MP: {5/9*(x-32):.1f}°C"}
    ```
    These transformations only affect columns in `subset` and `tooltip` and do not interfere with `style`.
* `selection=True` : bool  
    Enables the selection of molecules using a checkbox. Only usefull in the context of a Jupyter notebook. You can retrieve your selection of molecules (index and SMILES) through `mols2grid.get_selection()`
* `cache_selection=False` : bool
    Restores the selection from a previous grid with the same name if `True`
* `custom_css=None` : str or None  
    Custom CSS properties applied to the content of the HTML document
* `custom_header=None` : str or None  
    Custom libraries (CSS or JS) to be loaded in the header of the document
* `callback=None` : str or callable  
    JavaScript or Python callback to be executed when clicking on an image. A dictionnary containing the data for the full cell is directly available as `data` in JS. For Python, the callback function must have `data` as the first argument to the function. All the values in the `data` dict are parsed as strings, except "mols2grid-id" which is always an integer.
* `sort_by=None` : str or None  
    Sort the grid according to the following field (which must be present in `subset` or `tooltip`).
* `substruct_highlight=True` : bool
    Highlight the query when using a SMARTS substructure search. Only available when `prerender=False`.

Less options are available for the `table` template, you can check the complete list of arguments with `help(mols2grid.MolGrid.to_table)`

#### Output parameters

You can either:
* save the grid with `mols2grid.save(input, output="path/grid.html", ...)`. The file that is generated is a standalone HTML document that should work with most web browsers.
* display it directly in a Jupyter notebook with `mols2grid.display(...)` (optionnal argument: `width="100%"`, `height=None`)

## 🚀 Resources
---
* [Simple example](https://iwatobipen.wordpress.com/2021/06/13/draw-molecules-on-jupyter-notebook-rdkit-mols2grid/) by iwatobipen
* Creating a web app with Streamlit for filtering datasets:
  * [Blog post](https://blog.reverielabs.com/building-web-applications-from-python-scripts-with-streamlit/) by Justin Chavez
  * [Video tutorial](https://www.youtube.com/watch?v=0rqIwSeUImo) by Data Professor
* [Viewing clustered chemical structures](https://practicalcheminformatics.blogspot.com/2021/07/viewing-clustered-chemical-structures.html) and [Exploratory data analysis](https://practicalcheminformatics.blogspot.com/2021/10/exploratory-data-analysis-with.html) by Pat Walters
* [Advanced notebook (RDKit UGM 2021)](https://colab.research.google.com/github/rdkit/UGM_2021/blob/main/Notebooks/Bouysset_mols2grid.ipynb)

## 👏 Acknowledgments
---
* [@fredrikw](https://github.com/fredrikw) (contributor)
* [@JustinChavez](https://github.com/JustinChavez) (contributor)
* [@hadim](https://github.com/hadim) (conda feedstock maintainer)

## ⚖ License
---

Unless otherwise noted, all files in this directory and all subdirectories are distributed under the Apache License, Version 2.0:
```
    Copyright 2021-2022 Cédric BOUYSSET

    Licensed under the Apache License, Version 2.0 (the "License");
    you may not use this file except in compliance with the License.
    You may obtain a copy of the License at

        http://www.apache.org/licenses/LICENSE-2.0

    Unless required by applicable law or agreed to in writing, software
    distributed under the License is distributed on an "AS IS" BASIS,
    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and
    limitations under the License.
```
