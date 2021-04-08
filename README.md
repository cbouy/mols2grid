# 👀 mols2grid

[![Pypi Version](https://img.shields.io/pypi/v/mols2grid.svg)](https://pypi.python.org/pypi/mols2grid)
[![Build status](https://github.com/cbouy/mols2grid/workflows/build/badge.svg)](https://github.com/cbouy/mols2grid/actions/workflows/build.yml)
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/cbouy/mols2grid/blob/master/demo.ipynb)

mols2grid is an interactive chemical viewer for 2D structures of small molecules, based on RDKit.

![Demo showing mols2grid's integration in a Jupyter notebook](https://raw.githubusercontent.com/cbouy/mols2grid/master/demo.png)

➡️ Try the demo notebook on [Google Colab](https://colab.research.google.com/github/cbouy/mols2grid/blob/master/demo.ipynb)

## 🐍 Installation
---

mols2grid was developped for Python 3.6+ and requires rdkit (>=2019.09.1), pandas and jinja2 as dependencies.

To install mols2grid from a clean conda environment:
```shell
conda install -c conda-forge 'rdkit>=2019.09.1'
pip install mols2grid
```

It is compatible with Jupyter Notebook and Google Colab (Visual Studio notebooks and Jupyterlab are not supported) and can run on Streamlit.

## 📜 Usage
---

```python
import mols2grid

mols2grid.display("path/to/molecules.sdf",
                  # RDKit's MolDrawOptions parameters
                  fixedBondLength=25,
                  # rename fields for the output document
                  mapping={"SOL": "Solubility",
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
* or an **SDF file**

You can also rename each field of your input with the `mapping` parameter. Please note that 3 fields are automatically added regardless of your input: `mols2grid-id`, `SMILES` and `img`. If a "SMILES" field already exists, it will not be overwritten.

#### Parameters for the drawing of each molecule

* `useSVG=True`: use SVG images or PNG
* `coordGen=True`: use the coordGen library instead of the RDKit one to depict the molecules in 2D
* `size=(160, 120)`: size of each image
* `use_coords=True`: use the coordinates of the input molecules if available
* `MolDrawOptions=None`: RDKit's [MolDrawOptions](https://www.rdkit.org/docs/source/rdkit.Chem.Draw.rdMolDraw2D.html#rdkit.Chem.Draw.rdMolDraw2D.MolDrawOptions) class. Useful for making highly customized drawings. You can also leave this to `None`, and directly use the attributes of this class as parameters like `addStereoAnnotation=True`

#### Parameters for the grid
  
You can control the general look of the document through the `template` argument:
* `template="pages"` (default) which is displayed above. It integrates nicely with Jupyter notebooks and has a search bar
* `template="table"`, which displays the full list of molecules (no pages). Useful if you ever need to print the full list of molecules on paper (or print to PDF)

Both templates can be configured with the same parameters (a lot of which are [CSS](https://www.w3schools.com/cssref/) declarations):

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
* `border="1px solid #cccccc"`: str  
    Styling of the border around each cell (CSS)
* `gap=0`: int or str  
    Size of the margin around each cell (CSS)
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
    style={"Solubility": lambda x: "color: red" if x < -3 else "color: black"}
    ```
* `transform=None`: dict or None
    Functions applied to specific items in all cells. The dict must follow a `key: function` structure where the key must correspond to one of the columns in `subset`. The function takes the item's value as input and transforms it. For example, to round the "Solubility" to 2 decimals, and display the "Melting point" in Celsius instead of Fahrenheit with a single digit precision and some text before ("MP") and after ("°C") the value:
    ```python
    transform={"Solubility": lambda x: f"{x:.2f}",
               "Melting point": lambda x: f"MP: {5/9*(x-32):.1f}°C"}
    ```
    These transformations only affect columns in `subset` (not `tooltip`) and are applied independantly from `style`.

The `pages` template comes with additional parameters:

* `n_rows=3` : int  
    Number of rows per page
* `selection=True` : bool
    Enables the selection of molecules using a checkbox. Only usefull in the context of a Jupyter notebook. You can retrieve your selection of molecules (index and SMILES) through `mols2grid.selection`

The `pages` template also allows searching (by text or SMARTS) and sorting the grid.

#### Output parameters

You can either:
* save the grid with `mols2grid.save(output_path, ...)`. The file that is generated is a standalone HTML document that should work with most web browsers.
* display it directly in a Jupyter notebook with `mols2grid.display(...)` (optionnal argument: `width="100%"`, `height=None`)


## ⚖ License
---

Unless otherwise noted, all files in this directory and all subdirectories are distributed under the Apache License, Version 2.0:
```
    Copyright 2021 Cédric BOUYSSET

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
