# mols2grid

mols2grid is a Python chemical viewer for 2D structures of small molecules, based on RDKit.

![Demo showing mols2grid's integration in a Jupyter notebook](demo.png)

## Installation

mols2grid was developped for Python 3.6+ and requires rdkit, pandas and jinja2 as dependencies.

To install mols2grid from a clean conda environment:
```shell
conda install -c conda-forge rdkit
pip install git+https://github.com/cbouy/mols2grid.git
```

## Usage
---

You can setup the grid from various inputs:
* a pandas **DataFrame** (with a SMILES column),
* a list of **RDKit molecules** (with properties accessible through the `mol.GetPropsAsDict()` method),
* or an **SDF file**:

```python
from mols2grid import MolGrid

# pandas
mg = MolGrid(df)

# List of molecules
mg = MolGrid.from_mols([Chem.MolFromSmiles(smi) for smi in smiles])

# SDF
mg = MolGrid.from_sdf("/path/to/molecules.sdf")
```

This first step sets up the internal dataframe that will be used by mols2grid (which you can access and sort/rename with `mg.dataframe`) and creates the drawings of each molecule.
Here, you can set up how the drawings should look:

  * `useSVG`: use SVG images (default) or PNG
  * `coordGen`: use the coordGen library (default) instead of the RDKit one to depict the molecules in 2D
  * `size`: size of the images (160 by 120)
  * and all the arguments available in RDKit's [MolDrawOptions](https://www.rdkit.org/docs/source/rdkit.Chem.Draw.rdMolDraw2D.html#rdkit.Chem.Draw.rdMolDraw2D.MolDrawOptions), like `addStereoAnnotation=True`

Two columns are automatically added to the `mg.dataframe` after this operation: `img` and `SMILES`. If a `SMILES` column already existed in your molecules' properties or SDF file, it will not be overwritten.

You can then call the `mg.display` or `mg.save` methods to render and display/save the HTML document.
  
You can control which template is used through the `template` argument:
  * `pages` (default) which is displayed above. It integrates nicely with notebooks
  * `table`, which displays the full list of molecules (no pages). Useful if you ever need to print the full list of molecules on paper (or print to PDF)

Both templates can be configured with the same parameters (a lot of which are [CSS](https://www.w3schools.com/cssref/) declarations):

* `subset`: list or None  
    Columns to be displayed in each cell of the grid. Each column's value will be displayed from top to bottom in the same order given here. Use `"img"` for the image of the molecule. Default: all columns (with "img" in first position)
* `tooltip`: list or None  
    Columns to be displayed as a tooltip when hovering/clicking on the image of a cell. Use `None` for no tooltip.
* `tooltip_fmt`: str  
    Format string of each key/value pair in the tooltip (CSS)
* `tooltip_trigger`: str  
    Sequence of triggers for the tooltip: (click, hover, focus)
* `tooltip_placement`: str  
    Position of the tooltip: auto, top, bottom, left, right
* `cell_width`: int  
    Max width of each cell, in pixels
* `n_cols`: int  
    Number of columns per page
* `border`: str  
    Styling of the border around each cell (CSS)
* `gap`: int or str  
    Size of the margin around each cell (CSS)
* `fontsize`: str  
    Font size of the text displayed in each cell (CSS)
* `fontfamily`: str  
    Font used for the text in each cell (CSS)
* `textalign`: str  
    Alignment of the text in each cell (CSS)
* `hover_color`: str  
    Background color when hovering a cell (CSS)
* `style`: dict or None  
    CSS styling applied to each item in a cell. The dict must follow a `key: function` structure where the key must correspond to one of the columns in `subset` or `tooltip`. The function takes the item's value as input, and outputs a valid CSS styling. For example, if you want to color the text corresponding to the "Solubility"
    column in your dataframe:
    ```python
    style={"Solubility": lambda x: "color: red" if x < -5 else "color: black"}
    ```

The `pages` template comes with additional parameters:

* `n_rows` : int  
    Number of rows per page
    

## License
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