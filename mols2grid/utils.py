import re
from functools import wraps, partial
from importlib.util import find_spec
from jinja2 import Environment, FileSystemLoader
from pathlib import Path
from rdkit import Chem
import pandas as pd


env = Environment(loader=FileSystemLoader(Path(__file__).parent / 'templates'),
                  autoescape=False)

def requires(module):
    def inner(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            if find_spec(module):
                return func(*args, **kwargs)
            raise ModuleNotFoundError(
                f"The module {module!r} is required to use {func.__name__!r} "
                "but it is not installed!")
        return wrapper
    return inner

def tooltip_formatter(s, subset, fmt, style, transform):
    """Function to generate tooltips from a pandas Series

    Parameters
    ----------
    s : pandas.Series
        Row in the internal pandas DataFrame
    subset : list
        Subset of columns that are used for the tooltip
    fmt : str
        Format string for each key-value pair of the tooltip
    style : dict
        CSS styling applied to each item independently
    transform : dict
        Functions applied to each value before rendering
    """
    items = []
    for k, v in s[subset].to_dict().items():
        displayed = transform[k](v) if transform.get(k) else v
        v = (f'<span style="{style[k](v)}">{displayed}</span>'
             if style.get(k) else displayed)
        items.append(fmt.format(key=k, value=v))
    return "<br>".join(items)

def mol_to_smiles(mol):
    """Returns a SMILES from an RDKit molecule, or None if not an RDKit mol"""
    return Chem.MolToSmiles(mol) if mol else None

def mol_to_record(mol, mol_col="mol"):
    """Function to create a dict of data from an RDKit molecule"""
    return {**mol.GetPropsAsDict(includePrivate=True),
            mol_col: mol} if mol else {}

def sdf_to_dataframe(sdf_path, mol_col="mol"):
    """Returns a dataframe of molecules from an SDfile (.sdf or .sdf.gz)"""
    if sdf_path.endswith(".gz"):
        import gzip
        read_file = gzip.open
    else:
        read_file = partial(open, mode="rb")
    with read_file(sdf_path) as f:
        return pd.DataFrame([mol_to_record(mol, mol_col)
                             for mol in Chem.ForwardSDMolSupplier(f)])

def remove_coordinates(mol):
    """Removes the existing coordinates from the molecule. The molecule is
    modified inplace"""
    mol.RemoveAllConformers()
    return mol

def make_popup_callback(title, html, js="", style=""):
    """Creates a JavaScript callback that displays a popup window

    Parameters
    ----------
    title : str
        Title of the popup. Use `title='${data["Name"]}'` to use the value
        of the column "Name" as a title
    html : str
        Content of the popup window
    js : str
        JavaScript code executed before making the content of the popup window.
        This allows you to create variables and reuse them later in the `html`
        content of the popup, using the `${my_variable}` syntax
    style : str
        CSS style assigned to the popup window
    """
    return (env.get_template('js/popup.js')
               .render(js=js,
                       html=html,
                       title=title,
                       style=style))

def slugify(string):
    """Replaces whitespaces with hyphens"""
    return re.sub(r"\s+", "-", string)
