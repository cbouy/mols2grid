import gzip
import re
from ast import literal_eval
from functools import partial, wraps
from importlib.util import find_spec
from pathlib import Path

import pandas as pd
from jinja2 import Environment, FileSystemLoader
from rdkit import Chem

env = Environment(
    loader=FileSystemLoader(Path(__file__).parent / "templates"), autoescape=False
)


def requires(module):
    def inner(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            if find_spec(module):
                return func(*args, **kwargs)
            raise ModuleNotFoundError(
                f"The module {module!r} is required to use {func.__name__!r} "
                "but it is not installed!"
            )

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
        v = (
            f'<span class="copy-me" style="{style[k](v)}">{displayed}</span>'
            if style.get(k)
            else f'<span class="copy-me">{displayed}</span>'
        )
        items.append(fmt.format(key=k, value=v))
    return '<br>'.join(items)


def mol_to_smiles(mol):
    """Returns a SMILES from an RDKit molecule, or None if not an RDKit mol"""
    return Chem.MolToSmiles(mol) if mol else None


def mol_to_record(mol, mol_col="mol"):
    """Function to create a dict of data from an RDKit molecule"""
    return {**mol.GetPropsAsDict(includePrivate=True), mol_col: mol} if mol else {}


def sdf_to_dataframe(sdf_path, mol_col="mol"):
    """Creates a dataframe of molecules from an SDFile. All property fields in
    the SDFile are made available in the resulting dataframe

    Parameters
    ----------
    sdf_path : str, Path
        Path to the SDFile, ending with either ``.sdf`` or ``.sdf.gz``
    mol_col : str
        Name of the column containing the RDKit molecules in the dataframe

    Returns
    -------
    df : pandas.DataFrame
    """
    if str(sdf_path).endswith(".gz"):
        read_file = gzip.open
    else:
        read_file = partial(open, mode="rb")
    with read_file(sdf_path) as f:
        return pd.DataFrame(
            [mol_to_record(mol, mol_col)
             for mol in Chem.ForwardSDMolSupplier(f)]
        )


def remove_coordinates(mol):
    """Removes the existing coordinates from the molecule. The molecule is
    modified inplace"""
    mol.RemoveAllConformers()
    return mol


def slugify(string):
    """Replaces whitespaces with hyphens"""
    return re.sub(r"\s+", "-", string)


def callback_handler(callback, event):
    """Handler for applying the callback function on change"""
    data = literal_eval(event.new)
    callback(data)


def _get_streamlit_script_run_ctx():
    from streamlit.runtime.scriptrunner import get_script_run_ctx

    return get_script_run_ctx()


def is_running_within_streamlit():
    """
    Function to check whether python code is run within streamlit

    Returns
    -------
    use_streamlit : boolean
        True if code is run within streamlit, else False
    """
    try:
        ctx = _get_streamlit_script_run_ctx()
    except ImportError:
        return False
    else:
        return ctx is not None
