from functools import wraps
from importlib.util import find_spec
from rdkit import Chem
import pandas as pd

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
        v = transform[k](v) if transform.get(k) else v
        v = f'<span style="{style[k](v)}">{v}</span>' if style.get(k) else v
        items.append(fmt.format(key=k, value=v))
    return "<br>".join(items)

def mol_to_smiles(mol):
    """Returns a SMILES from an RDKit molecule, or None if not an RDKit mol"""
    return Chem.MolToSmiles(mol) if mol else None

def mol_to_record(mol, mol_col="mol"):
    """Function to create a dict of data from an RDKit molecule"""
    return {"SMILES": Chem.MolToSmiles(mol),
            **mol.GetPropsAsDict(includePrivate=True),
            mol_col: mol} if mol else {}

def sdf_to_dataframe(sdf_path, mol_col="mol"):
    """Returns a dataframe of molecules from an SDF file"""
    return pd.DataFrame([mol_to_record(mol, mol_col)
                         for mol in Chem.SDMolSupplier(sdf_path)])

def remove_coordinates(mol):
    """Removes the existing coordinates from the molecule. The molecule is
    modified inplace"""
    mol.RemoveAllConformers()
    return mol
