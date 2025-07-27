from pathlib import Path

import anywidget
from traitlets import Bool, List, Unicode

from mols2grid._version import __version__  # noqa: F401

BUNDLER_OUTPUT_DIR = Path(__file__).parent / "static"


class MolGridWidget(anywidget.AnyWidget):
    """A custom widget for the MolGrid class. Handles selections and callbacks.

    Attributes
    ----------
    selection : str
        JSON string containing the molecule selection as a dictionnary. Index
        are keys and SMILES string are values.
    callback_kwargs : str
        JSON string containing the keyword arguments with which to call the
        callback function.
    filter_mask : List[bool]
        List stating wether a molecule should be kept (True) or filtered out
        (False)
    """

    _esm = BUNDLER_OUTPUT_DIR / "widget.js"
    _css = BUNDLER_OUTPUT_DIR / "widget.css"

    options = Unicode().tag(sync=True)
    selection = Unicode("{}").tag(sync=True)
    callback_kwargs = Unicode("{}").tag(sync=True)
    filter_mask = List(Bool(), []).tag(sync=True)
