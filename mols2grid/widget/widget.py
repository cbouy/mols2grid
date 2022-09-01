from ipywidgets import DOMWidget, register
from traitlets import Unicode
from ._frontend import module_name, module_version


@register
class MolGridWidget(DOMWidget):
    """A custom widget for the MolGrid class. Handles selections and callbacks.

    Attributes
    ----------
    grid_id : str
        Name of the grid controlling the widget
    selection : str
        JSON string containing the molecule selection as a dictionnary. Index
        are keys and SMILES string are values.
    callback_kwargs : str
        JSON string containing the keyword arguments with which to call the
        callback function.
    """
    _model_name = Unicode('MolGridModel').tag(sync=True)
    _model_module = Unicode(module_name).tag(sync=True)
    _model_module_version = Unicode(module_version).tag(sync=True)
    _view_name = Unicode('MolGridView').tag(sync=True)
    _view_module = Unicode(module_name).tag(sync=True)
    _view_module_version = Unicode(module_version).tag(sync=True)

    grid_id = Unicode("default").tag(sync=True)
    selection = Unicode("{}").tag(sync=True)
    callback_kwargs = Unicode("{}").tag(sync=True)
