from ipywidgets import register, DOMWidget
from traitlets import Unicode

@register
class CommWidget(DOMWidget):
    """Custom widget to handle comm between the grid and the Python backend"""
    _view_name = Unicode('CommWidgetView').tag(sync=True)
    _model_name = Unicode('CommWidgetModel').tag(sync=True)
    _view_module = Unicode('mols2grid_widget').tag(sync=True)
    _model_module = Unicode('mols2grid_widget').tag(sync=True)
    _view_module_version = Unicode('^0.1.0').tag(sync=True)
    _model_module_version = Unicode('^0.1.0').tag(sync=True)

    grid_id = Unicode("default").tag(sync=True)
    selection = Unicode("{}").tag(sync=True)
    callback_args = Unicode("{}").tag(sync=True)
