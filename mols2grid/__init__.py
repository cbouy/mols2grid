from ._version import __version__
from .callbacks import make_popup_callback
from .dispatch import display, save
from .molgrid import MolGrid
from .select import get_selection, list_grids
from .utils import sdf_to_dataframe
from .widget import _jupyter_labextension_paths, _jupyter_nbextension_paths

try:
    from google.colab import output
except ImportError:
    pass
else:
    output.enable_custom_widget_manager()
