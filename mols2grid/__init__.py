from mols2grid import datafiles
from mols2grid._version import __version__
from mols2grid.callbacks import make_popup_callback
from mols2grid.dispatch import display, save
from mols2grid.molgrid import MolGrid
from mols2grid.select import get_selection, list_grids
from mols2grid.utils import sdf_to_dataframe

try:
    from google.colab import output
except ImportError:
    pass
else:
    output.enable_custom_widget_manager()
