from mols2grid import datafiles as datafiles
from mols2grid._version import __version__ as __version__
from mols2grid.callbacks import make_popup_callback as make_popup_callback
from mols2grid.dispatch import display as display
from mols2grid.dispatch import save as save
from mols2grid.molgrid import MolGrid as MolGrid
from mols2grid.select import get_selection as get_selection
from mols2grid.select import list_grids as list_grids
from mols2grid.utils import sdf_to_dataframe as sdf_to_dataframe

try:
    from google.colab import output
except ImportError:
    pass
else:
    output.enable_custom_widget_manager()
