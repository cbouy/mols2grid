# ruff: noqa: F401
from mols2grid import datafiles
from mols2grid._version import __version__
from mols2grid.callbacks import make_popup_callback
from mols2grid.dispatch import display, save
from mols2grid.molgrid import MolGrid
from mols2grid.select import get_selection, list_grids
from mols2grid.utils import is_running_within_streamlit, sdf_to_dataframe

try:
    from google.colab import output
except (ImportError, ModuleNotFoundError):
    pass
else:
    output.enable_custom_widget_manager()
    del output

if is_running_within_streamlit():
    import os

    os.environ["M2G_DEBUG"] = "1"

del is_running_within_streamlit
