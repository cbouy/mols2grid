from copy import deepcopy

selection = {}
current_selection = 0

def _update_current_grid(grid_id):
    global current_selection
    current_selection = grid_id

def _set_selection(_id, smiles):
    global selection
    global current_selection
    selection[current_selection][_id] = smiles

def _del_selection(_id):
    global selection
    global current_selection
    del selection[current_selection][_id]

def get_selection(grid_id=None):
    """Returns the selection for a specific MolGrid instance. If grid_id is
    `None`, the most recently updated grid is returned."""
    global selection
    global current_selection
    grid_id = current_selection if grid_id is None else grid_id
    sel = selection[grid_id]
    return deepcopy(sel)

try:
    from google import colab
except (ModuleNotFoundError, ImportError):
    pass
else:
    colab.output.register_callback('m2g._update_current_grid',
                                   _update_current_grid)
    colab.output.register_callback('m2g._set_selection', _set_selection)
    colab.output.register_callback('m2g._del_selection', _del_selection)
