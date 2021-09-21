import warnings
from collections import UserDict

class SelectionRegister:
    """Register for grid selections"""
    def __init__(self):
        self.SELECTIONS = {}
        self.current_selection = 0

    def _update_current_grid(self, grid_id):
        self.current_selection = grid_id

    def _set_selection(self, _id, smiles):
        self.SELECTIONS[self.current_selection][_id] = smiles

    def _del_selection(self, _id):
        del self.SELECTIONS[self.current_selection][_id]

    def get_selection(self, grid_id=None):
        """Returns the selection for a specific MolGrid instance

        Parameters
        ----------
        grid_id : int, str or None
            Name of the grid to fetch the selection from. If `None`, the most
            recently updated grid is returned
        """
        grid_id = self.current_selection if grid_id is None else grid_id
        return self.SELECTIONS.get(grid_id, {})

    def list_grids(self):
        """Returns a list of grid names"""
        return list(self.SELECTIONS.keys())


# deprecate old selection system
def warn(func):
    def wrapper(*args, **kwargs):
        warnings.warn(
            "Accessing the current grid's selection through "
            "`mols2grid.selection` is deprecated and will be removed soon. "
            "Please use `mols2grid.get_selection()` instead")
        return func(*args, **kwargs)
    return wrapper


class _OldSelection(UserDict):
    def __init__(self, register):
        self._sel = register.get_selection
        super().__init__()

    @warn
    def __repr__(self):
        return repr(self._sel())

    @warn
    def __iter__(self):
        return iter(self._sel())

    @warn
    def __len__(self):
        return len(self._sel())

    @warn
    def __missing__(self, key):
        return self._sel()[key]


register = SelectionRegister()
get_selection = register.get_selection
list_grids = register.list_grids
selection = _OldSelection(register)


# Callbacks for Google Colab
try:
    from google import colab
except (ModuleNotFoundError, ImportError):
    pass
else:
    colab.output.register_callback('m2g_sel._update_current_grid',
                                   register._update_current_grid)
    colab.output.register_callback('m2g_sel._set_selection',
                                   register._set_selection)
    colab.output.register_callback('m2g_sel._del_selection',
                                   register._del_selection)
