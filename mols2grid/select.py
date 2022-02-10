import warnings
from collections import UserDict
from threading import Lock

class SelectionRegister:
    """Register for grid selections

    Attributes
    ----------
    SELECTIONS : dict
        Stores each grid selection according to their name
    current_selection : str
        Name of the most recently updated grid
    """
    def __init__(self):
        self.SELECTIONS = {}
        self.lock = Lock()

    def _update_current_grid(self, name):
        self.current_selection = name

    def _init_grid(self, name):
        overwrite = self.SELECTIONS.get(name, False)
        if overwrite:
            warnings.warn(
                f"Overwriting non-empty {name!r} grid selection: {str(overwrite)}"
            )
        self.SELECTIONS[name] = {}
        self._update_current_grid(name)

    def _set_selection(self, _id, smiles):
        self.SELECTIONS[self.current_selection][_id] = smiles

    def _unset_selection(self, _id):
        del self.SELECTIONS[self.current_selection][_id]

    def add_selection(self, name, _id, smiles):
        """Add an entry to a grid

        Parameters
        ----------
        name : str
            Name of the grid to update
        _id : list[int]
            List of entry identifiers (`mols2grid-id`)
        smiles : list[str]
            List of SMILES

        Notes
        -----
        Manually using this function will not affect the checkbox displayed
        on the grid. It is only usefull when retrieving selections from Python
        with `mols2grid.get_selection()`.
        """
        with self.lock:
            self._update_current_grid(name)
            for i, s in zip(_id, smiles):
                self._set_selection(i, s)

    def del_selection(self, name, _id):
        """Remove an entry from a grid

         Parameters
        ----------
        name : str
            Name of the grid to update
        _id : list[int]
            List of entry identifiers (`mols2grid-id`)

        Notes
        -----
        Manually using this function will not affect the checkbox displayed
        on the grid. It is only usefull when retrieving selections from Python
        with `mols2grid.get_selection()`.
        """
        with self.lock:
            self._update_current_grid(name)
            for i in _id:
                self._unset_selection(i)

    def get_selection(self, name=None):
        """Returns the selection for a specific MolGrid instance

        Parameters
        ----------
        name : str or None
            Name of the grid to fetch the selection from. If `None`, the most
            recently updated grid is returned
        """
        name = self.current_selection if name is None else name
        return self.SELECTIONS[name]

    def list_grids(self):
        """Returns a list of grid names"""
        return list(self.SELECTIONS.keys())

    def _clear(self):
        """Clears all selections"""
        if hasattr(self, "current_selection"):
            del self.current_selection
        self.SELECTIONS.clear()


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
    colab.output.register_callback('_m2g_reg.add_selection',
                                   register.add_selection)
    colab.output.register_callback('_m2g_reg.del_selection',
                                   register.del_selection)
