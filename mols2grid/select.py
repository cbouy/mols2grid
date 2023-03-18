import warnings
from ast import literal_eval


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

    def selection_updated(self, name, event):
        self.SELECTIONS[name] = literal_eval(event.new)
        self._update_current_grid(name)

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


register = SelectionRegister()
get_selection = register.get_selection
list_grids = register.list_grids
