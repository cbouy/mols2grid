import warnings
from ast import literal_eval

from mols2grid.utils import is_running_within_marimo


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
        self.CALLBACKS = []
        self._current_selection = None

    @property
    def current_selection(self):
        """The name of the last updated grid (created or interacted with),
        or ``None`` if no grid has been displayed yet."""
        return self._current_selection

    @current_selection.setter
    def current_selection(self, name):
        if name is not None and name not in self.SELECTIONS:
            raise ValueError(
                f"The selection for {name} must be initialized "
                "before setting it as the current grid"
            )
        self._current_selection = name

    @current_selection.deleter
    def current_selection(self):
        self._current_selection = None

    def add_grid(self, name):
        """Adds a grid to track selections for."""
        overwrite = self.SELECTIONS.get(name, False)
        if overwrite and not is_running_within_marimo():
            warnings.warn(
                f"Overwriting non-empty {name!r} grid selection: {overwrite!s}",
                stacklevel=2,
            )
        self.SELECTIONS[name] = {}
        self.current_selection = name
        self.add_callback(self._store_selection)

    def selection_updated(self, name, event):
        """Callback function linked to the widget."""
        self.current_selection = name
        selection = literal_eval(event.new)
        for callback in self.CALLBACKS:
            callback(name, selection)

    def _store_selection(self, name, selection):
        """Makes the selection available to the register."""
        self.SELECTIONS[name] = selection

    def add_callback(self, callback):
        """Add a callback function to be called when the selection is updated.

        Parameters
        ----------
        callback : callable
            The function to execute when a the selection is updated.
        """
        self.CALLBACKS.append(callback)

    def get_selection(self, name=None):
        """Returns the selection for a specific MolGrid instance

        Parameters
        ----------
        name : str or None
            Name of the grid to fetch the selection from. If ``None``, the most
            recently updated grid is returned
        """
        name = self.current_selection if name is None else name
        return self.SELECTIONS[name]

    def link_marimo_state(self):
        """Link the register to marimo by initializing a ``state`` dict.
        When the selection on a grid is updated, the state setter function from marimo
        is called to update the state's value.

        Returns
        -------
        get_state, set_state:
            The state getter and setter returned by ``marimo.state({})``.
            ``get_state`` returns a dictionary of all selections where the keys are
            grid names and values are the indices selected for that grid.
        """
        if not is_running_within_marimo():
            raise RuntimeError("This method is only available in a marimo notebook.")

        import marimo as mo

        get_state, set_state = mo.state({})

        def marimo_callback(name, selection):
            set_state(lambda value: {**value, name: list(selection)})

        self.add_callback(marimo_callback)
        return get_state, set_state

    def list_grids(self):
        """Returns a list of grid names"""
        return list(self.SELECTIONS)

    def clear(self):
        """Clears all selections"""
        self.SELECTIONS.clear()
        self.CALLBACKS.clear()
        del self.current_selection


register = SelectionRegister()
get_selection = register.get_selection
list_grids = register.list_grids
link_marimo_state = register.link_marimo_state
