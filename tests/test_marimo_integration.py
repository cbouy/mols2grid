import sys
from unittest.mock import MagicMock, patch

import pandas as pd
import pytest

from mols2grid import MolGrid
from mols2grid.utils import is_running_within_marimo
from mols2grid.widget import MolGridWidget


@pytest.fixture
def mock_marimo_module():
    with patch.dict(sys.modules, {"marimo": MagicMock()}):
        yield


@pytest.mark.usefixtures("mock_marimo_module")
def test_is_running_within_marimo_true():
    assert is_running_within_marimo() is True


def test_is_running_within_marimo_false():
    # Ensure marimo is not in sys.modules for this test
    with patch.dict(sys.modules):
        if "marimo" in sys.modules:
            del sys.modules["marimo"]
        assert is_running_within_marimo() is False


@pytest.fixture
def grid_fixture():
    df = pd.DataFrame({"SMILES": ["C"]})
    return df, MolGrid(df, smiles_col="SMILES")


@pytest.mark.usefixtures("mock_marimo_module")
def test_display_in_marimo(grid_fixture):
    _, mg = grid_fixture

    # Mock marimo.Html and marimo.vstack
    with patch("marimo.ui.anywidget") as mock_anywidget:
        result = mg.display()

        # Verify that the widget was wrapped with marimo's UI
        mock_anywidget.assert_called_once()
        args, _ = mock_anywidget.call_args
        widget = args[0]
        assert isinstance(widget, MolGridWidget)

        # Ensure the result is the return value of marimo.vstack
        assert result == mock_anywidget.return_value


@pytest.mark.usefixtures("mock_marimo_module")
def test_get_marimo_selection_before_rendering_raises(grid_fixture):
    _, mg = grid_fixture
    with pytest.raises(RuntimeError, match="run the `display` method first"):
        mg.get_marimo_selection()


@pytest.mark.usefixtures("mock_marimo_module")
def test_get_selection_state_inside_marimo(grid_fixture):
    _, mg = grid_fixture
    mg.render()

    # Mock marimo.state
    mock_get_state = MagicMock()
    mock_set_state = MagicMock()
    with patch(
        "marimo.state", return_value=(mock_get_state, mock_set_state)
    ) as mock_state:
        # Call get_marimo_selection
        state_getter = mg.get_marimo_selection()

        # Check if marimo.state was called with empty list
        mock_state.assert_called_once_with([])

        # Check if _marimo_hooked is set
        assert getattr(mg.widget, "_marimo_hooked", False) is True

        # Verify return value
        assert state_getter == mock_get_state


def test_get_selection_state_outside_marimo(grid_fixture):
    _, mg = grid_fixture

    # Ensure marimo is not in sys.modules
    with patch.dict(sys.modules):
        if "marimo" in sys.modules:
            del sys.modules["marimo"]

        with pytest.raises(RuntimeError, match="only available in a marimo notebook"):
            mg.get_marimo_selection()


@pytest.mark.usefixtures("mock_marimo_module")
def test_selection_state_update_logic(grid_fixture):
    _, mg = grid_fixture
    mg.render()

    mock_set_state = MagicMock()
    with (
        patch("marimo.state", return_value=(MagicMock(), mock_set_state)),
        patch.object(mg.widget, "observe") as mock_observe,
    ):
        # Inspect the observe call to capture the callback
        mg.get_marimo_selection()

        # Verify observe was called
        mock_observe.assert_called()
        args, _ = mock_observe.call_args
        callback = args[0]

        # Simulate event with valid selection
        # The widget returns a string representation of a dict
        new_selection = {1: "C", 2: "CC"}
        event = {"new": str(new_selection)}

        callback(event)
        mock_set_state.assert_called_with([1, 2])

        # Test invalid input (should pass silently)
        mock_set_state.reset_mock()
        callback({"new": "invalid json"})
        mock_set_state.assert_not_called()
