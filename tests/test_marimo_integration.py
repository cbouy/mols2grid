import sys
from types import SimpleNamespace
from unittest.mock import MagicMock, patch

import pandas as pd
import pytest

from mols2grid import MolGrid
from mols2grid.select import register
from mols2grid.utils import is_running_within_marimo
from mols2grid.widget import MolGridWidget


@pytest.fixture
def mock_marimo_module(monkeypatch: pytest.MonkeyPatch):
    monkeypatch.setattr(sys, "modules", {**sys.modules, "marimo": MagicMock()})


@pytest.mark.usefixtures("mock_marimo_module")
def test_is_running_within_marimo_true():
    assert is_running_within_marimo() is True


def test_is_running_within_marimo_false():
    if "marimo" in sys.modules:
        assert isinstance(sys.modules["marimo"], MagicMock)
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
        state_getter, _ = register.link_marimo_state()

        # Check if marimo.state was called with empty list
        mock_state.assert_called_once_with({})

        # Verify return value
        assert state_getter == mock_get_state


def test_get_selection_state_outside_marimo():
    # Ensure marimo is not in sys.modules
    with patch.dict(sys.modules):
        if "marimo" in sys.modules:
            del sys.modules["marimo"]

        with pytest.raises(RuntimeError, match="only available in a marimo notebook"):
            register.link_marimo_state()


@pytest.mark.usefixtures("mock_marimo_module")
def test_selection_state_update_logic(grid_fixture, monkeypatch: pytest.MonkeyPatch):
    _, mg = grid_fixture
    mock_set_state = MagicMock()
    with (
        patch("marimo.state", return_value=(MagicMock(), mock_set_state)),
        patch.object(MolGridWidget, "observe") as mock_observe,
    ):
        # Inspect the observe call to capture the callback
        register.link_marimo_state()
        mock_callback = MagicMock(wraps=register.CALLBACKS[-1])
        monkeypatch.setattr(register, "CALLBACKS", [mock_callback])

        mg.render()
        # Verify observe was called
        mock_observe.assert_called()

        # Simulate event with valid selection
        event_values = {1: "C", 2: "CC"}
        event = SimpleNamespace(new=str(event_values))
        register.selection_updated("default", event)

        # check callback was called with expected selection
        mock_callback.assert_called_with("default", event_values)

        # check inner lambda works as expected: given a current state with
        # mol 42 selected, and new state where only 1 and 2 are selected,
        # 42 should disappear and only 1, 2 remain
        lambda_setter = mock_set_state.call_args[0][0]
        result = lambda_setter({"default": [42]})
        assert result == {"default": [1, 2]}
