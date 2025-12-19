import sys
from unittest.mock import MagicMock, patch

import pandas as pd
import pytest

from mols2grid import MolGrid
from mols2grid.utils import is_running_within_marimo


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


@pytest.mark.usefixtures("mock_marimo_module")
def test_init_in_marimo_does_not_display():
    df = pd.DataFrame({"SMILES": ["C"]})

    # Mock IPython.display.display which is imported as display in molgrid.py
    # We need to patch it where it is used, i.e., in mols2grid.molgrid
    with patch("mols2grid.molgrid.display") as mock_display:
        _ = MolGrid(df, smiles_col="SMILES")
        mock_display.assert_not_called()


@pytest.mark.usefixtures("mock_marimo_module")
def test_display_in_marimo():
    df = pd.DataFrame({"SMILES": ["C"]})
    mg = MolGrid(df, smiles_col="SMILES")

    # Mock marimo.Html and marimo.vstack
    with patch("marimo.Html") as mock_html, patch("marimo.vstack") as mock_vstack:
        result = mg.display()

        # Verify that an iframe is being rendered inside Html
        mock_html.assert_called_once()
        args, _ = mock_html.call_args
        html_content = args[0]
        assert "<iframe" in html_content
        assert 'class="mols2grid-iframe"' in html_content

        # Verify vstack was called with [widget, html]
        mock_vstack.assert_called_once()
        vstack_args = mock_vstack.call_args[0][0]
        assert len(vstack_args) == 2
        assert vstack_args[0] == mg.widget
        assert vstack_args[1] == mock_html.return_value

        # Ensure the result is the return value of marimo.vstack
        assert result == mock_vstack.return_value


@pytest.mark.usefixtures("mock_marimo_module")
def test_get_selection_state_inside_marimo():
    df = pd.DataFrame({"SMILES": ["C"]})
    mg = MolGrid(df, smiles_col="SMILES")

    # Mock marimo.state
    mock_get_state = MagicMock()
    mock_set_state = MagicMock()
    with patch("marimo.state", return_value=(mock_get_state, mock_set_state)) as mock_state:
        # Call get_selection_state
        state_getter = mg.get_selection_state()

        # Check if marimo.state was called with empty list
        mock_state.assert_called_once_with([])

        # Check if _marimo_hooked is set
        assert getattr(mg.widget, "_marimo_hooked", False) is True

        # Verify return value
        assert state_getter == mock_get_state


def test_get_selection_state_outside_marimo():
    df = pd.DataFrame({"SMILES": ["C"]})
    mg = MolGrid(df, smiles_col="SMILES")

    # Ensure marimo is not in sys.modules
    with patch.dict(sys.modules):
        if "marimo" in sys.modules:
            del sys.modules["marimo"]

        with pytest.raises(RuntimeError, match="only available in a marimo notebook"):
            mg.get_selection_state()


@pytest.mark.usefixtures("mock_marimo_module")
def test_selection_state_update_logic():
    df = pd.DataFrame({"SMILES": ["C"]})
    mg = MolGrid(df, smiles_col="SMILES")

    mock_set_state = MagicMock()
    with patch("marimo.state", return_value=(MagicMock(), mock_set_state)):
        # Inspect the observe call to capture the callback
        with patch.object(mg.widget, "observe") as mock_observe:
            mg.get_selection_state()

            # Verify observe was called
            assert mock_observe.called
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
