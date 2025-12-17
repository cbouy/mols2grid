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
