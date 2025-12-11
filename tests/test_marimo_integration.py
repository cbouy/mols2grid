
import sys
import pytest
from unittest.mock import patch, MagicMock
import pandas as pd
from mols2grid import MolGrid
from mols2grid.utils import is_running_within_marimo

@pytest.fixture
def mock_marimo_module():
    with patch.dict(sys.modules, {"marimo": MagicMock()}):
        yield

def test_is_running_within_marimo_true(mock_marimo_module):
    assert is_running_within_marimo() is True

def test_is_running_within_marimo_false():
    # Ensure marimo is not in sys.modules for this test
    with patch.dict(sys.modules):
        if "marimo" in sys.modules:
            del sys.modules["marimo"]
        assert is_running_within_marimo() is False

def test_display_in_marimo(mock_marimo_module):
    df = pd.DataFrame({"SMILES": ["C"]})
    mg = MolGrid(df, smiles_col="SMILES")
    
    # Mock marimo.Html
    with patch("marimo.Html") as mock_html:
        result = mg.display()
        mock_html.assert_called_once()
        
        # Verify that an iframe is being rendered
        args, _ = mock_html.call_args
        html_content = args[0]
        assert "<iframe" in html_content
        assert 'class="mols2grid-iframe"' in html_content
        
        # Ensure the result is the return value of marimo.Html
        assert result == mock_html.return_value

