from types import SimpleNamespace
from typing import cast
from unittest.mock import Mock, patch

import pandas as pd
import pytest

from mols2grid.utils import (
    callback_handler,
    is_running_within_streamlit,
    requires,
    slugify,
    tooltip_formatter,
)


def test_requires():
    @requires("_not_a_module")
    def func():
        pass

    with pytest.raises(
        ModuleNotFoundError, match="The module '_not_a_module' is required"
    ):
        func()


@pytest.mark.parametrize(
    ("subset", "fmt", "style", "transform", "exp"),
    [
        (
            ["SMILES", "ID"],
            "<strong>{key}</strong>: {value}",
            {},
            {},
            '<strong>SMILES</strong>: <span class="copy-me">CCO</span>'
            '<strong>ID</strong>: <span class="copy-me">0</span>',
        ),
        (["ID"], "foo-{value}", {}, {}, 'foo-<span class="copy-me">0</span>'),
        (
            ["ID"],
            "{value}",
            {"ID": lambda x: "color: red"},  # noqa: ARG005
            {},
            '<span class="copy-me" style="color: red">0</span>',
        ),
        (
            ["Activity"],
            "{value}",
            {},
            {"Activity": lambda x: f"{x:.2f}"},
            '<span class="copy-me">42.01</span>',
        ),
        (
            ["Activity"],
            "{key}: {value}",
            {"Activity": lambda x: "color: red" if x > 40 else ""},
            {"Activity": lambda x: f"{x:.2f}"},
            'Activity: <span class="copy-me" style="color: red">42.01</span>',
        ),
    ],
)
def test_tooltip_formatter(subset, fmt, style, transform, exp):
    row = pd.Series(
        {
            "ID": 0,
            "SMILES": "CCO",
            "Activity": 42.012345,
        }
    )
    tooltip = tooltip_formatter(row, subset, fmt, style, transform)
    assert tooltip.startswith(exp)


@pytest.mark.parametrize(
    ("string", "expected"),
    [
        ("Mol", "Mol"),
        ("mol name", "mol-name"),
        ("mol  name", "mol-name"),
        ("mol-name", "mol-name"),
        ("mol- name", "mol--name"),
        ("mol\tname", "mol-name"),
        ("mol\nname", "mol-name"),
        ("mol \t\n name", "mol-name"),
    ],
)
def test_slugify(string, expected):
    assert slugify(string) == expected


@pytest.mark.parametrize("value", [1, 2])
def test_callback_handler(value):
    def callback(x):
        return x + 1

    mock = Mock(side_effect=callback)
    event = SimpleNamespace(new=str(value))
    callback_handler(mock, event)
    mock.assert_called_once_with(value)


def test_is_running_within_streamlit():
    mock_streamlit = Mock(runtime=Mock(scriptrunner=Mock(get_script_run_ctx=Mock())))
    mocked_module = cast(Mock, mock_streamlit.runtime.scriptrunner)
    mocked_func = cast(Mock, mocked_module.get_script_run_ctx)
    with patch.dict("sys.modules", {"streamlit.runtime.scriptrunner": mocked_module}):
        mocked_func.side_effect = ImportError()
        assert not is_running_within_streamlit()
        mocked_func.side_effect = None
        mocked_func.return_value = object()
        assert is_running_within_streamlit()
        mocked_func.return_value = None
        assert not is_running_within_streamlit()
