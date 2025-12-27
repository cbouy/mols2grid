from types import SimpleNamespace

import pytest

import mols2grid as mg
from mols2grid.select import register


@pytest.fixture(autouse=True)
def clear_register_between_tests():
    register.clear()
    register.add_grid("foo")
    yield
    register.clear()


def test_clear_register():
    register.clear()
    assert register.current_selection is None
    assert register.SELECTIONS == {}


def test_update_current_grid(smiles_records):
    mg.MolGrid(smiles_records, name="bar").display()
    assert register.current_selection == "bar"


def test_add_grid():
    assert "foo" in register.SELECTIONS
    assert register.current_selection == "foo"


def test_overwrite_warning():
    event = SimpleNamespace(new='{0: "C"}')
    register.selection_updated("foo", event)
    with pytest.warns(
        UserWarning, match="Overwriting non-empty 'foo' grid selection: {0: 'C'}"
    ):
        register.add_grid("foo")
    assert register.get_selection() == {}


def test_update_and_get_selection():
    assert register.get_selection() == {}
    event = SimpleNamespace(new='{0: "CCO"}')
    register.selection_updated("foo", event)
    assert register.get_selection() == {0: "CCO"}
    event.new = "{}"
    register.selection_updated("foo", event)
    assert register.get_selection() == {}


def test_list_grids():
    assert register.list_grids() == ["foo"]
    register.add_grid("bar")
    assert register.list_grids() == ["foo", "bar"]
    register.add_grid("foo")
    assert register.list_grids() == ["foo", "bar"]
