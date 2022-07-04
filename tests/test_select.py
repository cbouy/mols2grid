from types import SimpleNamespace
import pytest
import mols2grid as mg
from mols2grid.select import register

data = [{"SMILES": "C"*i, "ID": i} for i in range(1, 5)]

@pytest.fixture(autouse=True)
def clear_register_between_tests():
    register._clear()
    register._init_grid("foo")
    yield
    register._clear()

def test_clear_register():
    register._clear()
    assert not hasattr(register, "current_selection")
    assert register.SELECTIONS == {}

def test_update_current_grid():
    mg.MolGrid(data, name="bar")
    assert register.current_selection == "bar"

def test_init_grid():
    assert "foo" in register.SELECTIONS.keys()
    assert register.current_selection == "foo"

def test_overwrite_warning():
    event = SimpleNamespace(new='{0: "C"}')
    register.selection_updated("foo", event)
    with pytest.warns(
        UserWarning,
        match=f"Overwriting non-empty 'foo' grid selection: {{0: 'C'}}"
    ):
        register._init_grid("foo")
    assert register.get_selection() == {}

def test_update_and_get_selection():
    assert register.get_selection() == {}
    event = SimpleNamespace(new='{0: "CCO"}')
    register.selection_updated("foo", event)
    assert register.get_selection() == {0: "CCO"}
    event.new = '{}'
    register.selection_updated("foo", event)
    assert register.get_selection() == {}

def test_list_grids():
    assert register.list_grids() == ["foo"]
    register._init_grid("bar")
    assert register.list_grids() == ["foo", "bar"]
    register._init_grid("foo")
    assert register.list_grids() == ["foo", "bar"]
