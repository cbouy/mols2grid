from threading import Thread
import pytest
import mols2grid as mg
from mols2grid.select import register

data = [{"SMILES": "C"*i, "ID": i} for i in range(1, 5)]

@pytest.fixture(autouse=True)
def clear_register_before_each_test():
    register._clear()
    register._init_grid("foo")
    yield

def test_clear_register():
    register._clear()
    assert not hasattr(register, "current_selection")
    assert register.SELECTIONS == {}

def test_update_current_grid():
    register._clear()
    mg.MolGrid(data, name="foo")
    assert register.current_selection == "foo"

def test_init_grid():
    assert "foo" in register.SELECTIONS.keys()
    assert register.current_selection == "foo"

def test_overwrite_warning():
    register._set_selection(0, "C")
    with pytest.warns(
        UserWarning,
        match=f"Overwriting non-empty 'foo' grid selection: {{0: 'C'}}"
    ):
        register._init_grid("foo")
    assert register.get_selection() == {}

def test_set_selection():
    for index, smi in [(0, "C"), (1, "CC")]:
        register._set_selection(index, smi)
        assert register.SELECTIONS["foo"][index] == smi

def test_unset_selection():
    for index, smi in [(0, "C"), (1, "CC")]:
        register._set_selection(index, smi)
    register._unset_selection(0)
    assert 0 not in register.SELECTIONS["foo"].keys()
    assert "C" not in register.SELECTIONS["foo"].values()

def test_add_selection():
    register._init_grid("bar")
    register.add_selection("foo", [0, 1], ["C", "CC"])
    assert register.get_selection()[0] == "C"
    assert register.get_selection()[1] == "CC"

def test_del_selection():
    register._init_grid("bar")
    register.add_selection("foo", [0, 1], ["C", "CC"])
    register.del_selection("foo", [1])
    assert register.get_selection() == {0: "C"}

def test_get_selection():
    assert register.get_selection() == {}
    d = {0: "C", 1: "CC"}
    for index, smi in d.items():
        register._set_selection(index, smi)
    assert register.get_selection() == d

def test_list_grids():
    assert register.list_grids() == ["foo"]
    register._init_grid("bar")
    assert register.list_grids() == ["foo", "bar"]
    register._init_grid("foo")
    assert register.list_grids() == ["foo", "bar"]

def test_warn_selection():
    with pytest.warns(UserWarning,
                      match="`mols2grid.selection` is deprecated"):
        len(mg.selection)

def test_lock(n=1e6):
    N = int(n)
    x = list(range(N))
    def func(name):
        register.add_selection(name, x, x)
    register._init_grid("bar")
    t1 = Thread(target=func, args=("foo",))
    t2 = Thread(target=func, args=("bar",))
    t1.start()
    t2.start()
    t1.join()
    t2.join()
    assert len(register.get_selection("foo")) == N
    assert len(register.get_selection("bar")) == N
