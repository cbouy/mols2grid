import atexit
from contextlib import ExitStack
from importlib import resources

_file_manager = ExitStack()
atexit.register(_file_manager.close)
_data_resource = resources.files("mols2grid") / "data/"
datapath = _file_manager.enter_context(resources.as_file(_data_resource))

SOLUBILITY_SDF = str(datapath / "solubility.test.sdf")
