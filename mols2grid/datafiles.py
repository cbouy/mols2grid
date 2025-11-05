from importlib import resources

datapath = resources.files("mols2grid").joinpath("data")

SOLUBILITY_SDF = str(datapath / "solubility.test.sdf")
