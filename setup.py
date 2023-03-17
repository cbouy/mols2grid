import os
from os.path import join as pjoin

from jupyter_packaging import (
    combine_commands,
    create_cmdclass,
    ensure_targets,
    get_version,
    install_npm,
    skip_if_exists,
)
from setuptools import setup

HERE = os.path.dirname(os.path.abspath(__file__))

# The name of the project
name = "mols2grid"

# Get the version
version = get_version(pjoin(HERE, name, "_version.py"))

# Representative files that should exist after a successful build
jstargets = [
    pjoin(HERE, name, "nbextension", "index.js"),
    pjoin(HERE, name, "labextension", "package.json"),
]

package_data_spec = {name: ["nbextension/**js*", "labextension/**"]}

data_files_spec = [
    ("share/jupyter/nbextensions/mols2grid", "mols2grid/nbextension", "**"),
    ("share/jupyter/labextensions/mols2grid", "mols2grid/labextension", "**"),
    ("share/jupyter/labextensions/mols2grid", ".", "install.json"),
    ("etc/jupyter/nbconfig/notebook.d", ".", "mols2grid.json"),
]

cmdclass = create_cmdclass(
    "jsdeps", package_data_spec=package_data_spec, data_files_spec=data_files_spec
)
npm_install = combine_commands(
    install_npm(HERE, build_cmd="build:prod"),
    ensure_targets(jstargets),
)
cmdclass["jsdeps"] = skip_if_exists(jstargets, npm_install)

setup_args = dict(
    version=version,
    cmdclass=cmdclass,
)

if __name__ == "__main__":
    setup(**setup_args)
