mols2grid_widget
===============================

A Custom Jupyter Widget Library

Installation
------------

To install use pip:

    $ pip install mols2grid_widget

For a development installation (requires [Node.js](https://nodejs.org) and [Yarn version 1](https://classic.yarnpkg.com/)),

    $ git clone https://github.com/cbouy/mols2grid_widget.git
    $ cd mols2grid_widget
    $ pip install -e .
    $ jupyter nbextension install --py --symlink --overwrite --sys-prefix mols2grid_widget
    $ jupyter nbextension enable --py --sys-prefix mols2grid_widget

When actively developing your extension for JupyterLab, run the command:

    $ jupyter labextension develop --overwrite mols2grid_widget

Then you need to rebuild the JS when you make a code change:

    $ cd js
    $ yarn run build

You then need to refresh the JupyterLab page when your javascript changes.
