# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.

# ruff: noqa: PTH100
import os
import sys
from datetime import datetime

from recommonmark.transform import AutoStructify

sys.path.insert(0, os.path.abspath("."))

# -- Project information -----------------------------------------------------

project = "mols2grid"
copyright = f"2021-{datetime.now().year}, Cédric Bouysset"  # noqa: A001
author = "Cédric Bouysset"


# -- General configuration ---------------------------------------------------

github_doc_root = "https://github.com/cbouy/mols2grid/tree/master/docs/"
needs_sphinx = "5.3.0"

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.intersphinx",
    "sphinx.ext.viewcode",
    "sphinx.ext.napoleon",
    "sphinx.ext.autosectionlabel",
    "recommonmark",
    "IPython.sphinxext.ipython_console_highlighting",
    "IPython.sphinxext.ipython_directive",
    "myst_nb",
    "sphinx_copybutton",
]

myst_enable_extensions = [
    "colon_fence",
]
nb_execution_allow_errors = False
nb_execution_raise_on_error = True
copybutton_exclude = ".linenos, .gp"

autosectionlabel_prefix_document = True
napoleon_google_docstring = False

source_suffix = {
    ".rst": "restructuredtext",
    ".md": "myst-nb",
    ".ipynb": "myst-nb",
    ".myst": "myst-nb",
}

# Add any paths that contain templates here, relative to this directory.
templates_path = []

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = "sphinx_book_theme"
pygments_style = "sphinx"
html_logo = "_static/mols2grid_logo.png"
html_theme_options = {
    "repository_url": "https://github.com/cbouy/mols2grid",
    "path_to_docs": "docs",
    "use_source_button": True,
    "use_download_button": True,
    "use_repository_button": True,
    "use_issues_button": True,
    "launch_buttons": {"colab_url": "https://colab.research.google.com"},
    "icon_links": [
        {
            "name": "GitHub",
            "url": "https://github.com/cbouy/mols2grid",
            "icon": "fa-brands fa-square-github",
            "type": "fontawesome",
        },
    ],
}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ["_static"]

###

intersphinx_mapping = {
    "python": ("https://docs.python.org/3/", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "rdkit": ("https://www.rdkit.org/docs/", None),
    "pandas": ("https://pandas.pydata.org/docs/", None),
    "ipython": ("https://ipython.readthedocs.io/en/stable/", None),
}


# app setup hook
def setup(app):
    app.add_config_value(
        "recommonmark_config",
        {
            # 'url_resolver': lambda url: github_doc_root + url,
            "auto_toc_tree_section": "Contents",
            "enable_math": False,
            "enable_inline_math": False,
            "enable_eval_rst": True,
        },
        True,
    )
    app.add_transform(AutoStructify)
    app.add_css_file("custom.css")
