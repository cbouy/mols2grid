import re
from ast import literal_eval
from functools import wraps
from importlib import import_module, resources
from importlib.util import find_spec
from typing import Any

from jinja2 import Environment, FileSystemLoader

templates = resources.files("mols2grid").joinpath("templates")
env = Environment(
    loader=FileSystemLoader(str(templates)),
    trim_blocks=True,
    lstrip_blocks=True,
)


def requires(module):
    def inner(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            if find_spec(module):
                return func(*args, **kwargs)
            raise ModuleNotFoundError(
                f"The module {module!r} is required to use {func.__name__!r} "
                "but it is not installed!"
            )

        return wrapper

    return inner


def import_object(import_path: str) -> Any:
    module_name, obj_name = import_path.rsplit(".", 1)
    return getattr(import_module(module_name), obj_name)


def tooltip_formatter(s, subset, fmt, style, transform):
    """Function to generate tooltips from a pandas Series

    Parameters
    ----------
    s : pandas.Series
        Row in the internal pandas DataFrame
    subset : list
        Subset of columns that are used for the tooltip
    fmt : str
        Format string for each key-value pair of the tooltip
    style : dict
        CSS styling applied to each item independently
    transform : dict
        Functions applied to each value before rendering
    """
    items = []
    for k, v in s[subset].to_dict().items():
        displayed = transform[k](v) if transform.get(k) else v
        value = (
            f'<span class="copy-me" style="{style[k](v)}">{displayed}</span>'
            if style.get(k)
            else f'<span class="copy-me">{displayed}</span>'
        )
        items.append(fmt.format(key=k, value=value))
    items.append("<div class='arrow'></div>")
    return "".join(items)


def slugify(string):
    """Replaces whitespaces with hyphens"""
    return re.sub(r"\s+", "-", string)


def callback_handler(callback, event):
    """Handler for applying the callback function on change"""
    data = literal_eval(event.new)
    callback(data)


def _get_streamlit_script_run_ctx():
    from streamlit.runtime.scriptrunner import get_script_run_ctx

    return get_script_run_ctx()


def is_running_within_streamlit():
    """
    Function to check whether python code is run within streamlit

    Returns
    -------
    use_streamlit : boolean
        True if code is run within streamlit, else False
    """
    try:
        ctx = _get_streamlit_script_run_ctx()
    except ImportError:
        return False
    else:
        return ctx is not None


def is_running_within_marimo():
    """
    Function to check whether python code is run within marimo

    Returns
    -------
    use_marimo : boolean
        True if code is run within marimo, else False
    """
    import sys

    return "marimo" in sys.modules
