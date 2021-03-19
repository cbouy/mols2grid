from functools import wraps
from importlib.util import find_spec

def requires(module):
    def inner(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            if find_spec(module):
                return func(*args, **kwargs)
            raise ModuleNotFoundError(
                f"The module {module!r} is required to use {func.__name__!r} "
                "but it is not installed!")
        return wrapper
    return inner

def tooltip_formatter(s, subset, fmt, style):
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
    """
    items = []
    for k, v in s[subset].to_dict().items():
        v = f'<span style="{style[k](v)}">{v}</span>' if style.get(k) else v
        items.append(fmt.format(key=k, value=v))
    return "<br>".join(items)