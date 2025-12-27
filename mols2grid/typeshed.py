from typing import Any, TypeAlias

try:
    from importlib.resources.abc import Traversable
except ImportError:
    # python 3.10
    from importlib.abc import Traversable

PathLike: TypeAlias = Traversable
Record: TypeAlias = dict[str, Any]
