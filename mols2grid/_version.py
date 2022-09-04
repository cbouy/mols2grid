version_info = (1, 0, 0, "dev")

__version__ = ".".join(map(str, version_info[:3]))
if len(version_info) > 3:
    __version__ += f"-{version_info[3]}"
