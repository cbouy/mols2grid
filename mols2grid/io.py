from abc import ABC, abstractmethod
from collections.abc import Callable, Iterator
from csv import Dialect, DictReader, Sniffer
from csv import Error as CSVError
from functools import partial
from itertools import accumulate
from pathlib import Path
from typing import (
    IO,
    Any,
    BinaryIO,
    ClassVar,
    Generic,
    NamedTuple,
    Optional,
    TextIO,
    TypeVar,
)

import pandas as pd
from rdkit import Chem
from upath import UPath

from mols2grid.chem import mol_to_record
from mols2grid.typeshed import PathLike, Record
from mols2grid.utils import import_object

T = TypeVar("T", TextIO, BinaryIO)
CSV_SAMPLE_SIZE = 4096
"""Sample size for CSV/SMILES dialect and header sniffing."""


class ReaderSpecs(NamedTuple):
    extensions: list[str]
    reader: "BaseReader"
    binary: bool


class CompressedSpecs(NamedTuple):
    extensions: list[str]
    import_path: str
    binary_mode: str
    text_mode: str


class ReaderFactory:
    SPECS: ClassVar[dict[str, ReaderSpecs]] = {}
    COMPRESSED_SPECS: ClassVar[dict[str, CompressedSpecs]] = {}
    READERS: ClassVar[dict[str, "BaseReader"]] = {}

    @classmethod
    def register(cls, *extensions: str, reader: "BaseReader", binary: bool) -> None:
        """Register a reader.

        Parameters
        ----------
        extensions:
            File extensions to register.
        reader:
            Reader to register.
        binary:
            Whether the reader processes files opened in binary or text mode.
        """
        for ext in extensions:
            cls.SPECS[ext] = ReaderSpecs(
                extensions=list(extensions), reader=reader, binary=binary
            )

    @classmethod
    def register_compressed(
        cls,
        *extensions: str,
        import_path: str,
        binary_mode: str = "rb",
        text_mode: str = "rt",
    ) -> None:
        """Used to register compressed file formats.

        Parameters
        ----------
        extensions:
            File extensions to register.
        import_path:
            Import path for an ``open``-like function.
        binary_mode:
            Mode to open the file in when binary.
        text_mode:
            Mode to open the file in when text.

        Examples
        --------
        >>> ReaderFactory.register_compressed(
        ...     "zstd", import_path="compression.zstd.open"
        ... )
        >>> read_records("foo.smi.zstd")

        """
        for ext in extensions:
            cls.COMPRESSED_SPECS[ext] = CompressedSpecs(
                extensions=list(extensions),
                import_path=import_path,
                binary_mode=binary_mode,
                text_mode=text_mode,
            )

    @classmethod
    def get(cls, extension: str) -> Optional["BaseReader"]:
        """Get a reader for a given file extension.

        Parameters
        ----------
        extension:
            File extension.

        Returns
        -------
        reader:
            Reader for the given file extension, or ``None`` if no reader is registered
            for that extension.
        """
        extension = extension.lstrip(".").lower()
        if extension in cls.READERS:
            return cls.READERS[extension]
        if extension in cls.SPECS:
            return cls.SPECS[extension].reader
        if "." in extension:
            ext, compressed_ext = extension.rsplit(".", 1)
            if compressed_ext in cls.COMPRESSED_SPECS and ext in cls.SPECS:
                compressed_specs = cls.COMPRESSED_SPECS[compressed_ext]
                base_specs = cls.SPECS[ext]
                mode = (
                    compressed_specs.binary_mode
                    if base_specs.binary
                    else compressed_specs.text_mode
                )
                compressed_func = import_object(compressed_specs.import_path)
                compressed_reader = partial(compressed_func, mode=mode)
                combined_reader = cls._compressed_wrapper(
                    compressed_reader, compressed_specs, base_specs.reader
                )
                cls.READERS[extension] = combined_reader
                return combined_reader
        return None

    @classmethod
    def _compressed_wrapper(
        cls,
        compressed_reader: Callable[..., IO],
        compressed_specs: CompressedSpecs,
        forward_to: "BaseReader",
    ) -> "BaseReader":
        """Forwards file-like objects between different readers.

        Parameters
        ----------
        compressed_reader:
            The reader that handles the compressed input file.
        forward_to:
            The reader that handles the file-like object returned by the
            ``compressed_reader``.
        """
        BaseCls = BinaryReader if isinstance(forward_to, BinaryReader) else TextReader

        class ForwardedReader(BaseCls):
            def __call__(self, fh: IO, **kwargs: Any) -> Iterator[Record]:
                with compressed_reader(fh) as fi:
                    yield from forward_to(fi, **kwargs)

        cname = compressed_specs.import_path.split(".", 1)[0].upper()
        name = f"{cname}Wraps{type(forward_to).__name__}"
        ForwardedReader.__name__ = ForwardedReader.__qualname__ = name

        return ForwardedReader()


ReaderFactory.register_compressed("gz", "gzip", import_path="gzip.open")
ReaderFactory.register_compressed("bz2", "bz", "bzip", "bzip2", import_path="bz2.open")
ReaderFactory.register_compressed("xz", "lz", "lzip", import_path="lzma.open")


class BaseReader(ABC, Generic[T]):
    """Base reader class.

    Notes
    -----
    Readers must add a ``mol`` column to the records they yield, and should yield an
    empty dictionary in case of an invalid entry.
    """

    def __init_subclass__(cls, extensions: list[str] | None = None) -> None:
        if extensions:
            binary = issubclass(cls, BinaryReader)
            ReaderFactory.register(*extensions, reader=cls(), binary=binary)

    @abstractmethod
    def __call__(self, fh: T, /, *args: Any, **kwargs: Any) -> Iterator[Record]: ...


class TextReader(BaseReader[TextIO], ABC):
    pass


class BinaryReader(BaseReader[BinaryIO], ABC):
    pass


class SDFReader(BinaryReader, extensions=["sdf", "mol", "ctab"]):
    """Reader for SDFiles.

    Parameters
    ----------
    fh: BinaryIO
        File-like object
    mol_col: str
        Molecule column in the output dict.
    kwargs: Any
        Parameters passed to :class:`~Chem.rdmolfiles.ForwardSDMolSupplier`,
        e.g. ``sanitize``
    """

    def __call__(
        self, fh: BinaryIO, /, *, mol_col: str = "mol", **kwargs: Any
    ) -> Iterator[Record]:
        for mol in Chem.ForwardSDMolSupplier(fh, **kwargs):
            yield mol_to_record(mol, mol_col=mol_col)


class CSVReader(TextReader, extensions=["csv"]):
    """Reader for CSV files containing SMILES.

    Parameters
    ----------
    fh: TextIO
        File-like object
    smiles_col: str
        SMILES column.
    mol_col: str
        Molecule column in the output dict.
    kwargs: Any
        Parameters passed to :class:`csv.DictReader`, e.g. ``delimiter``.
    """

    def __call__(
        self,
        fh: TextIO,
        /,
        *,
        smiles_col: str = "SMILES",
        mol_col: str = "mol",
        **kwargs: Any,
    ) -> Iterator[Record]:
        if kwargs.get("fieldnames"):
            kwargs["fieldnames"].append("mol")
        csv = DictReader(fh, **kwargs)
        for row in csv:
            smi = row.get(smiles_col, None)
            if smi is None:
                yield {}
                continue
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                yield {}
                continue
            # note: all row values are parsed as strings
            yield {**row, mol_col: mol}


class SmilesReader(TextReader, extensions=["smi", "txt", "tsv"]):
    """Reader for SMILES files.

    Parameters
    ----------
    fh:
        File-like object
    mol_col:
        Molecule column in the output dict.
    kwargs:
        Parameters passed to :class:`csv.DictReader`, e.g. `delimiter`.
    """

    def __call__(
        self, fh: TextIO, /, *, mol_col: str = "mol", **kwargs: Any
    ) -> Iterator[Record]:
        if "dialect" not in kwargs or "fieldnames" not in kwargs:
            sample = fh.read(CSV_SAMPLE_SIZE)
            fh.seek(0)
            # Detect the "dialect" (format) of the smiles file from a small sample
            # before passing it to the CSV reader
            if "dialect" not in kwargs:
                dialect = self._determine_smi_dialect(sample)
                if dialect:
                    kwargs["dialect"] = dialect
                    kwargs.setdefault("delimiter", dialect.delimiter)
                elif "delimiter" not in kwargs:
                    raise ValueError(
                        "No valid columns found, try specifying the `delimiter`"
                    )
            if "fieldnames" not in kwargs:
                kwargs["fieldnames"] = self._determine_smi_columns(
                    sample, kwargs.get("delimiter", " ")
                )
        reader = CSVReader()
        yield from reader(fh, mol_col=mol_col, **kwargs)

    @staticmethod
    def _determine_smi_dialect(sample: str) -> type[Dialect] | None:
        """Automatically determine the dialect from a sample"""
        try:
            return Sniffer().sniff(sample, delimiters=" \t")
        except CSVError:
            return None

    @staticmethod
    def _determine_smi_columns(sample: str, delimiter: str) -> list[str] | None:
        """Automaticall detect the column names from a sample. If there is a header in
        the file, returns ``None`` to let the :class:`csv.DictReader` handle it.

        Notes
        -----
        Assumes the first column in always the SMILES string.
        """
        if Sniffer().has_header(sample):
            return None
        first = sample.split("\n", 1)[0]
        n_fields = len(first.split(delimiter))
        match n_fields:
            case 1:
                return ["SMILES"]
            case _:
                return [
                    "SMILES",
                    "TITLE",
                    *[f"field_{i}" for i in range(n_fields - 2)],
                ]


def get_file_extensions(path: str | PathLike) -> list[str]:
    """Returns a list of possible file extensions from a file name.

    Examples
    --------
    >>> get_file_extensions("foo/bar.CSV")
    ['csv']
    >>> get_file_extensions("foo/bar.sdf.gz")
    ['sdf.gz', 'gz']
    >>> get_file_extensions("foo/bar.test.smi.gz")
    ['test.smi.gz', 'smi.gz', 'gz']
    """
    extensions = Path(str(path)).name.lower().split(".")[1:]
    return list(accumulate(extensions[::-1], lambda x, y: f"{y}.{x}"))[::-1]


def get_reader(path: str | PathLike, /, fmt: str | None = None) -> BaseReader:
    """
    Get the reader for a file based on its extension, or the provided ``fmt``
    format parameter, e.g. ``sdf``
    """
    file_formats = [fmt.lstrip(".").lower()] if fmt else get_file_extensions(path)
    reader = next(
        (reader for ext in file_formats if (reader := ReaderFactory.get(ext))), None
    )
    if reader is None:
        raise ValueError(f"No reader found for file {path!r}") from None
    return reader


def as_path(path: str | PathLike, /) -> PathLike:
    """
    Converts strings to :class:`upath.UPath`, and keeps other Path-like objects as is.
    """
    return UPath(path) if isinstance(path, str) else path


def read_records(path: str | PathLike, /, **kwargs: Any) -> Iterator[Record]:
    """Lazily reads molecules from a path as records (dict).

    Parameters
    ----------
    path : str, Path
        Path to an SDF, CSV, or SMI file (and their compressed equivalent).
    kwargs : Any
        Parameters passed to the underlying reader, e.g. ``mol_col``, ``delimiter``...

    Yields
    ------
    record: Record
        A dictionary containing the molecule object and data.
    """
    reader = get_reader(path)
    mode = "rb" if isinstance(reader, BinaryReader) else "r"
    with as_path(path).open(mode) as fh:
        yield from reader(fh, **kwargs)


def read_mols_to_df(path: str | PathLike, /, **kwargs: Any) -> pd.DataFrame:
    """Creates a dataframe of molecules from a path. All property fields in
    the file are made available in the resulting dataframe

    Parameters
    ----------
    path : str, Path
        Path to an SDF, CSV, or SMI file (and their compressed equivalent).
    kwargs : Any
        Parameters passed to the underlying reader, e.g. ``mol_col``, ``delimiter``...

    Returns
    -------
    df : pandas.DataFrame
    """
    records = read_records(path, **kwargs)
    return pd.DataFrame(records)
