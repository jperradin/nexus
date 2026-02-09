"""
Readers for molecular dynamics trajectory formats.

Provides a factory-based system for detecting file formats and parsing frames from
XYZ and LAMMPS trajectory files.
"""

from .base_reader import BaseReader
from .xyz_reader import XYZReader
from .lammps_reader import LAMMPSReader
from .reader_factory import ReaderFactory

__all__ = [
    BaseReader,
    ReaderFactory,
    XYZReader,
    LAMMPSReader
]