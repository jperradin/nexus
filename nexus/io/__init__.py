"""
Module: io
Description: This module provides various I/O functions for reading and writing data related to lattice properties, configurations, and results.
"""

from .read_lattices_properties import read_lattices_properties
from .read_number_of_configurations import count_configurations
from .read_and_create_system import read_and_create_system
from .write_list_of_files import write_list_of_files
from .result import *