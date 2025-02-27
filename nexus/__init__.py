"""
Module: nexus
-------------

This module initializes the main components of the application.

Submodules:
-----------
    - core: Core functionalities of the application.
    - data: Data handling functionalities.
    - io: Input/Output functionalities.
    - settings: Configuration settings for the application.
    - utils: Utility functions for the application.

Functions:
----------
    - main: The main entry point of the application.
"""

from . import core
from . import data
from . import io
from . import settings
from . import utils
from .main import main

__version__ = "1.0.0"

utils.print_title(__version__)
