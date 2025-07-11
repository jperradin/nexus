"""
REVE - Realistic Environment for Vitreous Exploration
(also verre en verlan)

A package for working with large-scale molecular dynamics trajectories
with a focus on memory efficiency and performance.
"""

# Import core components
from .core.node import Node
from .core.frame import Frame
from .core.system import System

# Import IO components
from .io.reader.reader_factory import ReaderFactory
from .io.parser.parser import Parser
# from .io.writer.writer_factory import WriterFactory # TODO: implement writers

# Import settings
from .config.settings import Settings, SettingsBuilder, AnalysisSettings

# Import main function
from .main import main

# Import version
from .version import __version__
