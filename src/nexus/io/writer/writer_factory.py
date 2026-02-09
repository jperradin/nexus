from typing import Optional
from .base_writer import BaseWriter
from .clusters_writer import ClustersWriter
from .logs_writer import LogsWriter
from .performance_writer import PerformanceWriter
from .multiple_files_summary_writer import MultipleFilesSummaryWriter
from ...config.settings import Settings


class WriterFactory:
    """
    Factory for creating output writers by class name.

    Registers all available writer classes on initialization. The ``get_writer()``
    method instantiates the appropriate writer with the current settings.

    Attributes:
        _writers (dict): Mapping of class names to writer classes.
        _settings (Settings): Configuration settings passed to writer constructors.
    """

    def __init__(self, settings: Settings):
        """
        Initialize the factory and register all available writer classes.

        Args:
            settings (Settings): Configuration settings for writer construction.
        """
        self._writers = {}
        self._settings: Settings = settings
        self.register_writer(ClustersWriter)
        self.register_writer(LogsWriter)
        self.register_writer(PerformanceWriter)
        self.register_writer(MultipleFilesSummaryWriter)

    def register_writer(self, writer: BaseWriter):
        """
        Register a writer class in the factory.

        Args:
            writer (BaseWriter): The writer class to register.
        """
        self._writers[writer.__class__.__name__] = writer

    def get_writer(self, name: str, mode: str = "all") -> Optional[BaseWriter]:
        """
        Instantiate and return a writer by class name.

        Args:
            name (str): Class name of the writer to create.
            mode (str): Output mode passed to writers that support it.

        Returns:
            Optional[BaseWriter]: The instantiated writer, or None if the name is
                not recognized.
        """
        if name == "ClustersWriter":
            return ClustersWriter(self._settings)
        elif name == "LogsWriter":
            return LogsWriter(self._settings)
        elif name == "PerformanceWriter":
            return PerformanceWriter(self._settings)
        elif name == "MultipleFilesSummaryWriter":
            return MultipleFilesSummaryWriter(self._settings, mode)
        else:
            return None
