from typing import Optional
from .base_writer import BaseWriter
from .clusters_writer import ClustersWriter
from .logs_writer import LogsWriter
from .performance_writer import PerformanceWriter
from .multiple_files_summary_writer import MultipleFilesSummaryWriter
from ...config.settings import Settings

# TODO: finish implementation of writers 
#       - add support for trajectory writers
#       - add support for system writers
#       - add support for configuration writers
#       - add support for summary writers
#       - add support for statistics writers
#       - add support for performance writers


class WriterFactory:
    """Factory for creating file writers based on file type."""

    def __init__(self, settings: Settings):
        self._writers = {}
        self._settings: Settings = settings
        self.register_writer(ClustersWriter)
        self.register_writer(LogsWriter)
        self.register_writer(PerformanceWriter)
        self.register_writer(MultipleFilesSummaryWriter)

    
    def register_writer(self, writer: BaseWriter):
        """Registers a new writer instance."""
        self._writers[writer.__class__.__name__] = writer

    def get_writer(self, name: str, mode: str = "all") -> Optional[BaseWriter]:
        """Returns the appropriate writer for a given file."""
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