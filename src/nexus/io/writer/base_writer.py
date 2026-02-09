from abc import ABC, abstractmethod
from ...config.settings import Settings

class BaseWriter(ABC):
    """
    Abstract base class for output file writers.

    Defines the ``write()`` interface that all writers must implement.

    Attributes:
        verbose (bool): Whether to print progress information.
        _settings (Settings): Configuration settings for output paths.
    """

    def __init__(self, settings: Settings) -> None:
        """
        Initialize the writer with settings.

        Args:
            settings (Settings): Configuration settings.
        """
        self.verbose: bool = True
        self._settings: Settings = settings

    @abstractmethod
    def write(self) -> None:
        """Write output data to file(s)."""
        pass

        