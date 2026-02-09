from abc import ABC, abstractmethod
from typing import List, Dict

from ...core.frame import Frame
from ...config.settings import Settings


class BaseAnalyzer(ABC):
    """
    Abstract base class defining the interface for all analyzers.

    Concrete analyzers implement ``analyze()`` for per-frame computation,
    ``finalize()`` for ensemble averaging, ``get_result()`` for programmatic
    access, and ``print_to_file()`` for CSV output.

    Attributes:
        frame_processed_count (int): Number of frames processed so far.
        _settings (Settings): Configuration settings for the analyzer.
    """

    def __init__(self, settings: Settings) -> None:
        """
        Initialize the analyzer.

        Args:
            settings (Settings): Configuration settings.
        """
        self.frame_processed_count: int = 0
        self._settings: Settings = settings

    @abstractmethod
    def analyze(self, frame: Frame, connectivities: List[str]) -> None:
        """
        Process a single frame for the given connectivity types.

        Args:
            frame (Frame): The frame to analyze.
            connectivities (List[str]): Connectivity labels to analyze.
        """
        pass

    def update_frame_processed(self) -> None:
        """Increment the processed-frame counter by one."""
        self.frame_processed_count += 1

    @abstractmethod
    def finalize(self) -> Dict:
        """
        Compute ensemble-averaged results over all processed frames.

        Returns:
            Dict: The finalized analysis results.
        """
        pass

    @abstractmethod
    def get_result(self) -> Dict:
        """
        Return the current analysis results.

        Returns:
            Dict: A dictionary containing analysis results.
        """
        pass

    @abstractmethod
    def print_to_file(self) -> None:
        """Write the analysis results to the export directory."""
        pass

    def __str__(self) -> str:
        """Return the class name."""
        return f"{self.__class__.__name__}"

    def __repr__(self) -> str:
        """Return a reproducible string representation."""
        return f"{self.__class__.__name__}()"
