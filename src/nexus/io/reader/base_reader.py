from abc import ABC, abstractmethod
from typing import Generator, List, Optional, TextIO

from ...config.settings import Settings
from ...core.frame import Frame

class BaseReader(ABC):
    """
    Abstract base class for trajectory file readers.

    Defines the interface that all readers must implement: ``detect()``, ``scan()``,
    and ``parse()``. Provides shared state for file indexing and verbose output.

    Attributes:
        verbose (bool): Whether to print progress information.
        filename (str): Path to the trajectory file.
        num_frames (int): Total number of frames found after scanning.
        frame_offsets (List[int]): Byte offset of each frame in the file.
        frame_sizes (List[int]): Byte size of each frame in the file.
        mmaped_file (Optional[memoryview]): Memory-mapped file handle, if used.
        is_indexed (bool): True once the file has been scanned and indexed.
    """

    def __init__(self, settings: Settings) -> None:
        """
        Initialize the reader with settings.

        Args:
            settings (Settings): Configuration settings providing file location and
                verbosity.
        """
        self.verbose: bool = settings.verbose
        self.filename: str = settings.file_location
        self._settings: Settings = settings
        self.num_frames: int = 0
        self.frame_offsets: List[int] = []
        self.frame_sizes: List[int] = []
        self.mmaped_file: Optional[memoryview] = None
        self.is_indexed: bool = False

    def set_verbose(self, verbose: bool) -> None:
        """
        Set the verbosity flag.

        Args:
            verbose (bool): Whether to print progress information.
        """
        self.verbose = verbose

    def seek_to_line(self, file_handle: TextIO, offset: int) -> None:
        """
        Advance an open file handle to a specific line number.

        Args:
            file_handle (TextIO): The open file handle to seek in.
            offset (int): The line number to seek to (0-based).
        """
        file_handle.seek(0)  # Reset to beginning of file
        for _ in range(offset):
            file_handle.readline()
        return

    @abstractmethod
    def detect(self, filepath: str) -> bool:
        """
        Check whether this reader supports the given file.

        Args:
            filepath (str): Path to the file to test.

        Returns:
            bool: True if the file format is supported by this reader.
        """
        pass

    @abstractmethod
    def scan(self) -> List[Frame]:
        """
        Scan the trajectory file and build an index of frame locations.

        Reads the file sequentially to record byte offsets, node counts, and lattice
        matrices for each frame, enabling fast random access during parsing.

        Returns:
            List[Frame]: Indexed frame metadata objects.
        """
        pass

    @abstractmethod
    def parse(self) -> Generator[Frame, None, None]:
        """
        Parse the trajectory file and yield frames one at a time.

        Yields:
            Frame: A frame containing raw node data and lattice information.
        """
        pass