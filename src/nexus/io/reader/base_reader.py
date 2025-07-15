from abc import ABC, abstractmethod
from typing import Generator, List, Optional, TextIO

from ...config.settings import Settings
from ...core.frame import Frame

class BaseReader(ABC):
    def __init__(self, settings: Settings) -> None:
        self.verbose: bool = settings.verbose
        self.filename: str = settings.file_location
        self._settings: Settings = settings
        self.num_frames: int = 0
        self.frame_offsets: List[int] = [] # byte offset of each frame
        self.frame_sizes: List[int] = [] # byte size of each frame
        self.mmaped_file: Optional[memoryview] = None # Mapped file
        self.is_indexed: bool = False # True if the file is indexed 

    def set_verbose(self, verbose: bool) -> None:
        self.verbose = verbose

    def seek_to_line(self, file_handle: TextIO, offset: int) -> None:
        """
        Seeks to the specified line in an already opened file.

        Args:
            file_handle: The open file handle
        offset (int): The line offset to seek to.
        """
        file_handle.seek(0)  # Reset to beginning of file
        for _ in range(offset):
            file_handle.readline()
        return

    @abstractmethod
    def detect(self, filepath: str) -> bool:
        """
        Detects if the file is supported by this reader.

        Returns:
            bool: True if the file is supported, False otherwise.
        """
        pass

    @abstractmethod
    def scan(self) -> List[Frame]:
        """
        Scans the trajectory file.
        Initializes Frame objects with the chunk locations of each frame.
        Parse the header to store the number of nodes, the lattice and other informations.

        Returns:
            List[Frame]: A list of Frame objects.
        """
        pass

    @abstractmethod
    def parse(self) -> Generator[Frame, None, None]:
        """
        Parses the trajectory file and yields frames.

        Yields:
            Frame: A data structure representing a frame.
        """
        pass