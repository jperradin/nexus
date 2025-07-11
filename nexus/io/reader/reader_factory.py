from typing import Optional
from .base_reader import BaseReader
from .xyz_reader import XYZReader
from .lammps_reader import LAMMPSReader
from ...config.settings import Settings
import os

class ReaderFactory:
    """Factory for creating file readers based on file type."""

    def __init__(self, settings: Settings) -> None:
        self._readers = {}
        self._settings = settings
        self.register_reader(XYZReader(settings))
        self.register_reader(LAMMPSReader(settings))  
        # Register other readers here

    def register_reader(self, reader: BaseReader):
        """Registers a new reader instance."""
        # Use a dummy filename with the correct extension to determine support
        for ext in ['.xyz', '.lammpstrj', '.other']: #add your extensions here.
            if reader.detect(f'dummy{ext}'):
                self._readers[ext] = reader
                break


    def get_reader(self) -> Optional[BaseReader]:
        """Returns the appropriate reader for a given file."""
        if os.path.exists(self._settings.file_location):
            for extension, reader in self._readers.items():
                if reader.detect(self._settings.file_location):
                    return reader
        else:
            raise ValueError(f"File {self._settings.file_location} does not exist.")
        return None