from typing import Optional
from .base_reader import BaseReader
from .xyz_reader import XYZReader
from .lammps_reader import LAMMPSReader
from ...config.settings import Settings
import os

class ReaderFactory:
    """
    Factory for creating file readers based on file extension.

    Registers all available reader instances on initialization and selects the
    appropriate one by testing each reader's ``detect()`` method against the
    configured file path.

    Attributes:
        _readers (dict): Mapping of file extensions to reader instances.
        _settings (Settings): Configuration settings containing the file location.
    """

    def __init__(self, settings: Settings) -> None:
        """
        Initialize the factory and register all available readers.

        Args:
            settings (Settings): Configuration settings providing the file location.
        """
        self._readers = {}
        self._settings = settings
        self.register_reader(XYZReader(settings))
        self.register_reader(LAMMPSReader(settings))

    def register_reader(self, reader: BaseReader):
        """
        Register a reader instance by probing supported file extensions.

        Args:
            reader (BaseReader): The reader instance to register.
        """
        # Use a dummy filename with the correct extension to determine support
        for ext in ['.xyz', '.lammpstrj', '.other']: #add your extensions here.
            if reader.detect(f'dummy{ext}'):
                self._readers[ext] = reader
                break


    def get_reader(self) -> Optional[BaseReader]:
        """
        Return the reader that supports the configured file.

        Returns:
            Optional[BaseReader]: The matching reader, or None if no reader supports
                the file format.

        Raises:
            ValueError: If the file does not exist.
        """
        if os.path.exists(self._settings.file_location):
            for extension, reader in self._readers.items():
                if reader.detect(self._settings.file_location):
                    return reader
        else:
            raise ValueError(f"File {self._settings.file_location} does not exist.")
        return None