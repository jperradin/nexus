from ...config.settings import Settings
from ...io.writer.base_writer import BaseWriter
from ...utils.aesthetics import print_title_to_file
from ...version import __version__

import os
from typing import TextIO

class LogsWriter(BaseWriter):
    """
    Writes the analysis configuration to a log file.

    Outputs the ASCII title banner followed by the full settings summary to
    ``log.txt`` in the export directory.
    """

    def __init__(self, settings: Settings) -> None:
        """
        Initialize the logs writer.

        Args:
            settings (Settings): Configuration settings.
        """
        super().__init__(settings)
        self._settings: Settings = settings

    def write(self) -> None:
        """Write the title banner and settings to the log file."""
        path = os.path.join(self._settings.export_directory, 'log.txt')
        print_title_to_file(__version__, path)
        with open(path, 'a') as f:
            f.write("\n")
            f.write(str(self._settings))
        