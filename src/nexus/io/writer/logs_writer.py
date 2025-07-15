from ...config.settings import Settings
from ...io.writer.base_writer import BaseWriter
from ...utils.aesthetics import print_title_to_file
from ...version import __version__

import os
from typing import TextIO

class LogsWriter(BaseWriter):
    def __init__(self, settings: Settings) -> None:
        super().__init__(settings)
        self._settings: Settings = settings

    def write(self) -> None:
        path = os.path.join(self._settings.export_directory, 'log.txt')
        print_title_to_file(__version__, path)
        with open(path, 'a') as f:
            f.write("\n")
            f.write(str(self._settings))
        