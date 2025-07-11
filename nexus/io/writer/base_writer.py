from abc import ABC, abstractmethod
from ...config.settings import Settings

class BaseWriter(ABC):
    def __init__(self, settings: Settings) -> None:
        self.verbose: bool = True
        self._settings: Settings = settings

    @abstractmethod
    def write(self) -> None:
        pass

        