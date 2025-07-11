from typing import Optional, List
from .finders.base_finder import BaseFinder
from .finders.general_distance_finder import GeneralDistanceFinder
from .finders.general_bond_finder import GeneralBondFinder
from .finders.coordination_based_finder import CoordinationBasedFinder
from .finders.shared_based_finder import SharedBasedFinder
from ..core.frame import Frame
from ..config.settings import Settings

class FinderFactory:
    def __init__(self, frame: Frame, settings: Settings) -> None:
        self._finders = {}
        # Register other finders here
        self.register_finder(GeneralDistanceFinder(frame, settings))
        self.register_finder(GeneralBondFinder(frame, settings))
        self.register_finder(CoordinationBasedFinder(frame, settings))
        self.register_finder(SharedBasedFinder(frame, settings))
        

    def register_finder(self, finder: BaseFinder) -> None:
        self._finders[finder.__class__.__name__] = finder

    def get_finder(self, settings: Settings) -> Optional[BaseFinder]:
        # get finder based on clustering settings
        config = settings.clustering
        
        if (config.with_coordination_number or config.with_alternating) and not config.with_number_of_shared:
            return self._finders.get("CoordinationBasedFinder")

        if config.with_number_of_shared:
            return self._finders.get("SharedBasedFinder")
        
        if not config.with_coordination_number and not config.with_alternating and config.criteria == "distance":
            return self._finders.get("GeneralDistanceFinder")
        
        if not config.with_coordination_number and not config.with_alternating and config.criteria == "bond":
            return self._finders.get("GeneralBondFinder")
            
        return self._finders.get("Not found.")