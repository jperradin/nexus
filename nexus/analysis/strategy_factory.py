from typing import Optional, List
from ..core.frame import Frame
from ..config.settings import Settings

from .strategies.base_strategy import BaseClusteringStrategy
from .strategies.distance_strategy import DistanceStrategy
from .strategies.bond_strategy import BondingStrategy
from .strategies.shared_strategy import SharedStrategy
from .strategies.coordination_strategy import CoordinationStrategy

class StrategyFactory:
    def __init__(self, frame: Frame, settings: Settings) -> None:
        self._strategies = {}
        # Register other strategies here
        self.register_strategy(DistanceStrategy(frame, settings))
        self.register_strategy(BondingStrategy(frame, settings))
        self.register_strategy(SharedStrategy(frame, settings))
        self.register_strategy(CoordinationStrategy(frame, settings))
        

    def register_strategy(self, strategy: BaseClusteringStrategy) -> None:
        self._strategies[strategy.__class__.__name__] = strategy

    def get_strategy(self, settings: Settings) -> Optional[BaseClusteringStrategy]:
        # get strategy based on clustering settings
        config = settings.clustering
        
        if (config.with_coordination_number or config.with_alternating) and not config.with_number_of_shared:
            return self._strategies.get("CoordinationStrategy")

        if config.with_number_of_shared:
            return self._strategies.get("SharedStrategy")
        
        if not config.with_coordination_number and not config.with_alternating and config.criteria == "distance":
            return self._strategies.get("DistanceStrategy")
        
        if not config.with_coordination_number and not config.with_alternating and config.criteria == "bond":
            return self._strategies.get("BondingStrategy")
            
        return self._strategies.get("Not found.")