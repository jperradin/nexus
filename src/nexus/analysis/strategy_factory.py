from typing import Optional
from ..core.frame import Frame
from ..config.settings import Settings

from .strategies.base_strategy import BaseClusteringStrategy
from .strategies.distance_strategy import DistanceStrategy
from .strategies.bond_strategy import BondingStrategy
from .strategies.shared_strategy import SharedStrategy
from .strategies.coordination_strategy import CoordinationStrategy


class StrategyFactory:
    """
    Factory class responsible for creating, registering, and retrieving clustering strategy instances.

    This factory implements the factory design pattern to manage a collection of clustering
    strategies. It handles the instantiation and registration of all available strategies,
    and selects the appropriate one based on the clustering settings flags.

    The selection priority is:
        1. CoordinationStrategy — when coordination number constraints are enabled
           (with_coordination_number or with_alternating) without shared neighbor analysis.
        2. SharedStrategy — when shared neighbor analysis is enabled (with_number_of_shared).
        3. DistanceStrategy — base criterion "distance" with a 2-element connectivity (e.g., [Si, Si]).
        4. BondingStrategy — base criterion "bond" with a 3-element connectivity (e.g., [Si, O, Si]).

    Attributes:
        _strategies (dict): Internal dictionary mapping strategy class names to strategy instances.
    """

    def __init__(self, frame: Frame, settings: Settings) -> None:
        """
        Initialize the StrategyFactory with a frame and settings, and register all available strategies.

        Creates and registers all predefined clustering strategy instances. Each strategy
        receives the current frame and settings so it can operate on the frame's nodes
        when selected.

        Args:
            frame (Frame): The simulation frame containing the nodes to be clustered.
            settings (Settings): Configuration settings object containing clustering parameters
                needed by the strategies for their initialization and operation.
        """
        self._strategies = {}
        # Register other strategies here
        self.register_strategy(DistanceStrategy(frame, settings))
        self.register_strategy(BondingStrategy(frame, settings))
        self.register_strategy(SharedStrategy(frame, settings))
        self.register_strategy(CoordinationStrategy(frame, settings))

    def register_strategy(self, strategy: BaseClusteringStrategy) -> None:
        """
        Register a strategy instance in the factory.

        Adds the provided strategy to the internal registry using its class name as the key.
        This allows for later retrieval of the strategy by its class name.

        Args:
            strategy (BaseClusteringStrategy): The strategy instance to register. Must be an
                instance of BaseClusteringStrategy or its subclasses.
        """
        self._strategies[strategy.__class__.__name__] = strategy

    def get_strategy(self, settings: Settings) -> Optional[BaseClusteringStrategy]:
        """
        Select and retrieve the appropriate clustering strategy based on the clustering settings.

        Evaluates the clustering configuration flags in a fixed priority order to determine
        which strategy matches the user's intent. The first matching condition wins.

        Args:
            settings (Settings): Configuration settings whose clustering flags drive the
                strategy selection.

        Returns:
            Optional[BaseClusteringStrategy]: The selected strategy instance if a match is
                found, None otherwise.
        """
        config = settings.clustering

        # Coordination-based clustering (includes pairwise, mixing, and alternating modes)
        # takes precedence unless shared neighbor analysis is also requested
        if (
            config.with_coordination_number or config.with_alternating
        ) and not config.with_number_of_shared:
            return self._strategies.get("CoordinationStrategy")

        # Shared neighbor analysis builds on coordination and adds a shared-count threshold
        if config.with_number_of_shared:
            return self._strategies.get("SharedStrategy")

        # Simple distance-based clustering: connects nodes within a cutoff distance
        if (
            not config.with_coordination_number
            and not config.with_alternating
            and config.criterion == "distance"
        ):
            return self._strategies.get("DistanceStrategy")

        # Bond-based clustering: connects nodes through a bridging atom (e.g., Si-O-Si)
        if (
            not config.with_coordination_number
            and not config.with_alternating
            and config.criterion == "bond"
        ):
            return self._strategies.get("BondingStrategy")

        # No matching strategy found for the given configuration
        return self._strategies.get("Not found.")

