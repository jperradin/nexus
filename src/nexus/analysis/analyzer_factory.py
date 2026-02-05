from typing import Optional
from ..config.settings import Settings
from .analyzers.base_analyzer import BaseAnalyzer
from .analyzers.average_cluster_size_analyzer import AverageClusterSizeAnalyzer
from .analyzers.concentration_analyzer import ConcentrationAnalyzer
from .analyzers.largest_cluster_size_analyzer import LargestClusterSizeAnalyzer
from .analyzers.spanning_cluster_size_analyzer import SpanningClusterSizeAnalyzer
from .analyzers.percolation_probability_analyzer import PercolationProbabilityAnalyzer
from .analyzers.order_parameter_analyzer import OrderParameterAnalyzer
from .analyzers.cluster_size_distribution_analyzer import (
    ClusterSizeDistributionAnalyzer,
)
from .analyzers.gyration_radius_analyzer import GyrationRadiusAnalyzer
from .analyzers.correlation_length_analyzer import CorrelationLengthAnalyzer


class AnalyzerFactory:
    """
    Module for factory pattern implementation of analyzer creation and registration.
    Factory class responsible for creating, registering, and retrieving analyzer instances.
    
    This factory implements the factory design pattern to manage a collection of analyzer
    objects. It handles the instantiation and registration of various analyzer types used
    for data analysis operations.
    
    Attributes:
        _analyzers (dict): Internal dictionary mapping analyzer class names to analyzer instances.
    """
    def __init__(self, settings: Settings, verbose: bool = True):
        """
        Initialize the AnalyzerFactory with settings and register all available analyzers.
        
        Creates and registers a collection of predefined analyzer instances that will be
        used for various analysis tasks on the dataset.
        
        Args:
            settings (Settings): Configuration settings object containing parameters needed
                by the analyzers for their initialization and operation.
            verbose (bool, optional): Flag to control verbosity of analyzer operations.
                Defaults to True.
        """
        self._analyzers = {}
        # Register other analyzers here
        self.register_analyzer(AverageClusterSizeAnalyzer(settings))
        self.register_analyzer(ConcentrationAnalyzer(settings))
        self.register_analyzer(LargestClusterSizeAnalyzer(settings))
        self.register_analyzer(SpanningClusterSizeAnalyzer(settings))
        self.register_analyzer(PercolationProbabilityAnalyzer(settings))
        self.register_analyzer(OrderParameterAnalyzer(settings))
        self.register_analyzer(ClusterSizeDistributionAnalyzer(settings))
        self.register_analyzer(GyrationRadiusAnalyzer(settings))
        self.register_analyzer(CorrelationLengthAnalyzer(settings))

    def register_analyzer(self, analyzer: BaseAnalyzer) -> None:
        """
        Register an analyzer instance in the factory.
        
        Adds the provided analyzer to the internal registry using its class name as the key.
        This allows for later retrieval of the analyzer by its class name.
        
        Args:
            analyzer (BaseAnalyzer): The analyzer instance to register. Must be an instance
                of BaseAnalyzer or its subclasses.
        """
        self._analyzers[analyzer.__class__.__name__] = analyzer

    def get_analyzer(self, analyzer_name: str) -> Optional[BaseAnalyzer]:
        """
        Retrieve a registered analyzer by its name.
        
        Fetches an analyzer from the registry using the provided analyzer name (typically
        the class name of the analyzer).
        
        Args:
            analyzer_name (str): The name of the analyzer to retrieve. Should correspond
                to the class name of a registered analyzer.
        
        Returns:
            Optional[BaseAnalyzer]: The analyzer instance if found, None otherwise.
        """
        return self._analyzers.get(analyzer_name)
