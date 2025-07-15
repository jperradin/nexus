from typing import Optional
from ..config.settings import Settings
from .analyzers.base_analyzer import BaseAnalyzer
from .analyzers.average_cluster_size_analyzer import AverageClusterSizeAnalyzer
from .analyzers.concentration_analyzer import ConcentrationAnalyzer
from .analyzers.largest_cluster_size_analyzer import LargestClusterSizeAnalyzer
from .analyzers.spanning_cluster_size_analyzer import SpanningClusterSizeAnalyzer
from .analyzers.percolation_probability_analyzer import PercolationProbabilityAnalyzer
from .analyzers.order_parameter_analyzer import OrderParameterAnalyzer
from .analyzers.cluster_size_distribution_analyzer import ClusterSizeDistributionAnalyzer
from .analyzers.gyration_radius_analyzer import GyrationRadiusAnalyzer
from .analyzers.correlation_length_analyzer import CorrelationLengthAnalyzer

class AnalyzerFactory:
    def __init__(self, settings: Settings, verbose: bool = True):
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
        self._analyzers[analyzer.__class__.__name__] = analyzer

    def get_analyzer(self, analyzer_name: str) -> Optional[BaseAnalyzer]:
        return self._analyzers.get(analyzer_name)