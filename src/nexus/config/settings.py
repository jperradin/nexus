import os
from dataclasses import dataclass, field
import numpy as np
from typing import Tuple, Optional, List


@dataclass
class GeneralSettings:
    """
    General configuration parameters for the analysis pipeline.

    Attributes:
        project_name (str): Name of the project, used for output directory naming.
        export_directory (str): Root directory for exported results.
        file_location (str): Path to the trajectory file.
        range_of_frames (Tuple[int, int]): Start and end frame indices to process.
            Use -1 as end to process all remaining frames.
        apply_pbc (bool): Whether to apply periodic boundary conditions.
        verbose (bool): Whether to print settings, progress bars, and other information.
        save_logs (bool): Whether to save log files.
        save_performance (bool): Whether to save performance metrics.
    """

    project_name: str = "Project"  # Name of the project
    export_directory: str = "exports"  # Directory to export results
    file_location: str = ""  # Path to the trajectory file
    range_of_frames: Tuple[int, int] = (
        0,
        -1,
    )  # Range of frames to process (0 to -1 = all frames)
    apply_pbc: bool = False  # Whether to apply periodic boundary conditions
    verbose: bool = (
        False  # Whether to print settings, progress bars and other information
    )
    save_logs: bool = False  # Whether to save logs
    save_performance: bool = False  # Whether to save performance


@dataclass
class Cutoff:
    """
    Distance cutoff between two node types for neighbor searching.

    Attributes:
        type1 (str): Symbol of the first node type.
        type2 (str): Symbol of the second node type.
        distance (float): Maximum distance for this pair to be considered neighbors.
    """

    type1: str
    type2: str
    distance: float

    def get_distance(self) -> float:
        """
        Return the cutoff distance.

        Returns:
            float: The cutoff distance for this type pair.
        """
        return self.distance

    def __str__(self) -> str:
        """Return a formatted string showing the type pair and its cutoff distance."""
        max_len = 5
        diff = max_len - len(self.type1) - 1 - len(self.type2)
        return f"{self.type1}-{self.type2}{' ' * diff} : distance = {self.distance}"

@dataclass
class ClusteringSettings:
    """
    Configuration parameters for the clustering algorithm.

    Controls which clustering strategy is selected, the neighbor search method,
    node filtering, cutoff distances, and optional coordination or shared-neighbor
    analysis modes.

    Attributes:
        criterion (str): Clustering criterion, either ``"distance"`` (2-element
            connectivity) or ``"bond"`` (3-element connectivity through a bridging node).
        neighbor_searcher (str): Spatial search algorithm (currently ``"kd_tree"``).
        node_types (List[str]): Node symbols to include in the analysis.
        node_masses (List[float]): Masses for each node type in reduced units.
        connectivity (List[str]): Connectivity pattern (e.g., ``["Si", "Si"]`` or
            ``["Si", "O", "Si"]``).
        cutoffs (List[Cutoff]): Distance cutoffs for each type pair.
        with_printed_unwrapped_clusters (bool): Whether to write unwrapped cluster
            coordinates to file.
        print_mode (str): Output mode for unwrapped clusters (``"all"``,
            ``"connectivity"``, ``"individual"``, or ``"none"``).
        with_coordination_number (bool): Whether to enable coordination-number-based
            clustering.
        coordination_mode (str): Which node types to count as coordination neighbors
            (``"all_types"``, ``"same_type"``, ``"different_type"``, or a specific
            node symbol).
        coordination_range (List[int]): Minimum and maximum coordination numbers to
            accept as a two-element list.
        with_pairwise (bool): Whether to compute pairwise coordination clusters.
        with_mixing (bool): Whether to compute mixing coordination clusters.
        with_alternating (bool): Whether to compute alternating coordination clusters.
        with_connectivity_name (str): Connectivity name for the default coordination mode.
        with_number_of_shared (bool): Whether to enable shared-neighbor analysis.
        shared_mode (str): Which node types to count for shared neighbors.
        shared_threshold (int): Minimum number of shared neighbors required.
        shared_threshold_mode (str): Threshold comparison mode (``"exact"`` or
            ``"at_least"``).
    """

    criterion: str = "distance"  # "distance" or "bond"
    neighbor_searcher: str = "kd_tree"  # "kd_tree", TODO : "cell_list"
    node_types: List[str] = field(default_factory=lambda: [])  # List of node types
    node_masses: List[float] = field(
        default_factory=lambda: []
    )  # List of node masses in reduced units
    connectivity: List[str] = field(default_factory=lambda: [])  # List of connectivity
    cutoffs: List[Cutoff] = field(
        default_factory=lambda: []
    )  # Cutoffs for distance and bond criterion
    with_printed_unwrapped_clusters: bool = (
        False  # Whether to print the unwrapped clusters
    )
    print_mode: str = "none"  # "all", "connectivity", "individual", "none"

    # Coordination number ie number of nearest neighbors
    # - all_types: all types of nodes are considered A-AB, B-AB
    # - same_type: only nodes of the same type are considered A-A, B-B
    # - different_type: only nodes of the different types are considered A-B, B-A

    # Calls clustering algorithm with coordination number
    with_coordination_number: bool = (
        False  # Whether to calculate the coordination number
    )
    coordination_mode: str = (
        "all_types"  # "all_types", "same_type", "different_type", "<node_type>"
    )
    coordination_range: List[int] = field(
        default_factory=lambda: []
    )  # Minimum and maximum coordination numbers to consider

    # Calls clustering algorithm with alternating clusters (with coordination number)
    # - with_pairwise: calculate pairwise coordination number ie A4-B4, B3-A3
    with_pairwise: bool = False
    # - with_mixing: calculate mixing coordination number ie A4-B5, B2-A3
    with_mixing: bool = False
    # - with_alternating: calculate alternating coordination number ie A4-B6, B2-A3
    with_alternating: bool = False
    # if with_coordination_number is True and not with_pairwise, with_mixing or with_alternating (default mode)
    with_connectivity_name: str = ""  # Name of the connectivity

    # Calls clustering algorithm with shared
    with_number_of_shared: bool = False  # Whether to calculate the number of shared
    shared_mode: str = (
        "all_types"  # "all_types", "same_type", "different_type", "<node_type>"
    )
    shared_threshold: int = 1  # Minimum shared threshold
    shared_threshold_mode: str = "exact"  # "exact", "minimum"

    def get_max_cutoff(self) -> float:
        """
        Return the largest cutoff distance across all type pairs.

        Returns:
            float: Maximum cutoff distance, or 0.0 if no cutoffs are defined.
        """
        max_cutoff = 0.0
        for cutoff in self.cutoffs:
            if cutoff.distance > max_cutoff:
                max_cutoff = cutoff.distance
        return max_cutoff

    def get_cutoff(self, type1: str, type2: str) -> float:
        """
        Look up the cutoff distance for a given type pair.

        The lookup is symmetric: ``get_cutoff("Si", "O")`` and ``get_cutoff("O", "Si")``
        return the same value.

        Args:
            type1 (str): Symbol of the first node type.
            type2 (str): Symbol of the second node type.

        Returns:
            float: The cutoff distance, or None if no matching pair is found.
        """
        for cutoff in self.cutoffs:
            if cutoff.type1 == type1 and cutoff.type2 == type2:
                return cutoff.distance
            elif cutoff.type1 == type2 and cutoff.type2 == type1:
                return cutoff.distance
        return None

    def __str__(self) -> str:
        """Return a formatted summary of active clustering settings."""
        lines = []
        for key, value in self.__dict__.items():
            if value is not None:
                if (
                    not self.with_coordination_number
                    and key == "with_coordination_number"
                ):
                    continue
                elif not self.with_coordination_number and key == "coordination_mode":
                    continue
                elif not self.with_coordination_number and key == "coordination_range":
                    continue
                elif not self.with_alternating and key == "with_alternating":
                    continue
                elif not self.with_number_of_shared and key == "with_number_of_shared":
                    continue
                elif not self.with_number_of_shared and key == "shared_mode":
                    continue
                elif not self.with_number_of_shared and key == "shared_threshold":
                    continue
                elif not self.with_number_of_shared and key == "shared_threshold_mode":
                    continue
                if key == "cutoffs":
                    line1 = f"\t\t|- {key:}:"
                    for cutoff in value:
                        line1 += f"\n\t\t\t{str(cutoff)}"
                    lines.append(line1)
                else:
                    lines.append(f"\t\t|- {key}: {value}")
        output = """
        Clustering Settings:
        -----------------
{}
        """.format("\n".join(lines))
        return output


@dataclass
class AnalysisSettings:
    """
    Configuration parameters controlling which analyzers are enabled.

    Each ``with_*`` flag enables the corresponding analyzer in the pipeline. Setting
    ``with_all`` enables every available analyzer at once.

    Attributes:
        overwrite (bool): Whether to overwrite existing output files. If False, results
            are appended.
        with_all (bool): Whether to enable all analyzers.
        with_average_cluster_size (bool): Whether to compute average cluster size.
        with_largest_cluster_size (bool): Whether to compute largest cluster size.
        with_concentration (bool): Whether to compute concentration.
        with_spanning_cluster_size (bool): Whether to compute spanning cluster size.
        with_gyration_radius (bool): Whether to compute gyration radius.
        with_correlation_length (bool): Whether to compute correlation length.
        with_percolation_probability (bool): Whether to compute percolation probability.
        with_order_parameter (bool): Whether to compute order parameter.
        with_cluster_size_distribution (bool): Whether to compute cluster size
            distribution.
    """

    overwrite: bool = True  # Whether to overwrite the existing file, if False, appends results to the file
    with_all: bool = False  # Whether to calculate all the properties
    with_average_cluster_size: bool = (
        False  # Whether to calculate the average cluster size
    )
    with_largest_cluster_size: bool = (
        False  # Whether to calculate the largest cluster size
    )
    with_concentration: bool = False  # Whether to calculate the concentration
    with_spanning_cluster_size: bool = (
        False  # Whether to calculate the spanning cluster size
    )
    with_gyration_radius: bool = False  # Whether to calculate the gyration radius
    with_correlation_length: bool = False  # Whether to calculate the correlation length
    with_percolation_probability: bool = (
        False  # Whether to calculate the percolation probability
    )
    with_order_parameter: bool = False  # Whether to calculate the order parameter
    with_cluster_size_distribution: bool = (
        False  # Whether to calculate the cluster size distribution
    )

    def get_analyzers(self) -> List[str]:
        """
        Return the list of enabled analyzer class names based on current flags.

        Returns:
            List[str]: Class names of enabled analyzers to be instantiated by the
                analyzer factory.
        """
        analyzers = []
        if self.with_average_cluster_size:
            analyzers.append("AverageClusterSizeAnalyzer")
        if self.with_largest_cluster_size:
            analyzers.append("LargestClusterSizeAnalyzer")
        if self.with_concentration:
            analyzers.append("ConcentrationAnalyzer")
        if self.with_spanning_cluster_size:
            analyzers.append("SpanningClusterSizeAnalyzer")
        if self.with_gyration_radius:
            analyzers.append("GyrationRadiusAnalyzer")
        if self.with_correlation_length:
            analyzers.append("CorrelationLengthAnalyzer")
        if self.with_percolation_probability:
            analyzers.append("PercolationProbabilityAnalyzer")
        if self.with_order_parameter:
            analyzers.append("OrderParameterAnalyzer")
        if self.with_cluster_size_distribution:
            analyzers.append("ClusterSizeDistributionAnalyzer")
        if self.with_all:
            analyzers.append("AverageClusterSizeAnalyzer")
            analyzers.append("ConcentrationAnalyzer")
            analyzers.append("LargestClusterSizeAnalyzer")
            analyzers.append("SpanningClusterSizeAnalyzer")
            analyzers.append("GyrationRadiusAnalyzer")
            analyzers.append("CorrelationLengthAnalyzer")
            analyzers.append("PercolationProbabilityAnalyzer")
            analyzers.append("OrderParameterAnalyzer")
            analyzers.append("ClusterSizeDistributionAnalyzer")
        return analyzers

    def __str__(self) -> str:
        """Return a formatted summary of active analysis settings."""
        lines = []
        for key, value in self.__dict__.items():
            if value is not None:
                if not self.with_all and key == "with_all":
                    continue
                elif (
                    not self.with_average_cluster_size
                    and key == "with_average_cluster_size"
                ):
                    continue
                elif not self.with_concentration and key == "with_concentration":
                    continue
                elif (
                    not self.with_largest_cluster_size
                    and key == "with_largest_cluster_size"
                ):
                    continue
                elif (
                    not self.with_spanning_cluster_size
                    and key == "with_spanning_cluster_size"
                ):
                    continue
                elif not self.with_gyration_radius and key == "with_gyration_radius":
                    continue
                elif (
                    not self.with_correlation_length
                    and key == "with_correlation_length"
                ):
                    continue
                elif (
                    not self.with_percolation_probability
                    and key == "with_percolation_probability"
                ):
                    continue
                elif not self.with_order_parameter and key == "with_order_parameter":
                    continue
                elif (
                    not self.with_cluster_size_distribution
                    and key == "with_cluster_size_distribution"
                ):
                    continue
                lines.append(f"\t\t|- {key}: {value}")
        output = """
        Analysis Settings:
        -----------------
{}
        """.format("\n".join(lines))
        return output


@dataclass
class LatticeSettings:
    """
    Configuration parameters for the simulation cell lattice.

    Controls whether a custom lattice is applied, the lattice matrix values, and
    whether the lattice is read from an external file.

    Attributes:
        apply_custom_lattice (bool): Whether to override the lattice read from the
            trajectory with a user-supplied matrix.
        custom_lattice (np.ndarray): User-supplied 3x3 lattice matrix.
        lattice (np.ndarray): Active 3x3 lattice matrix used during analysis.
        get_lattice_from_file (bool): Whether to read the lattice from a separate file.
        lattice_file_location (str): Path to the external lattice file.
        apply_lattice_to_all_frames (bool): Whether to apply the same lattice matrix
            to every frame in the trajectory.
    """

    apply_custom_lattice: bool = False
    custom_lattice: np.ndarray = field(
        default_factory=lambda: np.array(
            [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
        )
    )
    lattice: np.ndarray = field(
        default_factory=lambda: np.array(
            [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
        )
    )
    get_lattice_from_file: bool = False
    lattice_file_location: str = "./"
    apply_lattice_to_all_frames: bool = True

    def __str__(self) -> str:
        """Return a formatted summary of lattice settings."""
        lines = []
        for key, value in self.__dict__.items():
            if value is not None:
                if not self.apply_custom_lattice and key == "apply_custom_lattice":
                    lines.append(f"\t\t|- {key}: {value}")
                    break
                elif key == "custom_lattice":
                    line1 = f"\t\t|- {key}:"
                    lx = np.array2string(
                        value[0],
                        separator=", ",
                        formatter={"float_kind": lambda x: f"{x}"},
                    )
                    ly = np.array2string(
                        value[1],
                        separator=", ",
                        formatter={"float_kind": lambda x: f"{x}"},
                    )
                    lz = np.array2string(
                        value[2],
                        separator=", ",
                        formatter={"float_kind": lambda x: f"{x}"},
                    )
                    lines.append(
                        f"{line1}\n\t\t\tlx = {lx}\n\t\t\tly = {ly}\n\t\t\tlz = {lz}"
                    )
                else:
                    lines.append(f"\t\t|- {key}: {value}")
        output = """

        Lattice Settings:
        -----------------
{}
        """.format("\n".join(lines))
        return output


@dataclass
class Settings:
    """
    Composite configuration object for the entire analysis pipeline.

    Aggregates general, lattice, clustering, and analysis sub-settings into a single
    object. Constructed via ``SettingsBuilder`` which validates constraints between
    fields.

    Attributes:
        project_name (str): Name of the project, used for output directory naming.
        export_directory (str): Root directory for exported results.
        file_location (str): Path to the trajectory file.
        range_of_frames (Tuple[int, int]): Start and end frame indices to process.
        apply_pbc (bool): Whether to apply periodic boundary conditions.
        verbose (bool): Whether to print progress information.
        save_logs (bool): Whether to save log files.
        save_performance (bool): Whether to save performance metrics.
        general (GeneralSettings): General configuration sub-settings.
        lattice (LatticeSettings): Lattice configuration sub-settings.
        clustering (ClusteringSettings): Clustering configuration sub-settings.
        analysis (AnalysisSettings): Analysis configuration sub-settings.
    """

    project_name: str = "default"
    export_directory: str = "export"
    file_location: str = "./"
    range_of_frames: Tuple[int, int] = (0, -1)
    apply_pbc: bool = True
    verbose: bool = False
    save_logs: bool = False
    save_performance: bool = False
    general: GeneralSettings = field(default_factory=GeneralSettings)
    lattice: LatticeSettings = field(default_factory=LatticeSettings)
    clustering: ClusteringSettings = field(default_factory=ClusteringSettings)
    analysis: AnalysisSettings = field(default_factory=AnalysisSettings)

    @property
    def output_directory(self) -> str:
        """Return the full output path as ``export_directory/project_name``."""
        return os.path.join(self.export_directory, self.project_name)

    def set_range_of_frames(self, start: int, end: Optional[int] = None):
        """
        Set the range of frames to process.

        Args:
            start (int): First frame index (must be non-negative).
            end (Optional[int]): Last frame index, or None/-1 to process all remaining
                frames.

        Raises:
            ValueError: If ``start`` is negative or greater than ``end``.
        """
        if end is None:
            end = -1
        if start < 0:
            raise ValueError("Start frame cannot be negative")
        if end != -1 and start > end:
            raise ValueError("Start frame cannot be greater than end frame")
        self.range_of_frames = (start, end)

    def resolve_frame_range(self, num_frames: int) -> Tuple[int, int]:
        """
        Resolve the frame range into concrete start and end indices.

        Replaces sentinel values (-1) with actual bounds derived from the
        total number of frames. Both start and end are inclusive.

        Args:
            num_frames (int): Total number of frames in the trajectory.

        Returns:
            Tuple[int, int]: Concrete (start, end) indices with no sentinel values.
        """
        start, end = self.range_of_frames
        if start < 0:
            start = 0
        if end < 0:
            end = num_frames - 1
        return (start, end)

    def __str__(self) -> str:
        """Return a formatted summary of all settings."""
        lines = []
        for key, value in self.__dict__.items():
            if value is not None:
                if key == "general":
                    continue
                elif key == "lattice":
                    lines.append(f"\t{str(self.lattice)}")
                elif key == "analysis":
                    lines.append(f"\t{str(self.analysis)}")
                elif key == "clustering":
                    lines.append(f"\t{str(self.clustering)}")
                else:
                    lines.append(f"\t|- {key}: {value}")
        output = """
        General Settings:
        ----------------
{}
        """.format("\n".join(lines))
        return output


class SettingsBuilder:
    """
    Builder for constructing and validating a ``Settings`` object.

    Provides a fluent interface for setting sub-configurations. Each ``with_*`` method
    validates its input and returns ``self`` for chaining. Call ``build()`` to obtain
    the final ``Settings`` instance.

    Attributes:
        _settings (Settings): The settings object being constructed.
    """

    def __init__(self):
        """Initialize the builder with default settings."""
        self._settings = Settings()  # Start with default settings

    def with_lattice(self, lattice: LatticeSettings):
        """
        Set lattice configuration.

        Args:
            lattice (LatticeSettings): Lattice sub-settings to apply.

        Returns:
            SettingsBuilder: Self for method chaining.

        Raises:
            ValueError: If ``lattice`` is not a ``LatticeSettings`` instance.
        """
        if not isinstance(lattice, LatticeSettings):
            raise ValueError(f"Invalid lattice settings: {lattice}")
        self._settings.lattice = lattice
        return self

    def with_general(self, general: GeneralSettings):
        """
        Set general configuration and propagate top-level fields.

        Validates required fields and copies them to the top-level ``Settings`` attributes.

        Args:
            general (GeneralSettings): General sub-settings to apply.

        Returns:
            SettingsBuilder: Self for method chaining.

        Raises:
            ValueError: If any required field is missing or invalid.
        """
        if not isinstance(general, GeneralSettings):
            raise ValueError(f"Invalid general settings: {general}")
        if not general.project_name:
            raise ValueError(f"Invalid project name: {general.project_name}")
        if not general.export_directory:
            raise ValueError(f"Invalid export directory: {general.export_directory}")
        if not general.file_location:
            raise ValueError(f"Invalid file location: {general.file_location}")
        if not general.range_of_frames:
            raise ValueError(f"Invalid range of frames: {general.range_of_frames}")
        if general.apply_pbc is None:
            raise ValueError(f"Invalid apply pbc: {general.apply_pbc}")

        self._settings.project_name = general.project_name
        self._settings.export_directory = general.export_directory
        self._settings.file_location = general.file_location
        self._settings.range_of_frames = general.range_of_frames
        self._settings.apply_pbc = general.apply_pbc
        if general.verbose is not None:
            self._settings.verbose = general.verbose
        if general.save_logs is not None:
            self._settings.save_logs = general.save_logs
        if general.save_performance is not None:
            self._settings.save_performance = general.save_performance
        return self

    def with_analysis(self, analysis: AnalysisSettings):
        """
        Set analysis configuration.

        Args:
            analysis (AnalysisSettings): Analysis sub-settings to apply.

        Returns:
            SettingsBuilder: Self for method chaining.

        Raises:
            ValueError: If ``analysis`` is not an ``AnalysisSettings`` instance.
        """
        if not isinstance(analysis, AnalysisSettings):
            raise ValueError(f"Invalid analysis settings: {analysis}")
        self._settings.analysis = analysis
        return self

    def with_clustering(self, clustering: ClusteringSettings):
        """
        Set clustering configuration with full constraint validation.

        Validates criterion/connectivity compatibility, coordination number constraints,
        pairwise/alternating/mixing prerequisites, and shared-neighbor settings.

        Args:
            clustering (ClusteringSettings): Clustering sub-settings to apply.

        Returns:
            SettingsBuilder: Self for method chaining.

        Raises:
            ValueError: If any constraint is violated (see source for full list).
        """
        if not isinstance(clustering, ClusteringSettings):
            raise ValueError(f"Invalid clustering settings: {clustering}")

        if clustering.criterion not in ["bond", "distance"]:
            raise ValueError(f"Invalid criterion: {clustering.criterion}")

        if clustering.connectivity is None:
            raise ValueError(f"Invalid connectivity: {clustering.connectivity}")

        if clustering.criterion == "bond" and len(clustering.connectivity) != 3:
            raise ValueError(
                f"Invalid connectivity, connectivity must be a list of 3 elements, got {len(clustering.connectivity)}"
            )

        if clustering.criterion == "distance" and len(clustering.connectivity) != 2:
            raise ValueError(
                f"Invalid connectivity, connectivity must be a list of 2 elements, got {len(clustering.connectivity)}"
            )

        if clustering.with_coordination_number:
            modes = ["all_types", "same_type", "different_type"]
            if (
                clustering.coordination_mode not in modes
                and clustering.coordination_mode not in clustering.node_types
            ):
                raise ValueError(
                    f"Invalid coordination mode: {clustering.coordination_mode}"
                )
            if len(clustering.coordination_range) != 2:
                raise ValueError(
                    f"Invalid coordination range: {clustering.coordination_range}"
                )
            if clustering.coordination_range[0] < 1:
                raise ValueError(
                    f"Invalid coordination range: {clustering.coordination_range}"
                )
            if clustering.coordination_range[0] > clustering.coordination_range[1]:
                raise ValueError(
                    f"Invalid coordination range: {clustering.coordination_range}"
                )
            if clustering.coordination_mode is None:
                raise ValueError(
                    f"Invalid coordination mode: {clustering.coordination_mode} with with_coordination_number set to True"
                )

        if clustering.with_pairwise and not clustering.with_coordination_number:
            raise ValueError(f"Activate with_coordination_number before with_pairwise")

        if clustering.with_alternating and not clustering.with_coordination_number:
            raise ValueError(
                f"Activate with_coordination_number before with_alternating"
            )

        if clustering.with_mixing and not clustering.with_coordination_number:
            raise ValueError(f"Activate with_coordination_number before with_mixing")

        if (
            clustering.with_coordination_number
            and not clustering.with_pairwise
            and not clustering.with_alternating
            and not clustering.with_mixing
            and clustering.with_connectivity_name == ""
        ):
            raise ValueError(
                f"Default mode with_coordination_number requires a connectivity name"
            )

        if clustering.with_number_of_shared and not clustering.with_coordination_number:
            raise ValueError(
                f"Activate with_coordination_number before with_number_of_shared"
            )

        if (
            clustering.with_number_of_shared
            and clustering.shared_mode not in modes
            and clustering.shared_mode not in clustering.node_types
        ):
            raise ValueError(f"Invalid shared mode: {clustering.shared_mode}")

        if clustering.with_number_of_shared and clustering.shared_threshold < 1:
            raise ValueError(f"Invalid shared threshold: {clustering.shared_threshold}")

        if clustering.with_number_of_shared and clustering.shared_threshold is None:
            raise ValueError(f"Invalid shared threshold: {clustering.shared_threshold}")

        if (
            clustering.with_number_of_shared
            and clustering.shared_threshold_mode not in ["exact", "at_least"]
        ):
            raise ValueError(
                f"Invalid threshold mode: {clustering.shared_threshold_mode}"
            )

        if clustering.node_types is None:
            raise ValueError(f"Invalid node types: {clustering.node_types}")

        if clustering.node_masses is None:
            raise ValueError(f"Invalid node masses: {clustering.node_masses}")

        if len(clustering.node_types) != len(clustering.node_masses):
            raise ValueError(
                f"Invalid node types and masses: {clustering.node_types} and {clustering.node_masses}"
            )

        if clustering.with_printed_unwrapped_clusters and clustering.print_mode not in [
            "all",
            "connectivity",
            "individual",
            "none",
        ]:
            raise ValueError(f"Invalid print_mode: {clustering.print_mode}")

        self._settings.clustering = clustering
        return self

    def build(self) -> Settings:
        """
        Return the constructed settings object.

        Returns:
            Settings: The fully configured settings instance.
        """
        return self._settings


__all__ = [
    Settings,
    SettingsBuilder,
    AnalysisSettings,
    ClusteringSettings,
    LatticeSettings,
    Cutoff,
]

