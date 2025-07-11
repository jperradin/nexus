API Reference
=============

This section provides some detailed reference for the Nexus-CAT API.
For a complete list of classes and methods, refer to the source code.

## main.py

### `main` function

Analyzes molecular dynamics trajectories with performance tracking, cluster analysis, and multiple analytics capabilities.

**Parameters:**
- `settings (Settings)`: Settings object containing all configuration parameters.

**Functionality:**
- Performance tracking and metrics collection
- Trajectory reading and system initialization
- Frame-by-frame processing with progress visualization
- Neighbor finding and cluster identification
- Multiple analyzer application to each frame
- Results output and cluster visualization
- Performance metrics reporting

## Modules

### `core` module

#### `Cluster` class

Represents a cluster of nodes within a system.

**Attributes:**
- `nodes (List[Node])`: List of Node objects belonging to the cluster.
- `connectivity (str)`: Connectivity criteria of the cluster.
- `root_id (int)`: Node ID that is the root of the cluster.
- `size (int)`: Size of the cluster (number of nodes).
- `settings (Settings)`: Settings object containing configuration parameters.
- `lattice (np.ndarray)`: Lattice of the system containing the cluster.
- `center_of_mass (list)`: Center of mass of the cluster.
- `symbols (list)`: List of symbols of nodes in the cluster.
- `indices (list)`: List of indices of the nodes.
- `unwrapped_positions (list)`: List of unwrapped positions of the cluster.
- `percolation_probability (str)`: Percolation probability along different dimensions.
- `gyration_radius (float)`: Gyration radius of the cluster.
- `order_parameter (list)`: Order parameter values for the cluster.
- `total_nodes (int)`: Total number of nodes in the system.
- `concentration (float)`: Concentration of the cluster in the system.
- `is_percolating (bool)`: Whether the cluster percolates in any dimension.

**Methods:**
- `__init__(self, connectivity: str, root_id: int, size: int, settings: Settings, lattice: np.ndarray) -> None`: Initializes a Cluster object.
- `add_node(self, node: Node) -> None`: Adds a node to the cluster.
- `set_lattice(self, lattice: np.ndarray) -> None`: Sets the lattice for the cluster.
- `get_nodes(self) -> List[Node]`: Returns the list of Node objects belonging to the cluster.
- `get_connectivity(self) -> str`: Returns the connectivity criteria of the cluster.
- `get_size(self) -> int`: Returns the size of the cluster.
- `set_indices_and_positions(self, positions_dict) -> None`: Sets the array of unique indices and positions of nodes in the cluster.
- `calculate_center_of_mass(self) -> None`: Calculates the center of mass of the cluster.
- `calculate_gyration_radius(self) -> None`: Calculates the gyration radius of the cluster.
- `calculate_percolation_probability(self) -> None`: Calculates the percolation probability of the cluster.
- `calculate_order_parameter(self) -> None`: Calculates the order parameter of the cluster.
- `calculate_concentration(self) -> None`: Calculates the concentration of the cluster.
- `calculate_unwrapped_positions(self) -> None`: Calculates the unwrapped positions of atoms in the cluster.
- `unwrap_position(self, vector)`: Unwraps position considering periodic boundary conditions.
- `__str__(self) -> str`: Returns a string representation of the cluster.
- `__repr__(self) -> str`: Returns a detailed string representation of the cluster.

#### `Cutoff` class

Manages cutoff distances for pairs of elements.

**Attributes:**
- `cutoffs (dict)`: Dictionary containing the cutoffs for each pair of elements.
- `pairs (list)`: List of pairs of elements.
- `values (list)`: List of cutoff values.

**Methods:**
- `__init__(self, cutoffs) -> None`: Initializes a Cutoff object.
- `load_cutoffs(self) -> None`: Loads the cutoff values with their associated pair.
- `get_cutoff(self, element1, element2) -> float`: Returns the cutoff for the pair of elements.
- `get_max_cutoff(self) -> float`: Returns the maximum cutoff in the system.

#### `System` class

Manages the atomic system, trajectory data, and interaction with file readers.

**Attributes:**
- `reader (BaseReader)`: The file reader used to load data.
- `settings (Settings)`: The settings object containing configuration parameters.
- `current_frame (Optional[Frame])`: The currently loaded frame. None if no frame is loaded.
- `_current_frame_index (Optional[int])`: Index of the current frame.
- `_num_frames (Optional[int])`: Cache for number of frames.

**Methods:**
- `__init__(self, reader: BaseReader, settings: Settings) -> None`: Initializes a System object.
- `load_frame(self, frame_index: int) -> bool`: Loads a specific frame from the trajectory file.
- `get_frame(self, frame_index: int) -> Optional[Frame]`: Retrieves a specific frame, loading it if necessary.
- `get_num_frames(self) -> int`: Gets the total number of frames in the trajectory.
- `iter_frames(self) -> Generator[Frame, None, None]`: Iterates through the frames of the trajectory, respecting the range defined in settings.
- `__iter__(self) -> 'System'`: Makes the System object itself iterable.
- `__next__(self) -> Frame`: Returns the next frame during iteration.

#### `Frame` class

Representation of a frame of a trajectory.

**Attributes:**
- `frame_id (int)`: ID of the frame.
- `nodes (List[Node])`: List of nodes in the frame.
- `lattice (np.ndarray)`: Lattice of the frame.
- `clusters (Optional[List[Cluster]])`: List of clusters in the frame.
- `_data (Dict[str, np.ndarray])`: Internal data structure for node data (symbol, position).

**Methods:**
- `__post_init__(self) -> None`: Initializes the object after creation.
- `initialize_nodes(self) -> None`: Initializes the list of nodes in the frame.
- `set_lattice(self, lattice: np.ndarray) -> None`: Sets the lattice of the frame.
- `get_lattice(self) -> Optional[np.ndarray]`: Gets the lattice of the frame.
- `get_unique_elements(self) -> List[str]`: Gets the unique elements in the frame.
- `get_node_by_id(self, node_id: int) -> Optional[Node]`: Gets a node by its ID.
- `get_positions(self) -> np.ndarray`: Gets the positions of all nodes in the frame.
- `get_positions_by_element(self) -> Dict[str, np.ndarray]`: Gets the positions of all nodes in the frame grouped by element.
- `get_wrapped_positions(self) -> np.ndarray`: Gets the wrapped positions of all nodes in the frame.
- `get_wrapped_positions_by_element(self) -> Dict[str, np.ndarray]`: Gets the wrapped positions of all nodes in the frame grouped by element.
- `get_clusters(self) -> List[Cluster]`: Gets the clusters of the frame.
- `get_nodes(self) -> List[Node]`: Gets the nodes of the frame.
- `add_cluster(self, cluster: Cluster) -> None`: Adds a cluster to the frame.
- `set_clusters(self, clusters: List[Cluster]) -> None`: Sets the clusters of the frame.
- `get_concentration(self) -> float`: Gets the concentrations of each cluster connectivity in the frame.
- `__len__(self) -> int`: Gets the number of nodes in the frame.
- `__str__(self) -> str`: Returns a string representation of the frame.
- `__repr__(self) -> str`: Returns a detailed string representation of the frame.

#### `Node` class

Representation of a node within a system.

**Attributes:**
- `symbol (str)`: Symbol of the node.
- `node_id (int)`: ID of the node (unique identifier).
- `position (np.ndarray)`: Position of the node.
- `parent (Optional['Node'])`: Parent of the node (Optional).
- `neighbors (List['Node'])`: List of neighbors of the node.
- `cluster_id (Optional[int])`: Cluster ID the node belongs to.
- `distances (Optional[List[float]])`: Distances of the neighbors of the node.
- `indices (Optional[List[int]])`: Indices of the neighbors of the node.
- `mass (float)`: Mass of the node.
- `coordination (int)`: Coordination number of the node.
- `other (Optional[List[str]])`: Other attributes of the node.

**Methods:**
- `__post_init__(self) -> None`: Initializes the object after creation.
- `wrap_position(position: np.ndarray, lattice: np.ndarray) -> np.ndarray`: Wraps position in a periodic box defined by the lattice.
- `add_neighbor(self, node: 'Node') -> None`: Adds a node as a neighbor.
- `reset_parent(self) -> None`: Resets the parent of the node to itself.
- `set_coordination(self, coordination: int) -> None`: Sets the coordination number of the node.
- `__str__(self) -> str`: Returns a string representation of the node.
- `__repr__(self) -> str`: Returns a detailed string representation of the node.

### `io` module

The IO module handles file reading and writing operations for trajectory data and analysis results.

#### `reader` submodule

##### `BaseReader` class

Abstract base class for all file readers.

**Attributes:**
- `verbose (bool)`: Whether to print verbose output.
- `filename (str)`: Path to the file being read.
- `_settings (Settings)`: Settings object.
- `num_frames (int)`: Number of frames in the file.
- `frame_offsets (List[int])`: Byte offsets for each frame.
- `frame_sizes (List[int])`: Byte sizes of each frame.
- `mmaped_file (Optional[memoryview])`: Memory-mapped file for efficient access.
- `is_indexed (bool)`: Whether the file has been indexed.

**Methods:**
- `__init__(self, settings: Settings) -> None`: Initializes a BaseReader object.
- `set_verbose(self, verbose: bool) -> None`: Sets the verbosity.
- `seek_to_line(self, file_handle: TextIO, offset: int) -> None`: Seeks to a specific line in the file.
- `detect(self, filepath: str) -> bool`: Determines if this reader can process the file (abstract).
- `scan(self) -> List[Frame]`: Scans the file to collect metadata (abstract).
- `parse(self) -> Generator[Frame, None, None]`: Parses the file and yields frames (abstract).

##### `ReaderFactory` class

Factory for creating file readers based on file type.

**Methods:**
- `__init__(self, settings: Settings) -> None`: Initializes the factory and registers readers.
- `register_reader(self, reader: BaseReader)`: Registers a new reader instance.
- `get_reader(self) -> Optional[BaseReader]`: Returns the appropriate reader for a given file.

##### `XYZReader` class

Reader for XYZ format trajectory files.

##### `LAMMPSReader` class

Reader for LAMMPS trajectory (lammpstrj) format files.

#### `writer` submodule

##### `BaseWriter` class

Abstract base class for all file writers.

**Attributes:**
- `verbose (bool)`: Whether to print verbose output.
- `_settings (Settings)`: Settings object.

**Methods:**
- `__init__(self, settings: Settings) -> None`: Initializes a BaseWriter object.
- `write(self) -> None`: Writes data to a file (abstract).

##### `WriterFactory` class

Factory for creating file writers.

**Methods:**
- `__init__(self, settings: Settings)`: Initializes the factory and registers writers.
- `register_writer(self, writer: BaseWriter)`: Registers a new writer class.
- `get_writer(self, name: str, mode: str = "all") -> Optional[BaseWriter]`: Returns the appropriate writer instance.

##### `ClustersWriter` class

Writer for cluster data.

##### `LogsWriter` class

Writer for log files.

##### `PerformanceWriter` class

Writer for performance metrics.

##### `MultipleFilesSummaryWriter` class

Writer for summarizing data from multiple files.

### `settings` module

The settings module provides a configuration system for the Nexus package using dataclasses and the builder pattern.

#### `Settings` class

Main settings class that encapsulates all configuration parameters for the application.

**Attributes:**
- `project_name (str)`: Name of the project.
- `export_directory (str)`: Directory for exporting results.
- `file_location (str)`: Path to the trajectory file.
- `range_of_frames (Tuple[int, int])`: Range of frames to process (start, end).
- `apply_pbc (bool)`: Whether to apply periodic boundary conditions.
- `verbose (bool)`: Whether to print verbose information.
- `save_logs (bool)`: Whether to save log files.
- `save_performance (bool)`: Whether to save performance metrics.
- `general (GeneralSettings)`: General settings object.
- `lattice (LatticeSettings)`: Lattice settings object.
- `clustering (ClusteringSettings)`: Clustering settings object.
- `analysis (AnalysisSettings)`: Analysis settings object.

**Methods:**
- `output_directory(self) -> str`: Returns the full output directory path.
- `set_range_of_frames(self, start: int, end: Optional[int] = None) -> None`: Sets the range of frames to process.
- `__str__(self) -> str`: Returns a string representation of the settings.

#### `SettingsBuilder` class

Builder class for creating and validating Settings objects.

**Methods:**
- `__init__(self)`: Initializes the builder with default settings.
- `with_lattice(self, lattice: LatticeSettings)`: Sets the lattice settings.
- `with_general(self, general: GeneralSettings)`: Sets the general settings.
- `with_analysis(self, analysis: AnalysisSettings)`: Sets the analysis settings.
- `with_clustering(self, clustering: ClusteringSettings)`: Sets the clustering settings with validation.
- `build(self) -> Settings`: Builds and returns the final Settings object.

#### `GeneralSettings` class

Settings for general application configuration.

**Attributes:**
- `project_name (str)`: Name of the project.
- `export_directory (str)`: Directory to export results.
- `file_location (str)`: Path to the trajectory file.
- `range_of_frames (Tuple[int, int])`: Range of frames to process.
- `apply_pbc (bool)`: Whether to apply periodic boundary conditions.
- `verbose (bool)`: Whether to print settings and progress information.
- `save_logs (bool)`: Whether to save logs.
- `save_performance (bool)`: Whether to save performance metrics.

#### `LatticeSettings` class

Settings for lattice configuration.

**Attributes:**
- `apply_custom_lattice (bool)`: Whether to apply a custom lattice.
- `custom_lattice (np.ndarray)`: The custom lattice as a 3x3 matrix.
- `get_lattice_from_file (bool)`: Whether to get the lattice from a file.
- `lattice_file_location (str)`: Location of the lattice file.
- `apply_lattice_to_all_frames (bool)`: Whether to apply the lattice to all frames.

**Methods:**
- `__str__(self) -> str`: Returns a string representation of the lattice settings.

#### `ClusteringSettings` class

Settings for configuring cluster analysis.

**Attributes:**
- `criteria (str)`: Clustering criteria ("distance" or "bond").
- `node_types (List[str])`: List of node types.
- `node_masses (List[float])`: List of node masses in reduced units.
- `connectivity (List[str])`: List of connectivity criteria.
- `cutoffs (List[Cutoff])`: List of cutoff distances for pairs of elements.
- `with_printed_unwrapped_clusters (bool)`: Whether to print unwrapped clusters.
- `print_mode (str)`: Mode for printing clusters ("all", "connectivity", "individual", "none").
- `with_coordination_number (bool)`: Whether to calculate coordination numbers.
- `coordination_mode (str)`: Mode for determining coordination ("all_types", "same_type", "different_type", etc.).
- `coordination_range (List[int])`: Minimum and maximum coordination numbers to consider.
- `with_alternating (bool)`: Whether to calculate alternating clusters.
- `with_number_of_shared (bool)`: Whether to calculate shared nodes.
- `shared_mode (str)`: Mode for shared calculation.
- `shared_threshold (int)`: Minimum shared threshold.

**Methods:**
- `get_cutoff(self, type1: str, type2: str) -> float`: Returns the cutoff for a pair of elements.
- `__str__(self) -> str`: Returns a string representation of the clustering settings.

#### `AnalysisSettings` class

Settings for configuring various analyses.

**Attributes:**
- `overwrite (bool)`: Whether to overwrite existing output files.
- `with_all (bool)`: Whether to calculate all analyses.
- `with_average_cluster_size (bool)`: Whether to calculate average cluster size.
- `with_largest_cluster_size (bool)`: Whether to calculate largest cluster size.
- `with_spanning_cluster_size (bool)`: Whether to calculate spanning cluster size.
- `with_gyration_radius (bool)`: Whether to calculate gyration radius.
- `with_correlation_length (bool)`: Whether to calculate correlation length.
- `with_percolation_probability (bool)`: Whether to calculate percolation probability.
- `with_order_parameter (bool)`: Whether to calculate order parameter.
- `with_cluster_size_distribution (bool)`: Whether to calculate cluster size distribution.

**Methods:**
- `get_analyzers(self) -> List[str]`: Returns the list of enabled analyzers.
- `__str__(self) -> str`: Returns a string representation of the analysis settings.

#### `Cutoff` class

Class representing a cutoff distance between two types of nodes.

**Attributes:**
- `type1 (str)`: First node type.
- `type2 (str)`: Second node type.
- `distance (float)`: Cutoff distance.

**Methods:**
- `__str__(self) -> str`: Returns a string representation of the cutoff.
- `get_distance(self) -> float`: Returns the cutoff distance.

### `utils` module

The utils module provides various utility functions and classes used throughout the application.

#### `aesthetics` submodule

Provides functions for visual presentation and text processing.

##### `print_title` function

Prints the ASCII art title of the application with version information.

**Parameters:**
- `__version__ (str)`: The version of the package.

##### `print_title_to_file` function

Writes the ASCII art title of the application with version information to a file.

**Parameters:**
- `__version__ (str)`: The version of the package.
- `path (str)`: Path to the output file.

##### `generate_color_gradient` function

Generates a color gradient between two colors for visual styling.

**Parameters:**
- `num_iterations (int)`: Number of colors to generate in the gradient.

**Returns:**
- `List[Tuple[int, int, int]]`: A list of RGB color tuples.

##### `remove_duplicate_lines` function

Reads a file, removes duplicate lines, and rewrites the file with unique lines.

**Parameters:**
- `filepath (str)`: Path to the file to process.

#### `geometry` submodule

Provides optimized functions for geometric calculations in periodic systems, using numba for acceleration.

##### `wrap_position` function

Wraps a single position vector into the primary unit cell of a periodic lattice.

**Parameters:**
- `position (np.ndarray)`: The position to wrap.
- `lattice (np.ndarray)`: The lattice vectors.

**Returns:**
- `np.ndarray`: The wrapped position.

##### `wrap_positions` function

Wraps multiple position vectors into the primary unit cell of a periodic lattice.

**Parameters:**
- `positions (np.ndarray)`: The positions to wrap.
- `lattice (np.ndarray)`: The lattice vectors.

**Returns:**
- `np.ndarray`: The wrapped positions.

##### `calculate_direct_distance` function

Calculates the Euclidean distance between two positions.

**Parameters:**
- `position1 (np.ndarray)`: The first position.
- `position2 (np.ndarray)`: The second position.

**Returns:**
- `float`: The distance between the positions.

##### `calculate_pbc_distance` function

Calculates the minimum distance between two positions considering periodic boundary conditions.

**Parameters:**
- `position1 (np.ndarray)`: The first position.
- `position2 (np.ndarray)`: The second position.
- `lattice (np.ndarray)`: The lattice vectors.

**Returns:**
- `float`: The minimum distance between the positions.

##### `calculate_direct_angle` function

Calculates the angle formed by three positions.

**Parameters:**
- `position1 (np.ndarray)`: The first position.
- `position2 (np.ndarray)`: The second position (vertex).
- `position3 (np.ndarray)`: The third position.

**Returns:**
- `float`: The angle in degrees.

##### `calculate_pbc_angle` function

Calculates the angle formed by three positions considering periodic boundary conditions.

**Parameters:**
- `position1 (np.ndarray)`: The first position.
- `position2 (np.ndarray)`: The second position (vertex).
- `position3 (np.ndarray)`: The third position.
- `lattice (np.ndarray)`: The lattice vectors.

**Returns:**
- `float`: The angle in degrees.

##### `cartesian_to_fractional` function

Converts Cartesian coordinates to fractional coordinates in a periodic lattice.

**Parameters:**
- `position (np.ndarray)`: The Cartesian position.
- `lattice (np.ndarray)`: The lattice vectors.

**Returns:**
- `np.ndarray`: The fractional coordinates.

##### `fractional_to_cartesian` function

Converts fractional coordinates to Cartesian coordinates in a periodic lattice.

**Parameters:**
- `position (np.ndarray)`: The fractional position.
- `lattice (np.ndarray)`: The lattice vectors.

**Returns:**
- `np.ndarray`: The Cartesian coordinates.

#### `performance` submodule

Provides tools for tracking and analyzing performance metrics.

##### `Performance` class

A dataclass for recording performance metrics of operations.

**Attributes:**
- `id (str)`: Unique identifier for the performance record.
- `name (str)`: Name of the operation or component being measured.
- `timestamp (datetime)`: When the measurement was taken.
- `execution_time_ms (Optional[float])`: Execution time in milliseconds.
- `memory_usage_mb (Optional[float])`: Memory usage in megabytes.
- `cpu_usage_percent (Optional[float])`: CPU usage as a percentage.
- `metrics (Dict[str, Any])`: Additional custom metrics.
- `history (List[Dict[str, Any]])`: Historical performance data points.

**Methods:**
- `add_metric(self, name: str, value: Any) -> None`: Adds a custom metric.
- `record_history(self) -> None`: Records the current state in the history.
- `get_average_execution_time(self) -> Optional[float]`: Calculates the average execution time from history.
- `__str__(self) -> str`: Returns a string representation of the performance data.
