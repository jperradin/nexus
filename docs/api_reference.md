API Reference
=============

This section provides a detailed reference for the Nexus-CAT API. For a complete list of classes and methods, refer to the source code.

## main.py

### `main` function

Analyzes atomistic simulations trajectories with performance tracking, cluster analysis, and multiple analytics capabilities.

**Parameters:**
- `settings (Settings)`: Settings object containing all configuration parameters.

**Functionality:**
- Performance tracking and metrics collection.
- Trajectory reading and system initialization.
- Frame-by-frame processing with progress visualization.
- Neighbor finding and cluster identification using various strategies.
- Multiple analyzer application to each frame.
- Results output and cluster visualization.
- Performance metrics reporting.

## Modules

### `core` module

#### `Cluster` class

Represents a cluster of nodes within a system.

**Attributes:**
- `nodes (List[Node])`: List of Node objects belonging to the cluster.
- `connectivity (str)`: Connectivity criterion of the cluster.
- `root_id (int)`: Node ID that is the root of the cluster.
- `size (int)`: Size of the cluster (number of nodes).
- `settings (Settings)`: Settings object containing configuration parameters.
- `lattice (np.ndarray)`: Lattice of the system containing the cluster.
- `center_of_mass (list)`: Center of mass of the cluster.
- `unwrapped_positions (list)`: List of unwrapped positions of the cluster.
- `percolation_probability (str)`: Percolation probability along different dimensions.
- `gyration_radius (float)`: Gyration radius of the cluster.
- `order_parameter (list)`: Order parameter values for the cluster.
- `is_percolating (bool)`: Whether the cluster percolates in any dimension.
- `decoration_nodes (Dict)`: Dictionary of nodes that "decorate" the main cluster, such as bridging nodes.

**Methods:**
- `add_node(self, node: Node) -> None`: Adds a node to the cluster.
- `calculate_center_of_mass(self) -> None`: Calculates the center of mass of the cluster.
- `calculate_gyration_radius(self) -> None`: Calculates the gyration radius of the cluster.
- `calculate_percolation_probability(self) -> None`: Calculates the percolation probability of the cluster.
- `calculate_order_parameter(self) -> None`: Calculates the order parameter of the cluster.
- `calculate_unwrapped_positions(self) -> None`: Calculates the unwrapped positions of nodes in the cluster, handling periodic boundary conditions.

---
#### `System` class

Manages the atomic system, trajectory data, and interaction with file readers.

**Attributes:**
- `reader (BaseReader)`: The file reader used to load data.
- `settings (Settings)`: The settings object containing configuration parameters.
- `current_frame (Optional[Frame])`: The currently loaded frame.

**Methods:**
- `load_frame(self, frame_index: int) -> bool`: Loads a specific frame from the trajectory file.
- `get_frame(self, frame_index: int) -> Optional[Frame]`: Retrieves a specific frame, loading it if necessary.
- `get_num_frames(self) -> int`: Gets the total number of frames in the trajectory.
- `iter_frames(self) -> Generator[Frame, None, None]`: Iterates through the frames of the trajectory, respecting the range defined in settings.

---
#### `Frame` class

Representation of a frame of a trajectory.

**Attributes:**
- `frame_id (int)`: ID of the frame.
- `nodes (List[Node])`: List of nodes in the frame.
- `lattice (np.ndarray)`: Lattice of the frame.
- `clusters (Optional[List[Cluster]])`: List of clusters in the frame.
- `_data (Dict[str, np.ndarray])`: Internal data structure for node data (symbol, position).

**Methods:**
- `initialize_nodes(self) -> None`: Initializes the list of nodes in the frame.
- `set_clusters(self, clusters: List[Cluster]) -> None`: Sets the clusters of the frame.
- `get_concentration(self) -> float`: Gets the concentrations of each cluster connectivity in the frame.

---
#### `Node` class

Representation of a node within a system.

**Attributes:**
- `symbol (str)`: Symbol of the node.
- `node_id (int)`: ID of the node (unique identifier).
- `position (np.ndarray)`: Position of the node.
- `parent (Optional['Node'])`: Parent of the node (used in union-find algorithm).
- `neighbors (List['Node'])`: List of neighbors of the node.
- `coordination (int)`: Coordination number of the node.

**Methods:**
- `add_neighbor(self, node: 'Node') -> None`: Adds a node as a neighbor.
- `reset_parent(self) -> None`: Resets the parent of the node to itself.
- `set_coordination(self, coordination: int) -> None`: Sets the coordination number of the node.

### `io` module

The IO module handles file reading and writing operations for trajectory data and analysis results.

#### `reader` submodule

##### `ReaderFactory` class

Factory for creating file readers based on file type.

**Methods:**
- `get_reader(self) -> Optional[BaseReader]`: Returns the appropriate reader for a given file.

##### `XYZReader` and `LAMMPSReader` classes

Readers for XYZ and LAMMPS trajectory file formats, respectively. They both implement `detect`, `scan`, and `parse` methods for efficient file processing.

#### `writer` submodule

##### `WriterFactory` class

Factory for creating file writers.

**Methods:**
- `get_writer(self, name: str, mode: str = "all") -> Optional[BaseWriter]`: Returns the appropriate writer instance.

##### `ClustersWriter` class

Writer for cluster data, including unwrapped atomic positions and optional bond information.

##### `LogsWriter` and `PerformanceWriter` classes

Writers for log files and performance metrics.

### `settings` module

The settings module provides a configuration system for the Nexus package using dataclasses and the builder pattern.

#### `Settings` and `SettingsBuilder` classes

The `Settings` class encapsulates all configuration parameters, while the `SettingsBuilder` provides a fluent API for constructing and validating a `Settings` object.

#### `GeneralSettings`, `LatticeSettings`, `ClusteringSettings`, `AnalysisSettings`

These are dataclasses that hold the configuration for different aspects of the analysis, such as general project information, lattice parameters, clustering criterion, and which analyses to perform.

### `utils` module

The utils module provides various utility functions and classes used throughout the application.

#### `geometry` submodule

Provides optimized functions for geometric calculations in periodic systems, such as wrapping positions and calculating distances and angles, using `numba` for acceleration.