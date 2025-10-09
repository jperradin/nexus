API Reference
=============

This section provides a detailed reference for the Nexus-CAT API. For a complete list of classes and methods, refer to the source code.

## main.py

### `main` function

Main entry function to test the package by processing frames of data according to the specified settings. It initializes necessary components, executes frame-by-frame analysis with neighbor finding, clustering, and analysis strategies, tracks and records performance metrics, and saves relevant output and logs.

**Parameters:**
- `settings (Settings)`: Configuration object containing parameters for project setup, verbosity, frame ranges, analysis options, export paths, logging, and performance tracking.

**Descriptions**
- Initializes a performance tracker with a unique ID and timestamped run name if enabled.
- Sets up resource tracking including memory and CPU usage if `performance` is enabled.
- Prints version title and detailed settings if `verbosity` is enabled.
- Creates necessary export directories.
- Saves initial logs if `save_logs` is enabled.
- Initializes data reader and system objects; measures reader initialization time.
- Determines the total number of frames to process based on the settings.
- Creates analyzer instances based on the analysis configuration.
- Uses `tqdm` progress bar for visual feedback during frame processing.
- For each frame:
    - Applies custom lattice if configured.
    - Initializes nodes and performs neighbor finding using a strategy factory.
    - Builds clusters from the found neighbors.
    - Runs analysis on clusters.
    - Optionally prints unwrapped clusters.
    - Tracks and records time for each step and resource usage every 10 frames or on completion if `performance` is enabled.
- After processing frames:
    - Calls analyzers to print results to file.
    - Collects overall performance metrics including memory usage changes.
    - Saves overall performance metrics if enabled.

**Behavior and Side Effects**
- Creates and manages directories and files for output.
- Utilizes multiple factory classes (`ReaderFactory`, `StrategyFactory`, `AnalyzerFactory`, `WriterFactory`) to modularize component creation.
- Uses `psutil` to obtain CPU and memory usage stats.
- Tracks performance per frame and aggregates metrics like average times and node counts.
- Saves logs, cluster data, and performance data conditionally based on settings.

**Usage Example**
```python
from nexus import SettingsBuilder, main
import nexus.config.settings as c

# General settings
config_general = c.GeneralSettings(
    ...
)

# Lattice settings
config_lattice = c.LatticeSettings(
    ...
)

# Clustering settings
config_clustering = c.ClusteringSettings(
    ...
)

# Analysis settings
config_analysis = c.AnalysisSettings(
    ...
)

# Build Settings object
settings = (
    SettingsBuilder()
    .with_general(config_general)        # General settings \
    .with_lattice(config_lattice)        # Lattice settings \
    .with_clustering(config_clustering)  # Clustering settings \
    .with_analysis(config_analysis)      # Analysis settings \
    .build()                             # Build the settings object
)

# Run main processing function
main(settings)
```

## `config` Module (Settings)

### GeneralSettings

Represents general project-level settings used globally in the package.

#### Attributes
- `project_name` (str): Name of the project. Default: `"Project"`.
- `export_directory` (str): Directory path to export results. Default: `"exports"`.
- `file_location` (str): Path to the trajectory file. Default: `""`.
- `range_of_frames` (Tuple[int, int]): Frame range to process; `(0, -1)` means all frames. Default: `(0, -1)`.
- `apply_pbc` (bool): Whether to apply periodic boundary conditions. Default: `False`.
- `verbose` (bool): Whether to print detailed info and progress bars. Default: `False`.
- `save_logs` (bool): Whether to save logs during runs. Default: `False`.
- `save_performance` (bool): Whether to save performance data. Default: `False`.


### Cutoff

Represents a cutoff distance setting between two node types.

#### Attributes
- `type1` (str): First node type.
- `type2` (str): Second node type.
- `distance` (float): Cutoff distance between the node types.

#### Methods
- `__str__()`: Returns formatted string representation.
- `get_distance() -> float`: Returns the cutoff distance.


### ClusteringSettings

Holds settings related to clustering operations.

#### Attributes
- `criterion` (str): Clustering criterion, either `"distance"` or `"bond"`. Default: `"distance"`.
- `neighbor_searcher` (str): Method for neighbor searching, e.g., `"kd_tree"`. Default: `"kd_tree"`.
- `node_types` (List[str]): List of node types.
- `node_masses` (List[float]): List of node masses in reduced units.
- `connectivity` (List[str]): Connectivity specifications.
- `cutoffs` (List[Cutoff]): List of cutoff distances.
- `with_printed_unwrapped_clusters` (bool): Flag to print unwrapped clusters. Default: `False`.
- `print_mode` (str): Print mode among `"all"`, `"connectivity"`, `"individual"`, or `"none"`. Default: `"none"`.
- Additional flags and settings for coordination number calculations, pairwise/mixing/alternating coordination, shared neighbor settings, with relevant modes and thresholds.

#### Methods
- `get_cutoff(type1: str, type2: str) -> Optional[float]`: Returns the cutoff distance between two node types or `None`.
- `__str__()`: Returns a formatted multi-line string detailing the clustering settings.


### AnalysisSettings

Contains flags to specify which analyzers and analysis metrics to calculate.

#### Attributes
- Boolean flags per analysis type, e.g., `with_average_cluster_size`, `with_largest_cluster_size`, `with_concentration`, `with_gyration_radius`, etc.
- `overwrite` (bool): Whether to overwrite existing output files. Default: `True`.
- `with_all` (bool): Enables all analysis flags at once. Default: `False`.

#### Methods
- `get_analyzers() -> List[str]`: Returns list of analyzer class names based on active flags.
- `__str__()`: Returns formatted string reporting enabled analyzers.


### LatticeSettings

Settings for specifying lattice parameters and application.

#### Attributes
- `lattice` (np.ndarray): Base lattice matrix. Default zero 3x3 matrix.
- `apply_custom_lattice` (bool): Whether to apply a user-defined lattice. Default: `False`.
- `custom_lattice` (np.ndarray): Custom lattice matrix. Default zero 3x3 matrix.
- `get_lattice_from_file` (bool): Whether to load lattice from a file. Default: `False`.
- `lattice_file_location` (str): Path to lattice file. Default: `"./"`.
- `apply_lattice_to_all_frames` (bool): Whether to apply lattice across all frames. Default: `True`.

#### Methods
- `__str__()`: Returns detailed string representation including lattice matrices.


### Settings

Composite settings class aggregating general, lattice, clustering, and analysis settings.

#### Attributes
- `project_name` (str): Project identifier. Default: `"default"`.
- `export_directory` (str): Base directory for exports. Default: `"export"`.
- `file_location` (str): Path of trajectory file. Default: `"./"`.
- `range_of_frames` (Tuple[int, int]): Frame processing range. Default: `(0, -1)`.
- `apply_pbc` (bool): Apply periodic boundary conditions. Default: `True`.
- `verbose`, `save_logs`, `save_performance` (bool): Flags for verbosity and saving outputs.
- Sub-settings as dataclass fields:
  - `general` (GeneralSettings)
  - `lattice` (LatticeSettings)
  - `clustering` (ClusteringSettings)
  - `analysis` (AnalysisSettings)

#### Properties
- `output_directory` (str): Combines export directory and project name.

#### Methods
- `set_range_of_frames(start: int, end: Optional[int] = None)`: Sets the frame range with validation.
- `__str__()`: Returns formatted string summarizing major settings (excluding general).

### SettingsBuilder

Builder pattern class for incrementally constructing a `Settings` instance with validation.

#### Methods
- `with_lattice(lattice: LatticeSettings) -> SettingsBuilder`: Sets lattice settings.
- `with_general(general: GeneralSettings) -> SettingsBuilder`: Sets general settings with validation.
- `with_analysis(analysis: AnalysisSettings) -> SettingsBuilder`: Sets analysis settings.
- `with_clustering(clustering: ClusteringSettings) -> SettingsBuilder`: Sets clustering settings with validations on criteria, connectivities, coordination number modes, flags consistency, and list lengths.
- `build() -> Settings`: Returns the validated, constructed `Settings` instance.

### Usage Example

```python
from nexus import SettingsBuilder, main
import nexus.config.settings as c

general = c.GeneralSettings(
    project_name="MyProject",
    export_directory="results",
    file_location="trajectory.xyz",
    range_of_frames=(0, 100),
    apply_pbc=True,
    verbose=True,
    save_logs=True,
    save_performance=True
)

lattice = c.LatticeSettings(
    apply_custom_lattice=True,
    custom_lattice = np.eye(3)*10.0
)

clustering = c.ClusteringSettings(
    criterion="distance",
    node_types=["A", "B"],
    node_masses=[1.0, 1.5],
    connectivity=["A-B", "B-A"],
    cutoffs=[c.Cutoff("A", "B", 3.5)],
    with_printed_unwrapped_clusters=True,
    print_mode="all"
)

analysis = c.AnalysisSettings(
    with_average_cluster_size=True,
    with_concentration=True,
    with_all=False
)

settings = (SettingsBuilder()
    .with_general(general)
    .with_lattice(lattice)
    .with_clustering(clustering)
    .with_analysis(analysis)
    .build()
)
```

***

This reference summarizes all the key settings classes and the builder for configuring the package's runtime behavior.

## `core` Module

### Node

Represents a node in the system, encapsulating properties like position, neighbors, mass, and coordination. Utilizes Python dataclasses with slots and ordering for efficient and comparable instances.

#### Attributes
- `symbol` (str): Symbolic identifier of the node (i.e., "1", "2", "A", "B", "Si", "O", etc.).
- `node_id` (int): Unique identifier (auto-incremented if not provided).
- `position` (np.ndarray): 3D coordinates of the node.
- `parent` (Optional[Node]): Reference to parent node; defaults to self if not set.
- `neighbors` (List[Node]): List of neighboring nodes.
- `cluster_id` (Optional[int]): Identifier of the cluster the node belongs to; defaults to its own node_id.
- `distances` (Optional[List[float]]): Distances to each neighbor.
- `indices` (Optional[List[int]]): Indices corresponding to neighbors.
- `mass` (Optional[float]): Mass of the node; defaults to 0.0 if not set.
- `coordination` (Optional[int]): Coordination number; defaults to 0.
- `other` (Optional[List[str]]): List for additional arbitrary attributes.

#### Initialization
- Auto-assigns default values for optional attributes after creation:
  - Default position as zero vector if none provided.
  - Auto-incrementation of `node_id` if none assigned.
  - Defaults for mass, coordination, parent, cluster_id, neighbors, and other attributes.

#### Methods

- `wrap_position(position: np.ndarray, lattice: np.ndarray) -> np.ndarray`
  - Static method.
  - Wraps a given position vector inside a periodic box defined by the lattice matrix.
  
- `add_neighbor(node: Node) -> None`
  - Adds a node to the neighbors list.
  
- `reset_parent() -> None`
  - Resets the parent attribute to self.
  
- `set_coordination(coordination: int) -> None`
  - Sets the coordination number.
  
- `__str__() -> str`
  - Returns a readable string representation showing node id, symbol, coordination, neighbor count, and position.
  
- `__repr__() -> str`
  - Same as `__str__` for convenience.

### Usage Example

In the typical workflow of Necus-CAT, the node are created automatically by the `Frame` class during the frame parsing.
The following lines are not recommended and are only for illustration purposes if needed.

```python
import numpy as np
from nexus.core.node import Node

position = np.array([1.0, 2.0, 3.0])
node = Node(symbol="A", node_id=None, position=position)

print(node)  # e.g., Node 0 (A) | coordination: 0 | neighbors: 0 | position: [1. 2. 3.]

# Add neighbor
neighbor_node = Node(symbol="B", node_id=None, position=np.array([4.0, 5.0, 6.0]))
node.add_neighbor(neighbor_node)

# Set coordination
node.set_coordination(1)

print(node)
```

Here's the API reference for the `Cluster` class in the core module:

### Cluster

Represents a cluster of nodes with connectivity, spatial properties, and periodic boundary handling.

#### Initialization
```python
Cluster(connectivity: str, root_id: int, size: int, settings: Settings, lattice: np.ndarray)
```
Creates a cluster with a connectivity label, root node ID, size, settings, and lattice matrix.

#### Attributes
- `nodes` (List[Node]): Nodes in the cluster.
- `connectivity` (str): Connectivity descriptor.
- `root_id` (int): ID of the root node.
- `size` (int): Number of nodes.
- `settings` (Settings): Configuration settings.
- `lattice` (np.ndarray): Lattice matrix.
- `_inv_lattice` (np.ndarray): Inverse of lattice matrix.
- `center_of_mass` (np.ndarray): Cluster's center of mass.
- `symbols` (list): Symbols of cluster nodes.
- `indices` (list): Node IDs.
- `unwrapped_positions` (np.ndarray): Node positions unwrapped in periodic space.
- `percolation_probability` (str): Directions in which cluster percolates (e.g., "xyz").
- `gyration_radius` (float): Radius of gyration.
- `order_parameter` (list): Order parameters for each dimension.
- `total_nodes` (int): Total nodes considered in the system.
- `concentration` (float): Fractional concentration.
- `is_percolating` (bool): True if cluster spans all directions.
- `is_spanning` (bool): True if cluster spans boundaries.
- `linkages` (List[Tuple[int, int]]): Node linkage pairs.
- `decoration_atoms` (Dict[int, Dict]): Additional nodes decorating the cluster, e.g., bridging atoms.

#### Methods

- `add_node(node: Node)`: Add a node to the cluster and assign its cluster ID.
- `set_lattice(lattice: np.ndarray)`: Update lattice and its inverse.
- `get_nodes() -> List[Node]`: Return cluster nodes.
- `get_connectivity() -> str`: Return connectivity descriptor.
- `get_size() -> int`: Return cluster size.
- `set_indices_and_positions(positions_dict: Dict[int, np.ndarray])`: Sets nodes’ indices and unwrapped positions.
- `calculate_center_of_mass()`: Calculate cluster center of mass; wraps it in lattice and translates decoration atoms.
- `calculate_gyration_radius()`: Compute radius of gyration.
- `calculate_percolation_probability()`: Identify lattice directions in which cluster percolates and update flags.
- `calculate_order_parameter()`: Calculate order parameters based on percolation and cluster size.
- `calculate_concentration()`: Compute concentration relative to total nodes.
- `_unwrap_vector(vector: np.ndarray) -> np.ndarray`: Apply periodic unwrapping to a vector.
- `calculate_unwrapped_positions()`: Unwrap cluster node positions relative to root; unwraps decoration atoms if applicable, with progress bar.
- `__str__()` / `__repr__()`: String representation summarizing cluster attributes.


### Usage Example

In the typical workflow of Necus-CAT, the clusters are created automatically by a `finder` during the clustering process.
The following lines are not recommended and are only for illustration purposes if needed.


```python
import numpy as np
from nexus.core.cluster import Cluster, Node
from nexus.config.settings import Settings

# Dummy settings
lattice = np.eye(3)
settings = Settings(
    ...
)

cluster = Cluster(connectivity="A-B", root_id=0, size=5, settings=settings, lattice=lattice)

node1 = Node(symbol="A", node_id=0, position=np.array([0.0, 0.0, 0.0]))
cluster.add_node(node1)

cluster.calculate_unwrapped_positions()
cluster.calculate_center_of_mass()

print(cluster)
```

### Frame

Represents a trajectory frame containing nodes, lattice information, and cluster data.

#### Attributes
- `frame_id` (int): Unique identifier of the frame.
- `nodes` (List[Node]): List of nodes contained in this frame.
- `lattice` (np.ndarray): 3x3 lattice matrix defining the periodic box.
- `_data` (Dict[str, np.ndarray]): Internal data dictionary holding node properties such as symbols and positions.
- `_settings` (Settings): Settings object with configuration options affecting the frame.
- `clusters` (Optional[List[Cluster]]): List of clusters detected in the frame.
- `connectivities` (Optional[List[str]]): List of connectivity descriptors associated with clusters.

#### Initialization and Validation
- `__post_init__()`: Ensures that `nodes` is a list and `lattice` is a 3x3 numpy array.
  
#### Methods

- `initialize_nodes() -> None`
  - Initializes the node list from `_data`, filtering by node types specified in the clustering settings.
  - Assigns node IDs sequentially starting at zero.
  
- `set_lattice(lattice: np.ndarray) -> None`
  - Sets the frame lattice matrix.
  - Validates that input lattice is a 3x3 non-singular numpy array.
  
- `get_lattice() -> Optional[np.ndarray]`
  - Returns the lattice matrix.

- `get_unique_elements() -> List[str]`
  - Returns a list of unique node symbols present in the frame.

- `get_node_by_id(node_id: int) -> Optional[Node]`
  - Finds and returns a node by its unique node ID, or None if not found.

- `get_positions() -> np.ndarray`
  - Returns an array of all node positions in the frame.

- `get_positions_by_element() -> Dict[str, np.ndarray]`
  - Groups node positions by their element symbol and returns a dictionary mapping symbol to position arrays.

- `get_wrapped_positions() -> np.ndarray`
  - Returns positions wrapped inside the periodic lattice box using utility wrapping function.

- `get_wrapped_positions_by_element() -> Dict[str, np.ndarray]`
  - Returns wrapped positions grouped by element symbol.

- `get_clusters() -> List[Cluster]`
  - Returns the list of clusters in the frame.

- `get_nodes() -> List[Node]`
  - Returns the list of nodes in the frame.

- `get_networking_nodes() -> int`
  - Returns the total number of nodes participating in clusters by summing cluster sizes.

- `get_connectivities() -> List[str]`
  - Returns the list of connectivities associated with the frame.

- `set_connectivities(connectivities: List[str]) -> None`
  - Sets the frame's connectivities.

- `add_cluster(cluster: Cluster) -> None`
  - Adds a cluster to the frame's cluster list, initializing the list if necessary.

- `set_clusters(clusters: List[Cluster]) -> None`
  - Sets the frame's clusters, assigns frame ID to clusters, and sets their lattice matrix.

- `get_concentration() -> Dict[str, float]`
  - Calculates and returns concentrations for each connectivity type based on cluster sizes relative to total nodes, ensuring all connectivities have a value.

- `__len__() -> int`
  - Returns the number of nodes in the frame.

- `__str__() -> str`
  - Returns a readable string summary of the frame including frame ID, node count, and cluster count.

- `__repr__() -> str`
  - Returns detailed string representation including first node data and lattice matrix.

- `__del__() -> None`
  - Cleans up references to nodes, clusters, lattice, data, and connectivities.

***

### Usage Example

Similarly to the `Node` and `Cluster` classes, the `Frame` objects are created automatically by the `System` object, and the following lines are not recommended and are only for illustration purposes if needed.

```python
import numpy as np
from nexus.core.frame import Frame

# Example data
frame_data = {
    "symbol": np.array(["A", "B", "A"]),
    "position": np.array([[0,0,0],[1,1,1],[2,2,2]])
}

settings = ...  # Assume a Settings object is initialized appropriately
frame = Frame(frame_id=0, nodes=[], lattice=np.eye(3), _data=frame_data, _settings=settings)

# Initialize nodes filtering by clustering node types from settings
frame.initialize_nodes()

print(frame)  # e.g., Frame 0 (num_nodes=..., num_clusters=None)
```

### System

Manages trajectory data for atomic systems, interfacing with file readers, loading frames, and iterating through trajectory frames based on settings.

#### Initialization

```python
System(reader: BaseReader, settings: Settings)
```

Constructs a `System` instance with a file reader and configuration settings.

##### Parameters
- `reader` (`BaseReader`): File reader object used to load trajectory data.
- `settings` (`Settings`): Configuration object with parameters like file location and frame ranges.

***

#### Attributes

- `reader` (`BaseReader`): The trajectory file reader.
- `settings` (`Settings`): Settings controlling system and file paths.
- `current_frame` (`Optional[Frame]`): Currently loaded frame or `None` if none loaded.
- `_current_frame_index` (`Optional[int]`): Internal index tracking current frame position.
- `_num_frames` (`Optional[int]`): Cached total number of frames in the trajectory.

***

#### Methods

- `load_frame(frame_index: int) -> bool`
  - Loads a specific frame by index.
  - Validates index against settings range.
  - Uses reader’s `parse` method to load the frame.
  - Updates `current_frame` and current index if successful.
  - Returns `True` if frame loaded; `False` otherwise.

- `get_frame(frame_index: int) -> Optional[Frame]`
  - Retrieves a frame by index, loading if necessary.
  - Returns the `Frame` object or `None` if loading fails.

- `get_num_frames() -> int`
  - Returns total number of frames in trajectory.
  - Uses cached value if available.
  - Falls back to counting frames by iteration if needed.
  - Returns `0` on error or if no frames detected.

- `iter_frames() -> Generator[Frame, None, None]`
  - Generator yielding frames within settings-defined frame range.
  - Uses reader’s frame indices if available.
  - Iterates frame-by-frame loading each one on demand.
  - Handles loading exceptions gracefully by continuing.

- `__iter__() -> System`
  - Makes `System` iterable, resetting internal frame index to start.

- `__next__() -> Frame`
  - Returns next frame in iteration sequence respecting frame range.
  - Raises `StopIteration` when range or data is exhausted.
  - Loads each frame on-demand using `load_frame()`.

***

### Usage Example

```python
from io_module import SomeReader
from core_module import System
from config_module import Settings

settings = Settings(file_location="trajectory.dat", range_of_frames=(0, 100))
reader = SomeReader()
system = System(reader, settings)

# Iterate through frames using generator
for frame in system.iter_frames():
    print(frame)

# Use as iterator
for frame in system:
    print(frame)
```

## `analysis` Module

### StrategyFactory

Factory class responsible for managing and providing clustering strategies based on given frame data and settings.

#### Initialization
```python
StrategyFactory(frame: Frame, settings: Settings)
```
Creates a `StrategyFactory` instance by initializing and registering available clustering strategies using the provided frame and settings.

##### Parameters
- `frame` (`Frame`): The frame of trajectory data on which clustering strategies operate.
- `settings` (`Settings`): Configuration settings to determine the appropriate strategy.

***

#### Attributes
- `_strategies` (dict): Internal mapping of strategy class names to instantiated strategy objects.

***

#### Methods

- `register_strategy(strategy: BaseClusteringStrategy) -> None`
  - Registers a clustering strategy instance in the factory.
  - Adds the strategy to the internal dictionary keyed by its class name.

- `get_strategy(settings: Settings) -> Optional[BaseClusteringStrategy]`
  - Returns the appropriate clustering strategy instance based on clustering configuration flags in `settings`.
  - Strategy selection logic:
    - If the clustering has coordination number or alternating cluster enabled (but not number of shared neighbors), returns the `CoordinationStrategy`.
    - If number of shared neighbors is enabled, returns the `SharedStrategy`.
    - If coordination features disabled and criterion is `"distance"`, returns the `DistanceStrategy`.
    - If coordination features disabled and criterion is `"bond"`, returns the `BondingStrategy`.
    - Returns `None` if no matching strategy is found.

***

