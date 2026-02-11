API Reference
=============

This section provides a detailed reference for the Nexus-CAT API. For a complete list of classes and methods, refer to the source code.

## main.py

### `main` function

Run the full analysis pipeline. Executes the sequential workflow: scan trajectory file, iterate frames, find neighbors, build clusters via union-find, run enabled analyzers, and write results to the export directory. Performance metrics are optionally recorded at each step.

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

General configuration parameters for the analysis pipeline.

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

Distance cutoff between two node types for neighbor searching.

#### Attributes

- `type1` (str): Symbol of the first node type.
- `type2` (str): Symbol of the second node type.
- `distance` (float): Maximum distance for this pair to be considered neighbors.

#### Methods

- `__str__()`: Returns formatted string representation.
- `get_distance() -> float`: Returns the cutoff distance.


### ClusteringSettings

Configuration parameters for the clustering algorithm. Controls which clustering strategy is selected, the neighbor search method, node filtering, cutoff distances, and optional coordination or shared-neighbor analysis modes.

#### Attributes

- `criterion` (str): Clustering criterion, either `"distance"` or `"bond"`. Default: `"distance"`.
- `neighbor_searcher` (str): Method for neighbor searching, e.g., `"kd_tree"`. Default: `"kd_tree"`.
- `node_types` (List[str]): List of node types.
- `node_masses` (List[float]): List of node masses in reduced units.
- `connectivity` (List[str]): Connectivity specifications.
- `cutoffs` (List[Cutoff]): List of cutoff distances.
- `with_printed_unwrapped_clusters` (bool): Flag to print unwrapped clusters. Default: `False`.
- `print_mode` (str): Print mode among `"all"`, `"connectivity"`, `"individual"`, or `"none"`. Default: `"none"`.
- `with_coordination_number` (bool): Flag to enable coordination number constraints. Default: `False`.
- `coordination_mode` (str): Mode for coordination number calculation, either: 
  - `"all_types"`: counts all the neighbors regardless of type.
  - `"same_type"`: counts only neighbors of the same type as the central node.
  - `"different_type"`: counts only neighbors of different types than the central node.
  - `<node_type>`: counts only neighbors of the specified node type, i.e., 'O', 'Si' etc. Default: `"all_types"`.
- `coordination_range` (List[int]): Minimum and maximum coordination numbers to consider. Default: `[]`. (e.g., `[4, 6]`)
- `with_pairwise` (bool): Flag to enable pairwise coordination number constraints. Default: `False`.
- `with_mixing` (bool): Flag to enable mixing coordination number constraints. Default: `False`.
- `with_alternating` (bool): Flag to enable alternating coordination number constraints. Default: `False`.
- `with_connectivity_name` (str): Name of the connectivity to apply coordination number constraints to. Default: `""`.
- `with_number_of_shared` (bool): Flag to enable shared neighbor constraints. Default: `False`.
- `shared_mode` (str): Mode for shared neighbor calculation, either:
  - `"all_types"`: counts all shared neighbors regardless of type.
  - `"same_type"`: counts only shared neighbors of the same type as the central node.
  - `"different_type"`: counts only shared neighbors of different types than the central node.
  - `<node_type>`: counts only shared neighbors of the specified node type, i.e., 'O', 'Si' etc. Default: `"all_types"`.
- `shared_threshold` (int): Minimum number of shared neighbors required. Default: `1`.
- `shared_threshold_mode` (str): Mode for shared neighbor threshold behavior, either:
  - `"exact"`: requires exactly the specified number of shared neighbors to account the connectivity between the nodes.
  - `"minimum"`: requires at least the specified number of shared neighbors to account the connectivity between the nodes.
  - `"maximum"`: requires at most the specified number of shared neighbors to account the connectivity between the nodes. Default: `"exact"`.

#### Methods

- `get_cutoff(type1: str, type2: str) -> Optional[float]`: Returns the cutoff distance between two node types or `None`.
- `__str__()`: Returns a formatted multi-line string detailing the clustering settings.


### AnalysisSettings

Configuration parameters controlling which analyzers are enabled. Each ``with_*`` flag enables the corresponding analyzer in the pipeline. Setting ``with_all`` enables every available analyzer at once.

#### Attributes

- Boolean flags per analysis type, e.g., `with_average_cluster_size`, `with_largest_cluster_size`, `with_concentration`, `with_gyration_radius`, etc.
- `overwrite` (bool): Whether to overwrite existing output files. Default: `True`.
- `with_all` (bool): Enables all analysis flags at once. Default: `False`.

#### Methods

- `get_analyzers() -> List[str]`: Returns list of analyzer class names based on active flags.
- `__str__()`: Returns formatted string reporting enabled analyzers.


### LatticeSettings

Configuration parameters for the simulation cell lattice. Controls whether a custom lattice is applied, the lattice matrix values, and whether the lattice is read from an external file.

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

Composite configuration object for the entire analysis pipeline. Aggregates general, lattice, clustering, and analysis sub-settings into a single object. Constructed via ``SettingsBuilder`` which validates constraints between fields.

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

Builder for constructing and validating a ``Settings`` object. Provides a fluent interface for setting sub-configurations. Each ``with_*`` method validates its input and returns ``self`` for chaining. Call ``build()`` to obtain the final ``Settings`` instance.

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


This reference summarizes all the key settings classes and the builder for configuring the package's runtime behavior.

## `core` Module

### Node

Dataclass representing a node in the simulation frame. A Node is the fundamental unit of the cluster analysis pipeline. It stores properties (symbol, position, mass) and participates in the union-find algorithm through its parent field. Nodes self-parent by default, meaning each node starts as the root of its own cluster.

#### Attributes

- `symbol` (str): Chemical symbol of the atom (e.g., "Si", "O").
- `node_id` (int): Unique identifier (auto-incremented if not provided).
- `position` (np.ndarray): 3D coordinates of the node.
- `parent` (Optional[Node]): Reference to parent node; defaults to self if not set.
- `neighbors` (List[Node]): List of neighboring nodes.
- `cluster_id` (Optional[int]): Identifier of the cluster the node belongs to; defaults to its own node_id.
- `distances` (Optional[List[float]]): Distances to each neighbor.
- `indices` (Optional[List[int]]): Indices corresponding to neighbors.
- `mass` (Optional[float]): Mass of the node; defaults to 1.0 if not set.
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

Here's the API reference for the `Cluster` class in the core module:

### Cluster

Represents a group of connected nodes forming a single cluster. A cluster is produced by a clustering strategy and consumed by analyzers. It holds the set of member nodes and provides methods to compute physical properties such as unwrapped positions, center of mass, gyration radius, percolation probability, and order parameter.

#### Initialization

```python
Cluster(connectivity: str, root_id: int, size: int, settings: Settings, lattice: np.ndarray)
```
Creates a cluster with a connectivity label, root node ID, size, settings, and lattice matrix.

##### Parameters

- `connectivity` (str): Connectivity descriptor (e.g., `"Si-O-Si"` or `"Si4-Si4"`).
- `root_id` (int): ID of the root node defining the cluster.
- `size` (int): Number of nodes in the cluster.
- `settings` (`Settings`): Configuration settings for the analysis.
- `lattice` (np.ndarray): 3×3 lattice matrix defining the periodic simulation box.


#### Attributes

**Basic Properties:**

- `nodes` (List[Node]): List of nodes belonging to the cluster.
- `connectivity` (str): Connectivity descriptor identifying cluster type.
- `root_id` (int): Unique identifier of the cluster's root node.
- `size` (int): Number of nodes in the cluster.
- `settings` (`Settings`): Configuration settings.
- `lattice` (np.ndarray): 3×3 lattice matrix.
- `symbols` (list): Element symbols of cluster nodes.
- `indices` (list): Global node IDs.
- `_inv_lattice` (np.ndarray): Inverse of the lattice matrix for coordinate transformations.
- `_all_connectivities` (Set[str]): Set of all connectivity types in the analysis.

**Structural Properties:**

- `center_of_mass` (np.ndarray): Wrapped center of mass position in Cartesian coordinates.
- `unwrapped_positions` (np.ndarray): Node positions unwrapped across periodic boundaries.
- `gyration_radius` (float): Radius of gyration quantifying spatial extent.

**Percolation Properties:**

- `percolation_probability` (str): Directions in which cluster percolates (e.g., `"x"`, `"xy"`, `"xyz"`).
- `order_parameter` (list): Order parameters `[P∞_x, P∞_y, P∞_z]` for each dimension.
- `is_percolating` (bool): `True` if cluster percolates in all three dimensions.
- `is_spanning` (bool): `True` if cluster is the largest in its connectivity type.

**Connectivity Data:**

- `linkages` (List[Tuple[int, int]]): List of node ID pairs representing bonds between networking nodes.
- `_linkage_set` (Set[Tuple[int, int]]): Internal set for efficient linkage tracking during unwrapping.
- `period_vectors` (List[np.ndarray]): Period vectors detected during BFS unwrapping, representing displacements across periodic boundaries. Empty for non-periodic clusters.
- `decoration_atoms` (Dict[int, Dict]): Dictionary of decorating nodes (e.g., bridging atoms in bond-based clustering). Keys are node IDs, values contain `symbol`, `position`, and `coordination`.

**Statistical Properties:**

- `total_nodes` (int): Total number of networking nodes in the system (for normalization).
- `concentration` (float): Fractional concentration: `size / total_nodes`.

***

#### Methods

- `add_node(node: Node) -> None`
  - Adds a node to the cluster and assigns the cluster ID to the node.
  - **Parameters**:
    - `node` (Node): Node object to add.
  - **Returns**: None

- `set_lattice(lattice: np.ndarray) -> None`
  - Updates the lattice matrix and its inverse.
  - **Parameters**:
    - `lattice` (np.ndarray): New 3×3 lattice matrix.
  - **Returns**: None

- `get_nodes() -> List[Node]`
  - Returns the list of nodes in the cluster.
  - **Returns**: List of `Node` objects.

- `get_connectivity() -> str`
  - Returns the connectivity descriptor string.
  - **Returns**: Connectivity descriptor.

- `get_size() -> int`
  - Returns the number of nodes in the cluster.
  - **Returns**: Cluster size (integer).

- `set_indices_and_positions(positions_dict: Dict[int, np.ndarray]) -> None`
  - Sets node indices, symbols, and unwrapped positions from a dictionary.
  - **Parameters**:
    - `positions_dict` (Dict[int, np.ndarray]): Mapping of node IDs to unwrapped positions.
  - **Returns**: None
  - **Behavior**: Populates `symbols`, `indices`, and `unwrapped_positions` attributes.

- `calculate_center_of_mass() -> None`
  - Calculates the cluster's center of mass in unwrapped coordinates.
  - **Returns**: None
  - **Behavior**:
    - Computes mean of unwrapped positions.
    - Wraps center of mass back into primary periodic cell.
    - Translates all unwrapped positions and decoration atoms to center the cluster around the wrapped COM.
    - Updates `center_of_mass` attribute.

- `calculate_gyration_radius() -> None`
  - Computes the radius of gyration about the center of mass.
  - **Returns**: None
  - **Formula**: $R_g = \sqrt{\frac{1}{N} \sum_{i=1}^{N} |\mathbf{r}_i - \mathbf{r}_{\text{COM}}|^2}$
  - **Behavior**: Returns 0.0 for clusters with size ≤ 1.

- `calculate_percolation_probability() -> None`
  - Determines percolation using the period vector algorithm.
  - **Returns**: None
  - **Method**: Period vector detection following Livraghi et al. (2021) [J. Chem. Theory Comput.](https://pubs.acs.org/doi/10.1021/acs.jctc.1c00423)
  - **Behavior**:
    - Uses pre-computed `period_vectors` collected during `calculate_unwrapped_positions()`, avoiding a redundant BFS traversal.
    - Calculates percolation dimension from linear independence of period vectors.
    - Updates `percolation_probability` (e.g., `"xyz"`) and `is_percolating` flag.
  - **Percolation Theory Significance**: True percolation requires the cluster to connect to itself across periodic boundaries through actual bond connectivity, not just spatial extent. This method distinguishes genuine percolating clusters from large finite clusters.

- `_calculate_period_dimension() -> int`
  - Calculates the algebraic dimension of the period vector set using singular value decomposition.
  - **Returns**: Integer (0, 1, 2, or 3) indicating percolation dimensionality.
  - **Method**: Reads from `self.period_vectors`. Computes rank of period vector matrix via SVD with numerical tolerance.

- `_get_percolation_directions() -> str`
  - Determines which Cartesian directions (x, y, z) the cluster percolates in.
  - **Returns**: String containing percolating directions (e.g., `"xy"` for percolation in x and y).
  - **Behavior**: Reads from `self.period_vectors`. Checks which directions have significant fractional components (> 0.5 lattice units) in linearly independent period vectors.

- `_get_independent_periods(period_vectors: List[np.ndarray]) -> List[np.ndarray]`
  - Extracts linearly independent period vectors.
  - **Parameters**:
    - `period_vectors` (List[np.ndarray]): All detected period vectors.
  - **Returns**: List of up to 3 linearly independent period vectors.
  - **Method**: Uses incremental matrix rank calculation to select independent vectors.

- `calculate_order_parameter() -> None`
  - Calculates order parameters for each dimension based on percolation status.
  - **Returns**: None
  - **Formula**: $P_\infty = \frac{N_{\text{cluster}}}{N_{\text{total}}}$ for percolating dimensions, 0 otherwise.
  - **Behavior**: Assigns `[P∞, 0, 0]` for 1D percolation, `[P∞, P∞, 0]` for 2D, and `[P∞, P∞, P∞]` for 3D.

- `calculate_concentration() -> None`
  - Computes the cluster concentration relative to total networking nodes.
  - **Returns**: None
  - **Formula**: $c = \frac{N_{\text{cluster}}}{N_{\text{total}}}$
  - **Behavior**: Updates `concentration` attribute.

- `_unwrap_vector(vector: np.ndarray) -> np.ndarray`
  - Applies minimum image convention to unwrap a displacement vector.
  - **Parameters**:
    - `vector` (np.ndarray): Displacement vector to unwrap.
  - **Returns**: Unwrapped vector in Cartesian coordinates.
  - **Method**: Converts to fractional coordinates, subtracts nearest integer, converts back to Cartesian.

- `calculate_unwrapped_positions() -> None`
  - Unwraps all cluster node positions relative to the root node using breadth-first traversal, and simultaneously detects period vectors for percolation analysis.
  - **Returns**: None
  - **Behavior**:
    - Performs a single BFS starting from root node.
    - For each neighbor, calculates unwrapped position by adding minimum image displacement to current position.
    - When encountering a previously visited node through a different path, records the displacement as a period vector.
    - Records linkages between nodes.
    - For bond-based clustering, identifies and unwraps decorating atoms (e.g., bridging nodes).
    - Displays progress bar if verbose mode enabled.
    - Handles both `"distance"` and `"bond"` clustering criteria.
  - **Output**: Populates `unwrapped_positions`, `linkages`, `period_vectors`, and `decoration_atoms` attributes.

- `__str__() -> str`
  - Returns a concise string representation of the cluster.
  - **Format**: `"{root_id} {connectivity} {size} {is_percolating} {node_ids}"`
  - **Returns**: Human-readable cluster summary (truncated to 20 node IDs).

- `__repr__() -> str`
  - Returns the same representation as `__str__()`.

***

### Frame

Represents a single snapshot of a trajectory. A frame holds the raw node data, the simulation cell lattice, and the clusters produced by a clustering strategy. It is created by a reader, populated during the analysis pipeline, and consumed by analyzers.

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

### System

Manages trajectory data and provides frame-level access through a file reader. Wraps a reader with lazy frame iteration. On initialization it configures the reader's filename from settings and triggers a file scan to index frame offsets. Frames can then be accessed individually or iterated over as a generator.

#### Initialization

```python
System(reader: BaseReader, settings: Settings)
```

Constructs a `System` instance with a file reader and configuration settings.

##### Parameters

- `reader` (`BaseReader`): File reader object used to load trajectory data.
- `settings` (`Settings`): Configuration object with parameters like file location and frame ranges.


#### Attributes

- `reader` (`BaseReader`): The trajectory file reader.
- `settings` (`Settings`): Settings controlling system and file paths.
- `current_frame` (`Optional[Frame]`): Currently loaded frame or `None` if none loaded.
- `_current_frame_index` (`Optional[int]`): Internal index tracking current frame position.
- `_num_frames` (`Optional[int]`): Cached total number of frames in the trajectory.


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

## `analysis` Module

### StrategyFactory

Factory class responsible for creating, registering, and retrieving clustering strategy instances. Implements the factory design pattern to manage a collection of clustering strategies. Selects the appropriate one based on the clustering settings flags.

#### Initialization

```python
StrategyFactory(frame: Frame, settings: Settings)
```
Creates a `StrategyFactory` instance by initializing and registering available clustering strategies using the provided frame and settings.

##### Parameters

- `frame` (`Frame`): The frame of trajectory data on which clustering strategies operate.
- `settings` (`Settings`): Configuration settings to determine the appropriate strategy.


#### Attributes

- `_strategies` (dict): Internal mapping of strategy class names to instantiated strategy objects.


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


### DistanceStrategy

Clustering strategy based on a direct distance criterion. Connects any two nodes of specified types that are within a cutoff distance of each other. Requires a 2-element connectivity specification (e.g., ``["Si", "Si"]``).

#### Inheritance

Inherits from `BaseClusteringStrategy`.

#### Initialization

```python
DistanceStrategy(frame: Frame, settings: Settings)
```
Creates a distance-based clustering strategy for the given frame.

##### Parameters

- `frame` (`Frame`): The trajectory frame containing nodes to be clustered.
- `settings` (`Settings`): Configuration settings specifying clustering parameters, cutoffs, and node types.


#### Attributes

- `frame` (`Frame`): The frame being processed.
- `clusters` (`List[Cluster]`): List of identified clusters.
- `_lattice` (`np.ndarray`): Lattice matrix from the frame.
- `_nodes` (`List[Node]`): List of nodes in the frame.
- `_settings` (`Settings`): Configuration settings.
- `_counter` (int): Counter tracking number of clusters formed.
- `_neighbor_searcher` (`NeighborSearcher`): Helper object for finding neighbors based on distance.


#### Methods

- `find_neighbors() -> None`
  - Executes neighbor searching algorithm using the internal neighbor searcher.
  - Identifies all nodes within cutoff distances of each other.

- `get_connectivities() -> List[str]`
  - Returns list of connectivity descriptors for the clustering.
  - Validates that connectivity configuration is a two-element list.
  - Returns formatted connectivity string (e.g., `["A-B"]`).
  - Raises `ValueError` if connectivity format is invalid.

- `build_clusters() -> List[Cluster]`
  - Constructs clusters by applying union-find algorithm on nodes within distance cutoffs.
  - Process flow:
    1. Filters networking nodes based on configured node types.
    2. Iterates through nodes, connecting neighbors matching connectivity criteria using union operation.
    3. Groups nodes by their root to form clusters.
    4. Creates `Cluster` objects for groups with more than one node.
    5. Calculates cluster properties: unwrapped positions, center of mass, gyration radius, percolation probability, concentration, and order parameter.
  - Uses progress bars (tqdm) for visual feedback if verbosity enabled.
  - Returns list of formed clusters with computed properties.


#### Usage Example

```python
from nexus import SettingsBuilder, main
import nexus.config.settings as c

general_settings = c.GeneralSettings(
  ...
)
lattice_settings = c.LatticeSettings(
  ...
)
analysis_settings = c.AnalysisSettings(
  ...
)

clustering_settings = c.ClusteringSettings(
    criterion="distance",
    node_types=["A", "B"],
    connectivity=["A", "B"],
    cutoffs=[c.Cutoff("A", "B", 3.5)]
    with_connectivity_name='A-B',
)

settings = (
    SettingsBuilder()
    .with_general(config_general)  # General settings \
    .with_lattice(config_lattice)  # Lattice settings \
    .with_clustering(config_clustering)  # Clustering settings \
    .with_analysis(config_analysis)  # Analysis settings \
    .build()  # Don't forget to build the settings object
)

main(settings)

```

### BondingStrategy

Clustering strategy that connects networking nodes via a bridging node. Uses a 3-element connectivity pattern (e.g., ``["Si", "O", "Si"]``) to form clusters by linking networking nodes that share a common bridging node of the specified type.

#### Inheritance
Inherits from `BaseClusteringStrategy`.

#### Initialization

```python
BondingStrategy(frame: Frame, settings: Settings)
```
Creates a bond-based clustering strategy for the given frame.

##### Parameters

- `frame` (`Frame`): The trajectory frame containing nodes to be clustered.
- `settings` (`Settings`): Configuration settings specifying clustering parameters, cutoffs, and node types.


#### Attributes

- `frame` (`Frame`): The frame being processed.
- `clusters` (`List[Cluster]`): List of identified clusters.
- `_lattice` (`np.ndarray`): Lattice matrix from the frame.
- `_nodes` (`List[Node]`): List of nodes in the frame.
- `_settings` (`Settings`): Configuration settings.
- `_counter` (int): Counter tracking number of clusters formed.
- `_neighbor_searcher` (`NeighborSearcher`): Helper object for finding neighbors based on distance.


#### Methods

- `find_neighbors() -> None`
  - Executes neighbor searching algorithm using the internal neighbor searcher.
  - Identifies all nodes within cutoff distances of each other.

- `get_connectivities() -> List[str]`
  - Returns list of connectivity descriptors for the clustering.
  - Validates that connectivity configuration is a three-element list.
  - Returns formatted connectivity string (e.g., `["Si-O-Si"]`).
  - Raises `ValueError` if connectivity format is invalid.

- `build_clusters() -> List[Cluster]`
  - Constructs clusters by identifying nodes connected through bridging atoms.
  - Process flow:
    1. Filters networking nodes based on configured node types, excluding bridge atoms.
    2. Iterates through nodes of the first type, examining their neighbors for bridge atoms.
    3. For each bridge atom, checks its neighbors for nodes of the third type.
    4. Connects matching node pairs using union operation.
    5. Groups nodes by their root to form clusters.
    6. Creates `Cluster` objects for groups with more than one node.
    7. Calculates cluster properties: unwrapped positions, center of mass, gyration radius, percolation probability, concentration, and order parameter.
  - Uses progress bars (tqdm) for visual feedback if verbosity enabled.
  - Returns list of formed clusters with computed properties.


#### Usage Example

```python
from nexus import SettingsBuilder, main
import nexus.config.settings as c

general_settings = c.GeneralSettings(
  ...
)

lattice_settings = c.LatticeSettings(
  ...
)

analysis_settings = c.AnalysisSettings(
  ...
)

clustering_settings = c.ClusteringSettings(
    criterion="bond",
    node_types=["A", "B", "C"],
    connectivity=["A", "B", "C"],
    cutoffs=[
        c.Cutoff("A", "B", 2.0),
        c.Cutoff("B", "C", 2.0),
        c.Cutoff("C", "A", 2.0)
    ]
)

settings = (
    SettingsBuilder()
    .with_general(general_settings)
    .with_lattice(lattice_settings)
    .with_clustering(clustering_settings)
    .with_analysis(analysis_settings)
    .build()
)

main(settings)
```

### CoordinationStrategy

Clustering strategy that groups nodes by coordination number. Extends the bonding pattern by constraining the coordination number of networking nodes. Supports several pairing modes: pairwise, mixing, alternating, and default.

#### Inheritance

Inherits from `BaseClusteringStrategy`.

#### Initialization

```python
CoordinationStrategy(frame: Frame, settings: Settings)
```
Creates a coordination-based clustering strategy for the given frame.

##### Parameters

- `frame` (`Frame`): The trajectory frame containing nodes to be clustered.
- `settings` (`Settings`): Configuration settings specifying clustering parameters, coordination ranges, and node types.


#### Attributes

- `frame` (`Frame`): The frame being processed.
- `clusters` (`List[Cluster]`): List of identified clusters.
- `_lattice` (`np.ndarray`): Lattice matrix from the frame.
- `_nodes` (`List[Node]`): List of nodes in the frame.
- `_settings` (`Settings`): Configuration settings.
- `_counter` (int): Counter tracking number of clusters formed.
- `_search_mode` (str): Current search mode (`"default"`, `"pairwise"`, `"mixing"`, or `"alternating"`).
- `_neighbor_searcher` (`NeighborSearcher`): Helper object for finding neighbors based on distance.


#### Methods

- `find_neighbors() -> None`
  - Executes neighbor searching algorithm using the internal neighbor searcher.
  - Calculates coordination number for all nodes after neighbor identification.

- `calculate_coordination(idx: int) -> None`
  - Calculates and sets the coordination number for a node at the specified index.
  - Coordination modes:
    - `"all_types"`: Counts all neighboring nodes.
    - `"same_type"`: Counts only neighbors with the same element symbol.
    - `"different_type"`: Counts only neighbors with different element symbols.
    - `"<node_type>"`: Counts only neighbors matching specified node type.

- `find(node: Node) -> Node`
  - Finds the root node in the union-find structure with path compression.
  - Returns the root parent node.

- `union(node_1: Node, node_2: Node) -> None`
  - Merges two nodes into the same cluster by connecting their roots.

- `get_connectivities() -> List[str]`
  - Generates connectivity descriptors based on clustering criterion and coordination settings.
  - Supports multiple coordination modes:
    - **Pairwise**: Nodes with identical coordination numbers (e.g., `["Si4-Si4", "Si5-Si5"]`).
    - **Mixing**: All combinations of coordination numbers within range (e.g., `["Si4-Si4", "Si4-Si5", "Si5-Si5"]`).
    - **Alternating**: Consecutive coordination pairs (e.g., `["Si4-Si4", "Si4-Si5", "Si5-Si5", "Si5-Si6"]`).
    - **Default**: Single user-specified connectivity name.
  - Returns list of formatted connectivity strings.

- `build_clusters() -> List[Cluster]`
  - Constructs clusters by grouping nodes with matching coordination criteria.
  - Process flow:
    1. Filters networking nodes based on node types and criterion.
    2. Generates connectivities based on coordination number range and mode.
    3. For each connectivity, finds and connects matching node pairs using union-find algorithm.
    4. Groups nodes by their root to form clusters.
    5. Creates `Cluster` objects for groups with more than one node.
    6. Calculates cluster properties: unwrapped positions, center of mass, gyration radius, percolation probability, concentration, and order parameter.
  - Uses progress bars (tqdm) for visual feedback if verbosity enabled.
  - Returns list of formed clusters with computed properties.

- `_find_cluster(networking_nodes: List[Node], connectivity: str, z1: int, z2: int) -> None`
  - Internal helper method to find clusters for a specific connectivity pattern.
  - Parameters:
    - `networking_nodes`: List of nodes to cluster.
    - `connectivity`: Connectivity descriptor string.
    - `z1`: First coordination number constraint.
    - `z2`: Second coordination number constraint.
  - Applies union operations on node pairs matching coordination criteria.
  - Creates and stores resulting clusters.


#### Usage Examples

The following examples use the `CoordinationStrategy` to cluster SiOz polyhedra in a SiO2 system with various approaches.

##### Example 1

The first example clusters SiOz polyhedra with a pairwise coordination number constraint in the interval [4, 6].
It will find the following connectivities: `SiO4-SiO4`, `SiO5-SiO5`, `SiO6-SiO6`.
Note that if the criterion is `bond`, the connectivity must be specified with 3 elements, the first and last being the networking nodes while the one in the middle being the 'bridging' node, e.g., ["Si", "O", "Si"] for SiOz-SiOz with `z` the coordination number.
In the case of `distance`, the connectivity must be specified with 2 elements, the first and last being the networking nodes, e.g., ["Si", "Si"] for Siz-Siz with `z` the coordination number.

```python
from nexus import SettingsBuilder, main
import nexus.config.settings as c

general_settings = c.GeneralSettings(
    ...
)

lattice_settings = c.LatticeSettings(
    ...
)

analysis_settings = c.AnalysisSettings(
    ...
)

clustering_settings = c.ClusteringSettings(
    criterion="bond",
    node_types=["Si", "O"],
    connectivity=["Si", "O", "Si"],
    cutoffs=[
        c.Cutoff("Si", "O", 2.3),
        c.Cutoff("Si", "Si", 3.5),
        c.Cutoff("O", "O", 3.05),
    ],
    with_coordination_number=True,
    coordination_mode="different_type",
    coordination_range=[4, 6],
    with_pairwise=True
)

# OR

clustering_settings = c.ClusteringSettings(
    criterion="distance",
    node_types=["Si", "Si"],
    connectivity=["Si", "Si"],
    cutoffs=[
        c.Cutoff("Si", "Si", 3.5),
    ],
    with_coordination_number=True,
    coordination_mode="same_type", # "all_types" or "same_type" or "different_type" or "<node_type>"
    coordination_range=[4, 6],
    with_pairwise=True
)

settings = (
    SettingsBuilder()
    .with_general(general_settings)
    .with_lattice(lattice_settings)
    .with_clustering(clustering_settings)
    .with_analysis(analysis_settings)
    .build()
)

main(settings)
```

##### Example 2

The second example clusters SiOz polyhedra with a mixing coordination number constraint in the interval [4, 6].
It will find the following connectivities: `SiO4-SiO4`, `SiO4-SiO5`, `SiO5-SiO5`, `SiO4-SiO6`, `SiO5-SiO6`, `SiO6-SiO6`.

```python
...

clustering_settings = c.ClusteringSettings(
    criterion="bond",
    node_types=["Si", "O"],
    connectivity=["Si", "O", "Si"],
    cutoffs=[
        c.Cutoff("Si", "O", 2.3),
        c.Cutoff("Si", "Si", 3.5),
        c.Cutoff("O", "O", 3.05),
    ],
    with_coordination_number=True,
    coordination_mode="different_type",
    coordination_range=[4, 6],
    with_pairwise=False,
    with_mixing=True
)

...
```

##### Example 3

The third example clusters SiOz polyhedra with an alternating coordination number constraint in the interval [4, 6].
It will find the following connectivities: `SiO4-SiO4`, `SiO4-SiO5`, `SiO5-SiO5`, `SiO4-SiO6`, `SiO5-SiO6`, `SiO6-SiO6`.

```python
...

clustering_settings = c.ClusteringSettings(
    criterion="bond",
    node_types=["Si", "O"],
    connectivity=["Si", "O", "Si"],
    cutoffs=[
        c.Cutoff("Si", "O", 2.3),
        c.Cutoff("Si", "Si", 3.5),
        c.Cutoff("O", "O", 3.05),
    ],
    with_coordination_number=True,
    coordination_mode="different_type",
    coordination_range=[4, 6],
    with_pairwise=False,
    with_mixing=False,
    with_alternating=True
)

...
```

##### Example 4

The fourth example clusters all SiOz polyhedra with a coordination number constraint in the interval [5, 7].
In this settings, a custom name to the connectivity can be given, e.g. `HD`.

```python
...

clustering_settings = c.ClusteringSettings(
    criterion="bond",
    node_types=["Si", "O"],
    connectivity=["Si", "O", "Si"],
    cutoffs=[
        c.Cutoff("Si", "O", 2.3),
        c.Cutoff("Si", "Si", 3.5),
        c.Cutoff("O", "O", 3.05),
    ],
    with_coordination_number=True,
    coordination_mode="different_type",
    coordination_range=[5, 7],
    with_pairwise=False,
    with_mixing=False,
    with_alternating=False,
    with_number_of_shared=False,
    with_connectivity_name="HD"
)

...
```

### SharedStrategy

Clustering strategy based on a minimum or exact number of shared neighbors. Extends the coordination strategy by additionally requiring that two networking nodes share at least or exactly a threshold number of common bridging neighbors. This distinguishes polyhedral linkage types such as corner-sharing, edge-sharing, or face-sharing.

#### Inheritance

Inherits from `BaseClusteringStrategy`.

#### Initialization

```python
SharedStrategy(frame: Frame, settings: Settings)
```
Creates a shared-neighbor-based clustering strategy for the given frame.

##### Parameters

- `frame` (`Frame`): The trajectory frame containing nodes to be clustered.
- `settings` (`Settings`): Configuration settings specifying clustering parameters, coordination ranges, shared thresholds, and node types.


#### Attributes

- `frame` (`Frame`): The frame being processed.
- `clusters` (`List[Cluster]`): List of identified clusters.
- `_lattice` (`np.ndarray`): Lattice matrix from the frame.
- `_nodes` (`List[Node]`): List of nodes in the frame.
- `_settings` (`Settings`): Configuration settings.
- `_counter` (int): Counter tracking number of clusters formed.
- `_neighbor_searcher` (`NeighborSearcher`): Helper object for finding neighbors based on distance.


#### Methods

- `find_neighbors() -> None`
  - Executes neighbor searching algorithm using the internal neighbor searcher.
  - Calculates coordination number for all nodes after neighbor identification.

- `calculate_coordination(idx: int) -> None`
  - Calculates and sets the coordination number for a node at the specified index.
  - Coordination modes:
    - `"all_types"`: Counts all neighboring nodes.
    - `"same_type"`: Counts only neighbors with the same element symbol.
    - `"different_type"`: Counts only neighbors with different element symbols.
    - `"<node_type>"`: Counts only neighbors matching specified node type.

- `get_number_of_shared(node_1: Node, node_2: Node) -> int`
  - Calculates the number of shared neighbors between two nodes.
  - Respects the shared mode from settings:
    - `"all_types"`: Counts all shared neighbors.
    - `"same_type"`: Counts shared neighbors of the same type.
    - `"different_type"`: Counts shared neighbors of different types.
    - `"<node_type>"`: Counts shared neighbors matching specified node type.
  - Returns count of shared neighbors.

- `find(node: Node) -> Node`
  - Finds the root node in the union-find structure with path compression.
  - Returns the root parent node.

- `union(node_1: Node, node_2: Node) -> None`
  - Merges two nodes into the same cluster by connecting their roots.

- `get_connectivities() -> List[str]`
  - Generates connectivity descriptors based on clustering criterion and coordination settings.
  - Uses `=` separator to indicate shared neighbor requirement (e.g., `"Si4=Si4"`).
  - Supports multiple coordination modes:
    - **Pairwise**: Nodes with identical coordination numbers (e.g., `["Si4=Si4", "Si5=Si5"]`).
    - **Mixing**: All combinations of coordination numbers within range (e.g., `["Si4=Si4", "Si4=Si5", "Si5=Si5"]`).
    - **Alternating**: Consecutive coordination pairs (e.g., `["Si4=Si4", "Si4=Si5", "Si5=Si5", "Si5=Si6"]`).
    - **Default**: Single user-specified connectivity name.
  - Returns list of formatted connectivity strings.

- `build_clusters() -> List[Cluster]`
  - Constructs clusters by grouping nodes with matching coordination and shared neighbor criteria.
  - Process flow:
    1. Filters networking nodes based on node types and criterion.
    2. Generates connectivities based on coordination number range and mode.
    3. For each connectivity, finds and connects matching node pairs with sufficient shared neighbors.
    4. Groups nodes by their root to form clusters.
    5. Creates `Cluster` objects for groups with more than one node.
    6. Calculates cluster properties: unwrapped positions, center of mass, gyration radius, percolation probability, concentration, and order parameter.
  - Uses progress bars (tqdm) for visual feedback if verbosity enabled.
  - Returns list of formed clusters with computed properties.

- `_find_cluster(networking_nodes: List[Node], connectivity: str, z1: int, z2: int) -> None`
  - Internal helper method to find clusters for a specific connectivity pattern.
  - Parameters:
    - `networking_nodes`: List of nodes to cluster.
    - `connectivity`: Connectivity descriptor string.
    - `z1`: First coordination number constraint.
    - `z2`: Second coordination number constraint.
  - Applies union operations on node pairs matching coordination and shared neighbor threshold.
  - Only connects nodes if the number of shared neighbors meets or exceeds the configured threshold.
  - Creates and stores resulting clusters.

### Usage Example

```python
from nexus import SettingsBuilder, main
import nexus.config.settings as c

general_settings = c.GeneralSettings(
  ...
)

lattice_settings = c.LatticeSettings(
  ...
)

analysis_settings = c.AnalysisSettings(
  ...
)

clustering_settings = c.ClusteringSettings(
    criterion="bond",
    node_types=["Si", "O"],
    connectivity=["Si", "O", "Si"],
    cutoffs=[
        c.Cutoff("Si", "O", 2.3),
        c.Cutoff("O", "O", 3.05),
        c.Cutoff("Si", "Si", 3.5)
    ],
    with_coordination_number=True,
    coordination_mode="different_type",
    coordination_range=[6, 6],
    with_number_of_shared=True,
    shared_mode="different_type",
    shared_threshold=2,
    with_pairwise=False,
    with_mixing=False,
    with_alternating=False,
    with_connectivity_name="Stishovite"
)

settings = (
    SettingsBuilder()
    .with_general(general_settings)
    .with_lattice(lattice_settings)
    .with_clustering(clustering_settings)
    .with_analysis(analysis_settings)
    .build()
)

main(settings)
```

### AnalyzerFactory

Factory class responsible for creating, registering, and retrieving analyzer instances. Implements the factory design pattern to manage a collection of analyzer objects.

#### Initialization

```python
AnalyzerFactory(settings: Settings, verbose: bool = True)
```
Creates an `AnalyzerFactory` instance by initializing and registering all available analyzer types.

##### Parameters

- `settings` (`Settings`): Configuration settings used to initialize analyzers with appropriate parameters.
- `verbose` (bool): Flag controlling verbosity of output. Default: `True`.


#### Attributes

- `_analyzers` (dict): Internal mapping of analyzer class names to instantiated analyzer objects.


#### Methods

- `register_analyzer(analyzer: BaseAnalyzer) -> None`
  - Registers an analyzer instance in the factory.
  - Adds the analyzer to the internal dictionary keyed by its class name.

- `get_analyzer(analyzer_name: str) -> Optional[BaseAnalyzer]`
  - Returns the analyzer instance corresponding to the specified name.
  - Returns `None` if no analyzer with the given name is registered.
  - Available analyzer names:
    - `"AverageClusterSizeAnalyzer"`: Calculates average cluster sizes.
    - `"ConcentrationAnalyzer"`: Computes cluster concentrations.
    - `"LargestClusterSizeAnalyzer"`: Identifies largest cluster sizes.
    - `"SpanningClusterSizeAnalyzer"`: Measures spanning cluster properties.
    - `"PercolationProbabilityAnalyzer"`: Determines percolation probabilities.
    - `"OrderParameterAnalyzer"`: Calculates order parameters.
    - `"ClusterSizeDistributionAnalyzer"`: Generates cluster size distributions.
    - `"GyrationRadiusAnalyzer"`: Computes radii of gyration.
    - `"CorrelationLengthAnalyzer"`: Calculates correlation lengths.

### AverageClusterSizeAnalyzer

Computes the weight-average cluster size ``<S>``. Uses the formula ``<S> = sum(s^2 * n(s)) / sum(s * n(s))`` where *s* is the cluster size and *n(s)* the number of clusters of size *s*. Percolating clusters are excluded to focus on the finite-cluster distribution.

#### Inheritance

Inherits from `BaseAnalyzer`.

#### Initialization

```python
AverageClusterSizeAnalyzer(settings: Settings)
```
Creates an analyzer for computing weight-average cluster sizes across trajectory frames.

##### Parameters

- `settings` (`Settings`): Configuration settings specifying analysis parameters and output paths.


#### Attributes

- `_raw_average_sizes` (`Dict[str, List[float]]`): Raw per-frame average cluster sizes for each connectivity type.
- `_raw_concentrations` (`Dict[str, List[float]]`): Raw per-frame concentration values for each connectivity.
- `average_sizes` (`Dict[str, float]`): Final mean average cluster size for each connectivity type.
- `std` (`Dict[str, float]`): Standard deviation of average cluster sizes across frames.
- `error` (`Dict[str, float]`): Standard error of average cluster sizes across frames.
- `concentrations` (`Dict[str, float]`): Mean concentration for each connectivity type.
- `_finalized` (bool): Flag indicating whether final calculations have been performed.


#### Methods

- `analyze(frame: Frame, connectivities: List[str]) -> None`
  - Analyzes a single frame to compute the weight-average cluster size for each connectivity type.
  - **Calculation Methodology**:
    - For each connectivity, extracts all non-percolating clusters.
    - Computes the weight-average size using the formula:
      $$
      \langle S \rangle = \frac{\sum_s s^2 n_s}{\sum_s s n_s}
      $$
      where $s$ is the cluster size and $n_s$ is the number of clusters of size $s$.
    - Stores raw values for later aggregation across all frames.
  - **Percolation Theory Significance**:
    - The weight-average cluster size is a second moment of the cluster size distribution.
    - This quantity is proportional to the percolation susceptibility $\chi$, which measures the system's response to changes in connectivity.
    - Near the percolation threshold $p_c$, $\langle S \rangle$ diverges as a power law: $\langle S \rangle \sim |p - p_c|^{-\gamma}$, where $\gamma$ is a critical exponent.
    - Percolating clusters are excluded because they represent infinite clusters in the thermodynamic limit and would dominate the calculation.

- `finalize() -> Dict[str, Dict[str, float]]`
  - Computes final statistics by averaging over all processed frames.
  - **Calculation Details**:
    - Calculates mean average cluster size across all frames.
    - Computes standard deviation and standard error using Bessel's correction (ddof=1).

- `get_result() -> Dict[str, Dict[str, float]]`
  - Returns the finalized analysis results without recalculation.
  - Provides access to concentrations, average cluster sizes, standard deviations, and errors for all connectivity types.

- `print_to_file() -> None`
  - Writes finalized results to a CSV-formatted data file.
  - Output file: `average_cluster_size.dat` in the export directory.
  - Includes header with metadata: date, number of frames analyzed, and column descriptions.
  - Format: `Connectivity_type,Concentration,Average_cluster_size,Standard_deviation,Standard_error`
  - Automatically removes duplicate lines to prevent data redundancy.

- `_write_header() -> None`
  - Internal method to write the file header with analysis metadata.
  - Respects overwrite settings: creates new file or appends to existing based on configuration.
  - Includes timestamp and frame count for reproducibility.

### ClusterSizeDistributionAnalyzer

Computes the cluster size distribution ``n(s)`` for each connectivity type. Counts how many clusters of each size *s* exist per connectivity across all frames. Percolating clusters are excluded from the distribution. Results are written to one file per connectivity.

#### Inheritance

Inherits from `BaseAnalyzer`.

#### Initialization

```python
ClusterSizeDistributionAnalyzer(settings: Settings)
```
Creates an analyzer for computing cluster size distributions across trajectory frames.

##### Parameters

- `settings` (`Settings`): Configuration settings specifying analysis parameters and output paths.


#### Attributes

- `_raw_size_distribution` (`Dict[str, Dict[int, List[int]]]`): Raw per-frame counts of clusters for each size and connectivity type.
- `_raw_concentrations` (`Dict[str, List[float]]`): Raw per-frame concentration values for each connectivity.
- `size_distribution` (`Dict[str, Dict[int, float]]`): Final total number of clusters count per size for each connectivity type.
- `std` (`Dict[str, Dict[int, float]]`): Standard deviation of cluster counts for each size across frames.
- `concentrations` (`Dict[str, float]`): Mean concentration for each connectivity type.
- `_finalized` (bool): Flag indicating whether final calculations have been performed.


#### Methods

- `analyze(frame: Frame, connectivities: List[str]) -> None`
  - Analyzes a single frame to compute the cluster size distribution for each connectivity type.
  - **Calculation Methodology**:
    - For each connectivity, extracts all non-percolating clusters.
    - Counts the occurrence of each cluster size $s$ to obtain $n_s$.
    - Stores raw counts for later averaging across all frames.
  - **Percolation Theory Significance**:
    - The cluster size distribution $n_s$ describes the total number of clusters of size $s$ accross each frames.
    - Near the percolation threshold $p_c$, this distribution exhibits power-law scaling: $n_s \sim s^{-\tau}$, where $\tau$ is the Fisher exponent (approximately 2.18 for 3D systems).
    - Percolating clusters are excluded to focus on finite cluster statistics, as the infinite percolating cluster would dominate near and above the threshold.

- `finalize() -> Dict`
  - Computes final statistics by averaging cluster counts over all processed frames.
  - **Calculation Details**:
    - For each cluster size, calculates the total number of clusters in each frame.
    - Computes standard deviation accounting for frames where specific sizes may not appear (zero counts included).
    - Uses Bessel's correction (ddof=1) for unbiased standard deviation estimation.
  - Idempotent operation: can be called multiple times safely.
  - Returns dictionary containing concentrations, size distributions, and standard deviations.

- `get_result() -> Dict[str, Dict]`
  - Returns the finalized analysis results without recalculation.
  - Provides access to concentrations, cluster size distributions, and standard deviations for all connectivity types.

- `print_to_file() -> None`
  - Writes finalized results to separate CSV-formatted data files for each connectivity type.
  - Output files: `cluster_size_distribution-{connectivity}.dat` in the export directory.
  - Includes header with metadata: connectivity type, date, number of frames analyzed, and column descriptions.
  - Format: `Connectivity_type,Concentration,Cluster_size,N_clusters_per_frame,Standard_deviation`
  - Data sorted by cluster size in descending order for convenient plotting.
  - Automatically removes duplicate lines to prevent data redundancy.

- `_write_header(connectivity: str) -> None`
  - Internal method to write the file header with analysis metadata for a specific connectivity type.
  - Respects overwrite settings: creates new file or appends to existing based on configuration.
  - Includes timestamp and frame count for reproducibility.

### ConcentrationAnalyzer

Computes the node concentration for each connectivity type. The concentration is defined as the ratio of the number of nodes participating in clusters of a given connectivity to the total number of networking nodes.

#### Inheritance

Inherits from `BaseAnalyzer`.

#### Initialization

```python
ConcentrationAnalyzer(settings: Settings)
```
Creates an analyzer for computing cluster concentrations across trajectory frames.

##### Parameters

- `settings` (`Settings`): Configuration settings specifying analysis parameters and output paths.


#### Attributes

- `_raw_concentrations` (`Dict[str, List[float]]`): Raw per-frame concentration values for each connectivity type.
- `concentrations` (`Dict[str, float]`): Final mean concentration for each connectivity type.
- `std` (`Dict[str, float]`): Standard deviation of concentrations across frames.
- `error` (`Dict[str, float]`): Standard error of concentrations across frames.
- `_finalized` (bool): Flag indicating whether final calculations have been performed.


#### Methods

- `analyze(frame: Frame, connectivities: List[str]) -> None`
  - Analyzes a single frame to compute the concentration for each connectivity type.
  - **Calculation Methodology**:
    - Extracts concentration values from the frame for each connectivity type.
    - Concentration is calculated as the ratio of nodes in clusters to total nodes in the system.
    - Stores raw values for later aggregation across all frames.
  - **Percolation Theory Significance**:
    - Concentration represents the fraction of the system participating in the cluster network.
    - Above the percolation threshold, the concentration of the spanning cluster corresponds to the order parameter $P_\infty$, which characterizes the phase transition.
    - Below threshold, concentration measures the extent of finite cluster formation.

- `finalize() -> Dict[str, Dict[str, float]]`
  - Computes final statistics by averaging over all processed frames.
  - **Calculation Details**:
    - Calculates mean concentration across all frames.
    - Computes standard deviation and error using Bessel's correction (ddof=1).
  - Idempotent operation: can be called multiple times safely.
  - Returns dictionary containing concentrations, standard deviations, and error.

- `get_result() -> Dict[str, Dict[str, float]]`
  - Returns the finalized analysis results without recalculation.
  - Provides access to concentrations, standard deviations, and error for all connectivity types.

- `print_to_file() -> None`
  - Writes finalized results to a CSV-formatted data file.
  - Output file: `concentrations.dat` in the export directory.
  - Includes header with metadata: date, number of frames analyzed, and column descriptions.
  - Format: `Connectivity_type,Concentration,Standard_deviation,Standard_error`
  - Automatically removes duplicate lines to prevent data redundancy.

- `_write_header() -> None`
  - Internal method to write the file header with analysis metadata.
  - Respects overwrite settings: creates new file or appends to existing based on configuration.
  - Includes timestamp and frame count for reproducibility.

### CorrelationLengthAnalyzer

Computes the correlation length (xi) of the cluster size distribution. Defined as ``xi^2 = sum(2 * R_s^2 * s^2 * n_s) / sum(s^2 * n_s)`` where *R_s* is the gyration radius and *s* the cluster size. Only non-percolating clusters contribute.

#### Inheritance

Inherits from `BaseAnalyzer`.

#### Initialization

```python
CorrelationLengthAnalyzer(settings: Settings)
```
Creates an analyzer for computing correlation lengths across trajectory frames.

##### Parameters

- `settings` (`Settings`): Configuration settings specifying analysis parameters and output paths.


#### Attributes

- `_raw_correlation_lengths` (`Dict[str, List[float]]`): Raw per-frame correlation length values for each connectivity type.
- `_raw_concentrations` (`Dict[str, List[float]]`): Raw per-frame concentration values for each connectivity.
- `correlation_length` (`Dict[str, float]`): Final mean correlation length for each connectivity type.
- `std` (`Dict[str, float]`): Standard deviation of correlation lengths across frames.
- `error` (`Dict[str, float]`): Standard error of correlation lengths across frames.
- `concentrations` (`Dict[str, float]`): Mean concentration for each connectivity type.
- `_finalized` (bool): Flag indicating whether final calculations have been performed.


#### Methods

- `analyze(frame: Frame, connectivities: List[str]) -> None`
  - Analyzes a single frame to compute the correlation length for each connectivity type.
  - **Calculation Methodology**:
    - For each connectivity, extracts all non-percolating clusters.
    - Computes the correlation length using the weighted second moment of gyration radii:
      $$
      \xi^2 = \frac{\sum_s 2 R_{g,s}^2 s^2 n_s}{\sum_s s^2 n_s}
      $$
      where $R_{g,s}$ is the radius of gyration for clusters of size $s$, and $n_s$ is the number of such clusters.
    - Stores raw values for later aggregation across all frames.
  - **Percolation Theory Significance**:
    - The correlation length defines the characteristic size scale over which structural correlations persist in the system.
    - Near the percolation threshold $p_c$, $\xi$ diverges as a power law: $\xi \sim |p - p_c|^{-\nu}$, where $\nu$ is the correlation length critical exponent (approximately 0.88 for standard percolation in 3D systems).
    - Below threshold, $\xi$ represents the typical size of finite clusters and sets the cutoff for power-law scaling in the cluster size distribution.
    - At the percolation threshold, the divergence of $\xi$ reflects the emergence of correlations at all length scales, a hallmark of critical phenomena.
    - Percolating clusters are excluded as they represent infinite clusters whose spatial extent would dominate the calculation.

- `finalize() -> Dict[str, Dict[str, float]]`
  - Computes final statistics by averaging over all processed frames.
  - **Calculation Details**:
    - Calculates mean correlation length across all frames.
    - Computes standard deviation and error using Bessel's correction (ddof=1).
  - Idempotent operation: can be called multiple times safely.
  - Returns dictionary containing concentrations, correlation lengths, standard deviations, and errors.

- `get_result() -> Dict[str, Dict[str, float]]`
  - Returns the finalized analysis results without recalculation.
  - Provides access to concentrations, correlation lengths, standard deviations, and errors for all connectivity types.

- `print_to_file() -> None`
  - Writes finalized results to a CSV-formatted data file.
  - Output file: `correlation_length.dat` in the export directory.
  - Includes header with metadata: date, number of frames analyzed, and column descriptions.
  - Format: `Connectivity_type,Concentration,Correlation_length,Standard_deviation,Standard_error`
  - Automatically removes duplicate lines to prevent data redundancy.

- `_write_header() -> None`
  - Internal method to write the file header with analysis metadata.
  - Respects overwrite settings: creates new file or appends to existing based on configuration.
  - Includes timestamp and frame count for reproducibility.

### GyrationRadiusAnalyzer

Computes the mean gyration radius binned by cluster size for each connectivity. Collects the gyration radius of all non-percolating clusters, groups them by cluster size, and averages over all processed frames. Results are written to one file per connectivity.

#### Inheritance
Inherits from `BaseAnalyzer`.

#### Initialization

```python
GyrationRadiusAnalyzer(settings: Settings)
```
Creates an analyzer for computing gyration radius distributions across trajectory frames.

##### Parameters

- `settings` (`Settings`): Configuration settings specifying analysis parameters and output paths.


#### Attributes

- `_raw_gyration_radii` (`Dict[str, Dict[int, List[float]]]`): Raw per-frame gyration radii for each cluster size and connectivity type.
- `_raw_concentrations` (`Dict[str, List[float]]`): Raw per-frame concentration values for each connectivity.
- `gyration_radii` (`Dict[str, Dict[int, float]]`): Final averaged gyration radius for each cluster size and connectivity type.
- `std` (`Dict[str, Dict[int, float]]`): Standard deviation of gyration radii for each size across frames.
- `concentrations` (`Dict[str, float]`): Mean concentration for each connectivity type.
- `_finalized` (bool): Flag indicating whether final calculations have been performed.


#### Methods

- `analyze(frame: Frame, connectivities: List[str]) -> None`
  - Collects gyration radii of non-percolating clusters for each connectivity type, grouped by cluster size.
  - **Calculation Methodology**:
    - For each connectivity, extracts all non-percolating clusters.
    - Records the gyration radius $R_g$ for each cluster, binned by its size $s$.
    - The gyration radius is computed as:
      $$
      R_g^2 = \frac{1}{N} \sum_{i=1}^{N} |\mathbf{r}_i - \mathbf{r}_{\text{COM}}|^2
      $$
      where $N$ is the number of nodes in the cluster, $\mathbf{r}_i$ are node positions, and $\mathbf{r}_{\text{COM}}$ is the center of mass.
    - Stores raw values for later averaging across all frames.
  - **Percolation Theory Significance**:
    - The gyration radius quantifies the spatial extent of clusters and their degree of compactness.
    - Near the percolation threshold, clusters exhibit fractal geometry with scaling relation: $R_g \sim s^{1/d_f}$, where $d_f$ is the fractal dimension (approximately 2.53 for standard percolation in 3D systems).
    - For compact clusters, $d_f = d$ (the spatial dimension), while fractal clusters have $d_f < d$.
    - The scaling behavior of $R_g$ versus $s$ reveals whether clusters are ramified (fractal) or compact, providing insight into the mechanisms of cluster growth.
    - At criticality, clusters are self-similar fractals with characteristic power-law scaling that persists across multiple length scales.
    - Percolating clusters are excluded to focus on finite cluster geometry.

- `finalize() -> Dict[str, Dict]`
  - Computes average gyration radius and standard deviation for each cluster size across all processed frames.
  - **Calculation Details**:
    - For each cluster size, calculates the mean gyration radius across all observations.
    - Computes standard deviation using Bessel's correction (ddof=1) when multiple observations exist.
  - Idempotent operation: can be called multiple times safely.
  - Returns dictionary containing concentrations, gyration radii, and standard deviations.

- `get_result() -> Dict[str, Dict]`
  - Returns the finalized analysis results without recalculation.
  - Provides access to concentrations, gyration radii distributions, and standard deviations for all connectivity types.

- `print_to_file() -> None`
  - Writes finalized results to separate CSV-formatted data files for each connectivity type.
  - Output files: `gyration_radius_distribution-{connectivity}.dat` in the export directory.
  - Includes header with metadata: date, number of frames analyzed, and column descriptions.
  - Format: `Connectivity_type,Concentration,Cluster_size,Gyration_radius,Standard_deviation`
  - Data sorted by cluster size in descending order for convenient plotting.
  - Automatically removes duplicate lines to prevent data redundancy.

- `_write_header() -> None`
  - Internal method to write the file header with analysis metadata for each connectivity type.
  - Respects overwrite settings: creates new file or appends to existing based on configuration.
  - Includes timestamp and frame count for reproducibility.

### LargestClusterSizeAnalyzer

Computes the size of the largest cluster for each connectivity type. Tracks the maximum cluster size per connectivity across frames and provides ensemble-averaged results with standard deviation and error.

#### Inheritance

Inherits from `BaseAnalyzer`.

#### Initialization

```python
LargestClusterSizeAnalyzer(settings: Settings)
```
Creates an analyzer for computing largest cluster sizes across trajectory frames.

##### Parameters

- `settings` (`Settings`): Configuration settings specifying analysis parameters and output paths.


#### Attributes

- `_raw_sizes` (`Dict[str, List[float]]`): Raw per-frame largest cluster sizes for each connectivity type.
- `_raw_concentrations` (`Dict[str, List[float]]`): Raw per-frame concentration values for each connectivity.
- `largest_cluster_sizes` (`Dict[str, float]`): Final mean largest cluster size for each connectivity type.
- `std` (`Dict[str, float]`): Standard deviation of largest cluster sizes across frames.
- `error` (`Dict[str, float]`): Standard error of largest cluster sizes across frames.
- `concentrations` (`Dict[str, float]`): Mean concentration for each connectivity type.
- `_finalized` (bool): Flag indicating whether final calculations have been performed.


#### Methods

- `analyze(frame: Frame, connectivities: List[str]) -> None`
  - Analyzes a single frame to identify the largest cluster size for each connectivity type.
  - **Calculation Methodology**:
    - For each connectivity, examines all clusters (including percolating clusters).
    - Identifies the maximum cluster size: $S_{\text{max}} = \max\{s_1, s_2, \ldots, s_n\}$.
    - Stores raw values for later aggregation across all frames.
  - **Percolation Theory Significance**:
    - The largest cluster size is intimately connected to the percolation order parameter $P_\infty$.
    - Below the percolation threshold ($p < p_c$), the largest cluster is finite and grows slowly with system size.
    - At the threshold ($p = p_c$), the largest cluster exhibits critical scaling: $S_{\text{max}} \sim L^{d\beta/\nu}$, where $L$ is the system size, $d$ is dimensionality, and $\beta, \nu$ are critical exponents.
    - Above threshold ($p > p_c$), the largest cluster becomes the percolating cluster, scaling with system volume.
    - The largest cluster captures the dominant connectivity pathway and is the primary contributor to transport properties in percolating systems.
    - Unlike other metrics, this analyzer includes percolating clusters, as they represent the dominant structure above threshold.

- `finalize() -> Dict[str, Dict[str, float]]`
  - Computes final statistics by averaging over all processed frames.
  - **Calculation Details**:
    - Calculates mean largest cluster size across all frames.
    - Computes standard deviation and error using Bessel's correction (ddof=1).
  - Idempotent operation: can be called multiple times safely.
  - Returns dictionary containing concentrations, largest cluster sizes, standard deviations, and errors.

- `get_result() -> Dict[str, Dict[str, float]]`
  - Returns the finalized analysis results without recalculation.
  - Provides access to concentrations, largest cluster sizes, standard deviations, and errors for all connectivity types.

- `print_to_file() -> None`
  - Writes finalized results to a CSV-formatted data file.
  - Output file: `largest_cluster_size.dat` in the export directory.
  - Includes header with metadata: date, number of frames analyzed, and column descriptions.
  - Format: `Connectivity_type,Concentration,Largest_cluster_size,Standard_deviation,Standard_error`
  - Automatically removes duplicate lines to prevent data redundancy.

- `_write_header() -> None`
  - Internal method to write the file header with analysis metadata.
  - Respects overwrite settings: creates new file or appends to existing based on configuration.
  - Includes timestamp and frame count for reproducibility.

### OrderParameterAnalyzer

Computes the percolation order parameter (P_inf) for each connectivity type. The order parameter is the fraction of networking nodes that belong to the percolating cluster.

#### Inheritance

Inherits from `BaseAnalyzer`.

#### Initialization

```python
OrderParameterAnalyzer(settings: Settings)
```
Creates an analyzer for computing order parameters across trajectory frames.

##### Parameters

- `settings` (`Settings`): Configuration settings specifying analysis parameters and output paths.


#### Attributes

- `_raw_order_parameters` (`Dict[str, List[float]]`): Raw per-frame order parameter values for each connectivity type.
- `_raw_concentrations` (`Dict[str, List[float]]`): Raw per-frame concentration values for each connectivity.
- `order_parameters` (`Dict[str, float]`): Final mean order parameter for each connectivity type.
- `std` (`Dict[str, float]`): Standard deviation of order parameters across frames.
- `concentrations` (`Dict[str, float]`): Mean concentration for each connectivity type.
- `_finalized` (bool): Flag indicating whether final calculations have been performed.


#### Methods

- `analyze(frame: Frame, connectivities: List[str]) -> None`
  - Analyzes a single frame to compute the order parameter for each connectivity type.
  - **Calculation Methodology**:
    - For each connectivity, identifies all percolating clusters.
    - Selects the largest percolating cluster if multiple exist.
    - Extracts the first component of the order parameter vector (directional order parameter).
    - If no percolating cluster exists, assigns order parameter of zero.
    - Stores raw values for later aggregation across all frames.
  - **Percolation Theory Significance**:
    - The order parameter $P_\infty$ represents the probability that a randomly selected node belongs to the infinite percolating cluster.
    - Below the percolation threshold ($p < p_c$), no infinite cluster exists and $P_\infty = 0$.
    - Above threshold ($p > p_c$), the order parameter emerges continuously: $P_\infty \sim (p - p_c)^\beta$, where $\beta$ is the order parameter critical exponent (approximately 0.41 for standard percolation in 3D systems).
    - This metric is analogous to spontaneous magnetization in the Ising model and serves as the primary signature of the percolation phase transition.
    - The order parameter distinguishes the non-percolating phase (isolated finite clusters) from the percolating phase (system-spanning connectivity).
    - The directional component reflects percolation along a specific spatial direction, relevant for anisotropic systems or directional transport properties.

- `finalize() -> Dict[str, Dict[str, float]]`
  - Computes final statistics by averaging over all processed frames.
  - **Calculation Details**:
    - Calculates mean order parameter across all frames.
    - Computes standard deviation using Bessel's correction (ddof=1).
  - Idempotent operation: can be called multiple times safely.
  - Returns dictionary containing concentrations, order parameters, and standard deviations.

- `get_result() -> Dict[str, Dict[str, float]]`
  - Returns the finalized analysis results without recalculation.
  - Provides access to concentrations, order parameters, and standard deviations for all connectivity types.

- `print_to_file() -> None`
  - Writes finalized results to a CSV-formatted data file.
  - Output file: `order_parameter.dat` in the export directory.
  - Includes header with metadata: date, number of frames analyzed, and column descriptions.
  - Format: `Connectivity_type,Concentration,Order_parameter,Standard_deviation`
  - Automatically removes duplicate lines to prevent data redundancy.

- `_write_header() -> None`
  - Internal method to write the file header with analysis metadata.
  - Respects overwrite settings: creates new file or appends to existing based on configuration.
  - Includes timestamp and frame count for reproducibility.

### PercolationProbabilityAnalyzer

Computes the percolation probability (Pi) for each connectivity type. The percolation probability is the fraction of frames in which at least one cluster spans the simulation box *in all three dimensions*. A value of 1.0 indicates percolation in every frame; 0.0 indicates it never occurs.

#### Inheritance

Inherits from `BaseAnalyzer`.

#### Initialization

```python
PercolationProbabilityAnalyzer(settings: Settings)
```
Creates an analyzer for computing percolation probabilities across trajectory frames.

##### Parameters

- `settings` (`Settings`): Configuration settings specifying analysis parameters and output paths.


#### Attributes

- `_raw_percolation_prob_x` (`Dict[str, List[float]]`): Raw per-frame binary indicators (1.0 or 0.0) of percolation occurrence for each connectivity type.
- `_raw_concentrations` (`Dict[str, List[float]]`): Raw per-frame concentration values for each connectivity.
- `percolation_probabilities` (`Dict[str, float]`): Final mean percolation probability for each connectivity type.
- `std` (`Dict[str, float]`): Standard deviation of percolation probabilities across frames.
- `concentrations` (`Dict[str, float]`): Mea concentration for each connectivity type.
- `_finalized` (bool): Flag indicating whether final calculations have been performed.


#### Methods

- `analyze(frame: Frame, connectivities: List[str]) -> None`
  - Check whether percolation occurs for each connectivity in the frame.
  - **Calculation Methodology**:
    - For each connectivity, examines all clusters to detect percolation.
    - A cluster is considered percolating if it spans the simulation box in all three dimensions (percolation_probability contains "xyz").
    - Records binary outcome: 1.0 if percolation occurs, 0.0 otherwise.
    - Stores raw values for later averaging across all frames.
  - **Percolation Theory Significance**:
    - The percolation probability $\Pi(p)$ represents the ensemble probability of observing system-spanning connectivity at a given occupation probability $p$.
    - At the percolation threshold $p_c$, $\Pi$ undergoes a sharp transition from 0 to 1 in the thermodynamic limit.
    - For finite systems, the transition is smoothed over a narrow range around $p_c$, with width scaling as $L^{-1/\nu}$ where $L$ is system size and $\nu$ is the correlation length exponent.
    - The inflection point of $\Pi(p)$ provides an accurate estimator of the percolation threshold.
    - Unlike the order parameter $P_\infty$ (which measures the fraction of nodes in the percolating cluster), $\Pi$ is a Boolean property indicating the mere existence of percolation.
    - Directional percolation probability is particularly relevant for anisotropic systems and transport phenomena where connectivity along specific axes determines material properties.

- `finalize() -> Dict[str, Dict[str, float]]`
  - Computes final statistics by averaging over all processed frames.
  - **Calculation Details**:
    - Calculates mean percolation probability: $\Pi = \frac{1}{N_{\text{frames}}} \sum_{i=1}^{N_{\text{frames}}} \delta_i$, where $\delta_i = 1$ if frame $i$ percolates and 0 otherwise.
    - Computes standard deviation using Bessel's correction (ddof=1).
  - Idempotent operation: can be called multiple times safely.
  - Returns dictionary containing concentrations, percolation probabilities, and standard deviations.

- `get_result() -> Dict[str, Dict[str, float]]`
  - Returns the finalized analysis results without recalculation.
  - Provides access to concentrations, percolation probabilities, and standard deviations for all connectivity types.

- `print_to_file() -> None`
  - Writes finalized results to a CSV-formatted data file.
  - Output file: `percolation_probability.dat` in the export directory.
  - Includes header with metadata: date, number of frames analyzed, and column descriptions.
  - Format: `Connectivity_type,Concentration,Percolation_probability,Standard_deviation`
  - Automatically removes duplicate lines to prevent data redundancy.

- `_write_header() -> None`
  - Internal method to write the file header with analysis metadata.
  - Respects overwrite settings: creates new file or appends to existing based on configuration.
  - Includes timestamp and frame count for reproducibility.
  - Analysis checks percolation in all three dimensions.

### SpanningClusterSizeAnalyzer

Computes the size of the largest finite (non-percolating) cluster. Tracks the maximum non-percolating cluster size per connectivity across frames.

#### Inheritance

Inherits from `BaseAnalyzer`.

#### Initialization

```python
SpanningClusterSizeAnalyzer(settings: Settings)
```
Creates an analyzer for computing spanning cluster sizes across trajectory frames.

##### Parameters

- `settings` (`Settings`): Configuration settings specifying analysis parameters and output paths.


#### Attributes

- `_raw_spanning_sizes` (`Dict[str, List[float]]`): Raw per-frame largest non-percolating cluster sizes for each connectivity type.
- `_raw_concentrations` (`Dict[str, List[float]]`): Raw per-frame concentration values for each connectivity.
- `spanning_cluster_sizes` (`Dict[str, float]`): Final mean spanning cluster size for each connectivity type.
- `std` (`Dict[str, float]`): Standard deviation of spanning cluster sizes across frames.
- `error` (`Dict[str, float]`): Standard error of spanning cluster sizes across frames.
- `concentrations` (`Dict[str, float]`): Mean concentration for each connectivity type.
- `_finalized` (bool): Flag indicating whether final calculations have been performed.


#### Methods

- `analyze(frame: Frame, connectivities: List[str]) -> None`
  - Analyzes a single frame to identify the largest non-percolating cluster size for each connectivity type.
  - **Calculation Methodology**:
    - For each connectivity, extracts all non-percolating (finite) clusters.
    - Identifies the maximum cluster size among finite clusters: $S_{\text{span}} = \max\{s_1, s_2, \ldots, s_n\}$ where clusters do not span the system.
    - If no finite clusters exist (all clusters percolate), records size as zero.
    - Stores raw values for later aggregation across all frames.
  - **Percolation Theory Significance**:
    - In the sub-critical regime ($p < p_c$), the spanning cluster represents the largest structural feature and governs finite-size scaling behavior.
    - As the system approaches the percolation threshold from below, $S_{\text{span}}$ grows rapidly, following a power law: $S_{\text{span}} \sim |p_c - p|^{-1/\sigma}$, where $\sigma$ is related to critical exponents and the fractal dimension in percolation theory.
    - This metric is complementary to the largest cluster size analyzer—while that includes percolating clusters, the spanning cluster size focuses exclusively on the finite cluster distribution.
    - Above the threshold, $S_{\text{span}}$ measures the size of the largest finite cluster coexisting with the percolating network, providing insight into the residual finite cluster population.
    - The spanning cluster often dominates transport and structural properties in the sub-critical regime before system-wide percolation emerges.

- `finalize() -> Dict[str, Dict[str, float]]`
  - Computes final statistics by averaging over all processed frames.
  - **Calculation Details**:
    - Calculates mean spanning cluster size across all frames.
    - Computes standard deviation and error using Bessel's correction (ddof=1).
  - Idempotent operation: can be called multiple times safely.
  - Returns dictionary containing concentrations, spanning cluster sizes, standard deviations, and errors.

- `get_result() -> Dict[str, Dict[str, float]]`
  - Returns the finalized analysis results without recalculation.
  - Provides access to concentrations, spanning cluster sizes, standard deviations, and errors for all connectivity types.

- `print_to_file() -> None`
  - Writes finalized results to a CSV-formatted data file.
  - Output file: `spanning_cluster_size.dat` in the export directory.
  - Includes header with metadata: date, number of frames analyzed, and column descriptions.
  - Format: `Connectivity_type,Concentration,Spanning_cluster_size,Standard_deviation,Standard_error`
  - Automatically removes duplicate lines to prevent data redundancy.

- `_write_header() -> None`
  - Internal method to write the file header with analysis metadata.
  - Respects overwrite settings: creates new file or appends to existing based on configuration.
  - Includes timestamp and frame count for reproducibility.


## Parser module (`io` module)

Utility class for discovering trajectory files in a directory and parsing associated metadata from an `info.csv` file. Automatically sorts files and validates metadata consistency.

***

#### Initialization

```python
Parser(file_location: str, format: str)
```
Creates a parser instance that scans for trajectory files and loads metadata.

##### Parameters

- `file_location` (str): Path to either a trajectory file or directory containing trajectory files.
- `format` (str): File extension format to search for (e.g., `"xyz"`, `"dat"`).

##### Raises

- `ValueError`: If the specified `file_location` does not exist.
- `ValueError`: If `info.csv` is empty or has mismatched line count with trajectory files.
- `ValueError`: If `info.csv` is not found in the directory.

***

#### Attributes

- `file_location` (str): Path to the file or directory being parsed.
- `format` (str): File extension format to filter.
- `files` (List[str]): List of discovered trajectory file paths, sorted alphabetically.
- `infos` (Dict[str, List[float]]): Dictionary mapping column names from `info.csv` to their values.

***

#### Methods

- `parse() -> List[str]`
  - Discovers and collects all trajectory files matching the specified format.
  - **Behavior**:
    - If `file_location` is a file: Scans the parent directory for matching files.
    - If `file_location` is a directory: Scans that directory for matching files.
    - Sorts files alphabetically for consistent ordering.
  - **Returns**: List of full file paths (also stored in `self.files`).

- `parse_infos() -> None`
  - Reads and parses the `info.csv` metadata file from the directory.
  - **File Format Requirements**:
    - Must be named `info.csv` and located in `file_location` directory.
    - First line: comma-separated column names (header).
    - Subsequent lines: comma-separated values, one per trajectory file.
    - Must contain a `project_name` column (stored as strings).
    - All other columns stored as floats.
    - Number of data lines must match number of trajectory files.
  - **Returns**: None (populates `self.infos` dictionary).
  - **Use Case**: Associates metadata like probability, density, temperature, or pressure with each trajectory file for plotting and analysis.
  - **Raises**:
    - `ValueError`: If `info.csv` is empty.
    - `ValueError`: If line count doesn't match number of trajectory files.
    - `ValueError`: If `info.csv` is not found.

- `get_files() -> List[str]`
  - Returns the list of discovered trajectory files.
  - **Returns**: List of full file paths.
  - If `files` list is empty, automatically calls `parse()` first.

- `get_infos() -> Dict[str, List[float]]`
  - Returns the parsed metadata dictionary.
  - **Returns**: Dictionary mapping column names to lists of values.
  - If `infos` dictionary is empty, automatically calls `parse_infos()` first.

***

#### Usage Example

```python
from nexus.io.parser import Parser

# Parse trajectory files and metadata
parser = Parser(
    file_location="/path/to/trajectories/",
    format="xyz"
)

# Get sorted list of trajectory files
trajectory_files = parser.get_files()
print(f"Found {len(trajectory_files)} trajectory files")

# Get associated metadata
metadata = parser.get_infos()
print(f"Metadata columns: {list(metadata.keys())}")

# Access specific metadata
temperatures = metadata.get("temperature", [])
project_names = metadata.get("project_name", [])

# Iterate through files with metadata
for i, (file, temp, name) in enumerate(zip(trajectory_files, temperatures, project_names)):
    print(f"{i}: {name} at {temp}K - {file}")
```

#### info.csv Format Example

```csv
project_name,temperature,pressure,density
system_001,300.0,1.0,0.85
system_002,350.0,1.0,0.82
system_003,400.0,1.0,0.79
```

## Writer module (`io` module)

### WriterFactory

Factory class for creating and managing file writer instances based on writer type. Implements the factory pattern to provide centralized writer instantiation with shared settings.

#### Initialization

```python
WriterFactory(settings: Settings)
```
Creates a writer factory instance and registers all available writer types.

##### Parameters

- `settings` (`Settings`): Configuration settings shared across all writer instances.


#### Attributes

- `_writers` (dict): Internal registry mapping writer class names to writer classes.
- `_settings` (`Settings`): Configuration settings used to initialize writers.


#### Registered Writers

The factory automatically registers the following writer types:
- **ClustersWriter**: Writes cluster configuration and unwrapped position data.
- **LogsWriter**: Writes execution logs and run information.
- **PerformanceWriter**: Writes performance metrics and timing data.
- **MultipleFilesSummaryWriter**: Writes summary statistics across multiple trajectory files.


#### Methods

- `register_writer(writer: BaseWriter) -> None`
  - Registers a new writer class in the factory.
  - **Parameters**:
    - `writer` (BaseWriter): Writer class (not instance) to register.
  - **Returns**: None
  - Stores the writer class using its class name as the key.

- `get_writer(name: str, mode: str = "all") -> Optional[BaseWriter]`
  - Retrieves an instantiated writer of the specified type.
  - **Parameters**:
    - `name` (str): Name of the writer class to instantiate. Options:
      - `"ClustersWriter"`: For cluster data output.
      - `"LogsWriter"`: For execution logs.
      - `"PerformanceWriter"`: For performance metrics.
      - `"MultipleFilesSummaryWriter"`: For multi-file summaries.
    - `mode` (str): Output mode for `MultipleFilesSummaryWriter`. Default: `"all"`.
  - **Returns**: Instantiated writer object or `None` if writer name not recognized.
  - Each call creates a new writer instance with the factory's settings.

### ClustersWriter

Writes cluster data to XYZ and bond files, including unwrapped atomic positions for both networking nodes and decorating atoms (e.g., bridging oxygens). Supports multiple output modes for flexible data organization.

#### Inheritance

Inherits from `BaseWriter`.

#### Initialization

```python
ClustersWriter(settings: Settings)
```
Creates a cluster file writer with specified configuration settings.

##### Parameters

- `settings` (`Settings`): Configuration settings controlling output behavior and file paths.

***

#### Attributes

- `_settings` (`Settings`): Configuration settings.
- `_clusters` (`List[Cluster]`): List of clusters to write, sorted by size (largest first).

***

#### Methods

- `set_clusters(clusters: List[Cluster]) -> None`
  - Sets the clusters to be written and sorts them by size in descending order.
  - **Parameters**:
    - `clusters` (List[Cluster]): List of cluster objects to write.
  - **Returns**: None

- `write() -> None`
  - Writes cluster data according to the configured print mode.
  - **Print Modes** (from `settings.clustering.print_mode`):
    - `"none"`: No output written.
    - `"all"`: All clusters combined into single XYZ and bonds files.
    - `"connectivity"`: Clusters grouped by connectivity type into separate files.
    - `"individual"`: Each cluster written to its own file.
  - **Returns**: None

- `_write_all() -> None`
  - Writes all clusters to a single combined XYZ and bonds file.
  - **Output Files**:
    - `unwrapped_clusters/all_unwrapped_clusters-frame_{frame_id}.xyz`
    - `unwrapped_clusters/all_unwrapped_clusters-frame_{frame_id}.bonds`
  - **Behavior**:
    - Counts total unique atoms (networking + decorating nodes).
    - Builds global node ID to local file index mapping.
    - Writes all cluster atoms sequentially.
    - Writes bond information using local indices.

- `_write_connectivity() -> None`
  - Writes clusters grouped by connectivity type to separate files per connectivity.
  - **Output Files**:
    - `unwrapped_clusters/{connectivity}/{connectivity}_unwrapped_clusters-frame_{frame_id}.xyz`
    - `unwrapped_clusters/{connectivity}/{connectivity}_unwrapped_clusters-frame_{frame_id}.bonds`
  - **Behavior**:
    - Groups clusters by connectivity descriptor.
    - Marks first (largest) cluster in each group as spanning.
    - Creates subdirectories for each connectivity type.

- `_write_individual() -> None`
  - Writes each cluster to its own separate XYZ and bonds files.
  - **Output Files**:
    - `unwrapped_clusters/{connectivity}/cluster-frame_{frame_id}-id_{root_id}.xyz`
    - `unwrapped_clusters/{connectivity}/cluster-frame_{frame_id}-id_{root_id}.bonds`
  - **Behavior**:
    - Creates separate files for every cluster.
    - Organizes files by connectivity type in subdirectories.

- `_write_header_comment(f: TextIO, cluster: Cluster) -> None`
  - Writes the extended XYZ format header comment line with lattice and properties.
  - **Parameters**:
    - `f` (TextIO): Open file handle.
    - `cluster` (Cluster): Cluster object containing lattice information.
  - **Header Format**:
    ```
    Lattice="{lxx} {lxy} {lxz} {lyx} {lyy} {lyz} {lzx} {lzy} {lzz}" Properties=species:S:1:index:I:1:pos:R:3:cluster_id:I:1:coordination:I:1:percolating:I:1:spanning:I:1
    ```
  - **Properties**:
    - `species`: Element symbol (string)
    - `index`: Global node ID (integer)
    - `pos`: Unwrapped position (3 real numbers)
    - `cluster_id`: Root cluster ID (integer)
    - `coordination`: Coordination number (integer)
    - `percolating`: Percolation dimensionality (integer)
    - `spanning`: Spanning flag (integer, 1 or 0)

- `_write_cluster_atoms(f: TextIO, cluster: Cluster, start_index: int, unique_ids: List|None = None) -> Dict[int, int]`
  - Writes all atoms (networking and decorating nodes) for a cluster to file.
  - **Parameters**:
    - `f` (TextIO): Open file handle.
    - `cluster` (Cluster): Cluster object to write.
    - `start_index` (int): Starting local index for this cluster's atoms.
    - `unique_ids` (List|None): Set of already-written node IDs to prevent duplicates. Default: `None`.
  - **Returns**: Dictionary mapping global node IDs to local file indices, and updated unique IDs set.
  - **Atom Line Format**:
    ```
    {symbol} {global_id} {x} {y} {z} {cluster_id} {coordination} {percolation_dims} {spanning}
    ```
  - **Behavior**:
    - Writes primary networking nodes first.
    - Writes decorating atoms (e.g., bridging oxygens) second.
    - Tracks unique IDs to avoid writing shared decorating atoms multiple times.
    - Builds mapping from global node IDs to sequential local indices.

- `_write_cluster_bonds(f: TextIO, cluster: Cluster, id_map: Dict[int, int]) -> None`
  - Writes bond connectivity information using local file indices.
  - **Parameters**:
    - `f` (TextIO): Open file handle for bonds file.
    - `cluster` (Cluster): Cluster object containing linkage information.
    - `id_map` (Dict[int, int]): Mapping from global node IDs to local file indices.
  - **Returns**: None
  - **Bond Line Format**:
    ```
    {symbol1}({local_index1})-{symbol2}({local_index2})
    ```
  - **Behavior**:
    - Skips bond writing if clustering criterion is `"distance"`.
    - Only writes bonds between networking nodes (uses cluster linkages).
    - Translates global node IDs to local file indices for visualization compatibility.

***

#### Output File Formats

**XYZ Format**:
- Line 1: Total atom count
- Line 2: Extended XYZ header with lattice and properties
- Remaining lines: Atom data (symbol, ID, position, properties)

**Bonds Format**:
- One bond per line in format: `Symbol1(index1)-Symbol2(index2)`
- Uses local file indices (1-based) for visualization software compatibility

***

#### Usage Example

```python
config_clustering = c.ClusteringSettings(
  ...
  with_printed_unwrapped_clusters=True,
  print_mode="connectivity",  # "all", "connectivity", "individual"
  ...
)
```

### LogsWriter

Writes execution logs containing version information and configuration settings to a text file. Provides a record of analysis runs for reproducibility and debugging.

#### Inheritance

Inherits from `BaseWriter`.

#### Initialization

```python
LogsWriter(settings: Settings)
```
Creates a logs file writer with specified configuration settings.

##### Parameters

- `settings` (`Settings`): Configuration settings to be logged.

***

#### Attributes

- `_settings` (`Settings`): Configuration settings to write to log file.

***

#### Methods

- `write() -> None`
  - Writes the log file containing version banner and complete settings configuration.
  - **Output File**: `{export_directory}/log.txt`
  - **File Contents**:
    - ASCII art title banner (from `print_title_to_file`)
    - Package version number
    - Complete string representation of all configuration settings
  - **Returns**: None
  - **Use Case**: Creates a permanent record of analysis parameters for each run, facilitating result reproducibility and parameter tracking.

***

#### Usage Example

```python
# General settings
config_general = c.GeneralSettings(
  ...
  save_logs=True,  # Save logs    (save logs to export_directory/logs.txt)
  ...
)
```

### MultipleFilesSummaryWriter

Aggregates analysis results from multiple data files across directory structures and generates summary files. Supports combining data across different trajectory runs or parameter sweeps for comparative analysis.
** Note that this writer is not stable and may not work as expected. **

#### Inheritance

Inherits from `BaseWriter`.

#### Initialization

```python
MultipleFilesSummaryWriter(settings: Settings, mode: str = "all")
```
Creates a summary writer for aggregating multiple analysis result files.

##### Parameters

- `settings` (`Settings`): Configuration settings specifying export directory and paths.
- `mode` (str): Summary output mode. Options:
  - `"all"`: Combines all connectivity types into single summary files.
  - `"connectivity"`: Creates separate summary files for each connectivity type.
  - Default: `"all"`.


#### Attributes

- `_settings` (`Settings`): Configuration settings.
- `_mode` (str): Output mode controlling summary file organization.

#### Methods

- `write() -> None`
  - Scans the export directory for analysis result files and generates summary files.
  - **Behavior**:
    - Recursively walks through `export_directory` to find `.dat` files.
    - Categorizes files by metric type based on filename patterns.
    - Excludes files containing "unwrapped_clusters" or "summary" in their names.
    - Calls `_get_results()` to parse each file category.
    - Calls `_write_summary()` to generate summary files.
  - **Supported Metrics**:
    - Average cluster size
    - Correlation length
    - Largest cluster size
    - Order parameter
    - Percolation probability
    - Spanning cluster size
  - **Output Files**: Summary files named `{metric}_summary.dat` or `{connectivity}_{metric}_summary.dat` depending on mode.
  - **Returns**: None

- `_get_results(files: List[str]) -> Dict[str, List[Tuple[float, float, float]]]`
  - Parses multiple data files and extracts results for each connectivity type.
  - **Parameters**:
    - `files` (List[str]): List of file paths to parse.
  - **Returns**: Dictionary mapping connectivity types to lists of result tuples.
  - **Result Tuple Format**: `(concentration, metric_value, standard_deviation)`
  - **Behavior**:
    - Reads each file line by line.
    - Skips comment lines starting with `#`.
    - Parses CSV format: `connectivity_type,concentration,metric_value,std_dev,...`
    - Groups results by connectivity type.
    - Handles parsing errors gracefully with error messages.

- `_write_summary(results: Dict[str, List[Tuple[float, float, float]]], n_data: int, file_name: str, mode: str = "all") -> None`
  - Writes aggregated summary data to file(s) in the specified format.
  - **Parameters**:
    - `results` (Dict): Parsed results dictionary from `_get_results()`.
    - `n_data` (int): Number of data points (files) aggregated.
    - `file_name` (str): Base filename for the summary file.
    - `mode` (str): Output mode (`"all"` or `"connectivity"`).
  - **Returns**: None
  - **Mode Behaviors**:
    - **"all" mode**:
      - Creates single file with all connectivity types.
      - First columns: metric values for each connectivity.
      - Subsequent columns: standard deviations for each connectivity.
      - Header lists connectivity type indices.
      - Format: Space-separated values, one row per data point.
    - **"connectivity" mode**:
      - Creates separate file for each connectivity type.
      - Filename: `{connectivity}_{metric}_summary.dat`
      - Format: Three columns (concentration, metric_value, std_dev).
      - Includes header identifying connectivity type.

***

#### Output File Formats

**"all" Mode**:
```
# 1 : connectivity_type_1
# 2 : connectivity_type_2
# 3 : std connectivity_type_1
# 4 : std connectivity_type_2
value1_conn1 value1_conn2 std1_conn1 std1_conn2
value2_conn1 value2_conn2 std2_conn1 std2_conn2
...
```

**"connectivity" Mode**:
```
# Connectivity_type : Si-O-Si
# Concentration Average_cluster_size Standard_deviation
0.25 15.3 2.1
0.30 18.7 2.5
0.35 22.1 3.2
...
```

#### Usage Example

```python
from io.writer import MultipleFilesSummaryWriter
from config.settings import Settings

settings = Settings(
    export_directory="./batch_analysis_results"
)

# Create summary writer for all connectivities combined
writer_all = MultipleFilesSummaryWriter(settings, mode="all")
writer_all.write()

# Or create separate summaries per connectivity
writer_by_conn = MultipleFilesSummaryWriter(settings, mode="connectivity")
writer_by_conn.write()
```

#### Use Cases

- **Parameter Sweeps**: Aggregate results from analyses at different temperatures, pressures, or concentrations.
- **Trajectory Ensembles**: Combine statistics from multiple independent simulation runs.
- **Comparative Studies**: Compare metrics across different connectivity types or system compositions.
- **Data Visualization**: Generate files formatted for easy plotting with external tools.

### PerformanceWriter

Writes performance metrics to JSON files with datetime serialization support. Stores execution time, memory usage, CPU utilization, and custom metrics for analysis runs.

#### Inheritance

Inherits from `BaseWriter`.

#### Initialization

```python
PerformanceWriter(settings: Settings)
```
Creates a performance data writer with specified configuration settings.

##### Parameters

- `settings` (`Settings`): Configuration settings specifying export directory and output paths.

***

#### Attributes

- `_settings` (`Settings`): Configuration settings controlling file output location.

***

#### Methods

- `write(performance: Performance) -> None`
  - Writes a `Performance` object to a JSON file with proper datetime serialization.
  - **Parameters**:
    - `performance` (`Performance`): Performance metrics object to serialize and write.
  - **Returns**: None
  - **Output File**: `{export_directory}/performance_{performance.name}.json`
  - **Behavior**:
    - Converts the `Performance` dataclass to a dictionary using `asdict()`.
    - Serializes all `datetime` objects to ISO 8601 format strings.
    - Handles datetime conversion in both main timestamp and history entries.
    - Writes formatted JSON with 2-space indentation for readability.
  - **JSON Structure**: Mirrors the `Performance` dataclass structure with converted datetime fields.

- `serialize_datetime(obj)` (internal helper)
  - Converts `datetime` objects to ISO 8601 string format.
  - **Parameters**:
    - `obj`: Object to potentially convert.
  - **Returns**: ISO format string if `obj` is `datetime`, otherwise returns `obj` unchanged.
  - **Use**: Enables JSON serialization of datetime objects.

#### File Format

- **Format**: JSON with 2-space indentation
- **Datetime Serialization**: ISO 8601 format (`YYYY-MM-DDTHH:MM:SS.ffffff`)
- **Encoding**: UTF-8
- **Extension**: `.json`

## Reader module (`io` module)

### ReaderFactory

Factory class for creating and managing file reader instances based on file type detection. Implements the factory pattern to provide automatic reader selection based on file extensions.

***

#### Initialization

```python
ReaderFactory(settings: Settings)
```
Creates a reader factory instance and registers all available reader types.

##### Parameters

- `settings` (`Settings`): Configuration settings containing file location and other parameters.

***

#### Attributes

- `_readers` (dict): Internal registry mapping file extensions to reader instances.
- `_settings` (`Settings`): Configuration settings shared across reader instances.

***

#### Registered Readers

The factory automatically registers the following reader types:
- **XYZReader**: Handles extended XYZ format trajectory files (`.xyz`).
- **LAMMPSReader**: Handles LAMMPS dump format trajectory files (`.lammpstrj`).

***

#### Methods

- `register_reader(reader: BaseReader) -> None`
  - Registers a new reader instance in the factory.
  - **Parameters**:
    - `reader` (BaseReader): Reader instance to register.
  - **Returns**: None
  - **Behavior**:
    - Tests reader with dummy filenames for common extensions (`.xyz`, `.lammpstrj`, `.other`).
    - Associates the first matching extension with the reader instance.
    - Allows dynamic reader registration for extensibility.

- `get_reader() -> Optional[BaseReader]`
  - Returns the appropriate reader instance for the file specified in settings.
  - **Returns**: Reader instance supporting the file type, or `None` if no compatible reader found.
  - **Behavior**:
    - Validates that the file path in settings exists.
    - Iterates through registered readers testing file detection.
    - Returns the first reader that successfully detects the file format.
  - **Raises**:
    - `ValueError`: If the file specified in `settings.file_location` does not exist.


#### Extensibility

To add support for new file formats:

1. Create a new reader class inheriting from `BaseReader`
2. Implement the `detect()` method to identify compatible files
3. Register the reader in `ReaderFactory.__init__()`:

```python
class CustomReader(BaseReader):
    def detect(self, filename: str) -> bool:
        return filename.endswith('.custom')
    
    # Implement other BaseReader methods...

# In ReaderFactory.__init__():
self.register_reader(CustomReader(settings))
```

***

#### Design Pattern

The factory pattern provides:
- **Automatic format detection**: No manual reader selection required
- **Centralized management**: All readers registered in one location
- **Easy extension**: New readers added without modifying client code
- **Type safety**: Returns `Optional[BaseReader]` for type checking

### XYZReader

Reader for parsing extended XYZ format trajectory files. Efficiently scans and indexes frames for random access, supporting large trajectory files with frame-level seeking capabilities.

#### Inheritance

Inherits from `BaseReader`.

#### Initialization

```python
XYZReader(settings: Settings)
```
Creates an XYZ file reader with specified configuration settings.

##### Parameters

- `settings` (`Settings`): Configuration settings containing file location, lattice parameters, and other options.

***

#### Attributes

- `frame_indices` (List[FrameIndex]): List of frame metadata objects containing frame IDs, node counts, lattice matrices, and byte offsets.
- `num_frames` (int): Total number of frames indexed in the trajectory file.
- `is_indexed` (bool): Flag indicating whether the file has been scanned and indexed.

***

#### Named Tuple: FrameIndex

```python
FrameIndex = namedtuple("FrameIndex", ["frame_id", "num_nodes", "lattice", "byte_offset"])
```

Stores metadata for fast frame access:
- `frame_id` (int): Sequential frame identifier.
- `num_nodes` (int): Number of atoms in the frame.
- `lattice` (np.ndarray): 3×3 lattice matrix from frame header.
- `byte_offset` (int): File position (in bytes) where the frame begins.

***

#### Methods

- `detect(filepath: str) -> bool`
  - Determines if the reader supports the specified file format.
  - **Parameters**:
    - `filepath` (str): Path to the file to check.
  - **Returns**: `True` if file has `.xyz` extension (case-insensitive), `False` otherwise.
  - **Use**: Enables automatic reader selection by `ReaderFactory`.

- `scan() -> List[FrameIndex]`
  - Scans the entire trajectory file to build an index of all frames.
  - **Returns**: List of `FrameIndex` objects for all frames in the file.
  - **Behavior**:
    - Opens file with buffered I/O for efficient sequential reading.
    - For each frame:
      - Records byte offset of frame start.
      - Parses number of atoms from first line.
      - Extracts lattice matrix from extended XYZ header.
      - Skips atomic data lines to reach next frame.
      - Creates `FrameIndex` object with metadata.
    - Updates `frame_indices`, `num_frames`, and `is_indexed` attributes.
    - Prints summary message if verbose mode enabled.
  - **Extended XYZ Format Requirements**:
    - Line 1: Number of atoms (integer).
    - Line 2: Comment line containing `Lattice="lxx lxy lxz lyx lyy lyz lzx lzy lzz"`.
    - Subsequent lines: Atomic data (symbol x y z ...).
  - **Raises**:
    - `IOError`: If frame header is malformed or parsing fails.
    - `FileNotFoundError`: If trajectory file does not exist.
  - **Performance**: Single-pass sequential scan with buffered I/O for large files.

- `parse(frame_id: int) -> Generator[Frame, None, None]`
  - Parses a specific frame from the trajectory file and yields a `Frame` object.
  - **Parameters**:
    - `frame_id` (int): Zero-based index of the frame to parse.
  - **Yields**: `Frame` object containing atomic symbols, positions, and lattice information.
  - **Behavior**:
    - Calls `scan()` if file not yet indexed.
    - Retrieves frame metadata from `frame_indices[frame_id]`.
    - Seeks directly to frame's byte offset in file.
    - Parses lattice matrix from header.
    - Applies custom lattice if `settings.lattice.apply_custom_lattice` is `True`, otherwise uses lattice from file.
    - Skips 2 header lines.
    - Reads atomic data:
      - Extracts symbol (element) and x, y, z coordinates.
      - Stores as lists in data dictionary.
    - Constructs and yields `Frame` object with parsed data.
  - **Raises**:
    - `ValueError`: If atomic data line format is invalid (must have at least 4 fields: symbol, x, y, z).
  - **Generator Pattern**: Uses generator for memory efficiency with large trajectories.

***

#### Extended XYZ Format

Expected file structure for each frame:
```
{num_atoms}
Lattice="{lxx} {lxy} {lxz} {lyx} {lyy} {lyz} {lzx} {lzy} {lzz}" Properties=...
{symbol} {x} {y} {z} [additional properties...]
{symbol} {x} {y} {z} [additional properties...]
...
```

Example:
```
3
Lattice="10.0 0.0 0.0 0.0 10.0 0.0 0.0 0.0 10.0" Properties=species:S:1:pos:R:3
Si 1.234 2.345 3.456
O 4.567 5.678 6.789
Si 7.890 8.901 9.012
```

***

#### Usage Example

```python
# Import necessary modules
from nexus import SettingsBuilder, main
import nexus.config.settings as c

# Path to the trajectory file which ends with the extension `xyz`
path = "./examples/inputs/example-SiO2-27216at.xyz"

# General settings
config_general = c.GeneralSettings(
    project_name="example-SiO2",  # Project name
    export_directory="./examples/outputs",  # Export directory
    file_location=path,  # File location
    range_of_frames=(0, 0),  # Range of frames
    apply_pbc=True,  # Apply periodic boundary conditions
    verbose=True,  # Verbose mode (if True, print title, progress bars, etc.)
    save_logs=True,  # Save logs    (save logs to export_directory/logs.txt)
    save_performance=True,  # Save performance (save performance data to export_directory/performance...json)
)

...
```

#### Performance Characteristics

- **Indexing**: O(n) single-pass scan where n is the number of frames.
- **Frame Access**: O(1) direct seek using byte offsets.
- **Memory**: Minimal memory footprint—stores only frame metadata, not atomic data.
- **File I/O**: Buffered reading for efficient disk access.
- **Scalability**: Handles trajectories with thousands of frames efficiently.

### LAMMPSReader

Reader for parsing LAMMPS trajectory dump files. Efficiently scans and indexes frames for random access, supporting multiple LAMMPS file extensions and custom column layouts.

#### Inheritance

Inherits from `BaseReader`.

#### Initialization

```python
LAMMPSReader(settings: Settings)
```
Creates a LAMMPS file reader with specified configuration settings.

##### Parameters

- `settings` (`Settings`): Configuration settings containing file location, lattice parameters, and other options.

***

#### Attributes

- `frame_indices` (List[FrameIndex]): List of frame metadata objects containing frame IDs, node counts, lattice matrices, and byte offsets.
- `num_frames` (int): Total number of frames indexed in the trajectory file.
- `is_indexed` (bool): Flag indicating whether the file has been scanned and indexed.
- `columns` (dict): Mapping of column names to their indices in atomic data lines.

***

#### Named Tuple: FrameIndex

```python
FrameIndex = namedtuple("FrameIndex", ["frame_id", "num_nodes", "lattice", "byte_offset"])
```

Stores metadata for fast frame access:
- `frame_id` (int): Sequential frame identifier.
- `num_nodes` (int): Number of atoms in the frame.
- `lattice` (np.ndarray): 3×3 lattice matrix derived from box bounds.
- `byte_offset` (int): File position (in bytes) where the frame begins.

***

#### Methods

- `detect(filepath: str) -> bool`
  - Determines if the reader supports the specified file format.
  - **Parameters**:
    - `filepath` (str): Path to the file to check.
  - **Returns**: `True` if file has `.lammpstrj`, `.lammps`, or `.data` extension (case-insensitive), `False` otherwise.
  - **Use**: Enables automatic reader selection by `ReaderFactory`.

- `scan() -> List[FrameIndex]`
  - Scans the entire LAMMPS trajectory file to build an index of all frames.
  - **Returns**: List of `FrameIndex` objects for all frames in the file.
  - **Behavior**:
    - Opens file with buffered I/O for efficient sequential reading.
    - For each frame:
      - Records byte offset of frame start.
      - Validates `ITEM: TIMESTEP` header to identify frame boundaries.
      - Parses number of atoms from `ITEM: NUMBER OF NODES` section.
      - Extracts box bounds from `ITEM: BOX BOUNDS` section.
      - Constructs orthogonal 3×3 lattice matrix from box dimensions: $l_{ii} = \text{high} - \text{low}$.
      - Parses column header from `ITEM: ATOMS` line to build column index mapping.
      - Skips atomic data lines to reach next frame.
      - Creates `FrameIndex` object with metadata.
    - Updates `frame_indices`, `num_frames`, and `is_indexed` attributes.
    - Prints summary message if verbose mode enabled.
  - **LAMMPS Format Requirements**:
    - `ITEM: TIMESTEP` marker.
    - `ITEM: NUMBER OF NODES` followed by atom count.
    - `ITEM: BOX BOUNDS` followed by 3 lines with `low high` bounds.
    - `ITEM: ATOMS` followed by column names (e.g., `id type x y z`).
  - **Raises**:
    - `IOError`: If frame header is malformed or parsing fails.
    - `FileNotFoundError`: If trajectory file does not exist.
  - **Limitations**: Assumes orthogonal (non-triclinic) simulation boxes.

- `parse(frame_id: int) -> Generator[Frame, None, None]`
  - Parses a specific frame from the LAMMPS trajectory file and yields a `Frame` object.
  - **Parameters**:
    - `frame_id` (int): Zero-based index of the frame to parse.
  - **Yields**: `Frame` object containing atomic types, positions, and lattice information.
  - **Behavior**:
    - Calls `scan()` if file not yet indexed.
    - Retrieves frame metadata from `frame_indices[frame_id]`.
    - Seeks directly to frame's byte offset in file.
    - Skips 9 header lines (timestep, number of atoms, box bounds, and atoms header).
    - Uses `columns` dictionary to locate type, x, y, z fields in atomic data.
    - Reads atomic data:
      - Extracts atom type and x, y, z coordinates using column indices.
      - Stores as lists in data dictionary.
    - Constructs and yields `Frame` object with parsed data.
  - **Raises**:
    - `ValueError`: If atomic data line format is invalid.
  - **Generator Pattern**: Uses generator for memory efficiency with large trajectories.
  - **Column Flexibility**: Adapts to different column orderings in LAMMPS dump files.

***

#### LAMMPS Dump Format

Expected file structure for each frame:
```
ITEM: TIMESTEP
{timestep}
ITEM: NUMBER OF ATOMS
{num_atoms}
ITEM: BOX BOUNDS pp pp pp
{xlo} {xhi}
{ylo} {yhi}
{zlo} {zhi}
ITEM: ATOMS id type x y z [additional properties...]
{id} {type} {x} {y} {z} ...
{id} {type} {x} {y} {z} ...
...
```

Example:
```
ITEM: TIMESTEP
0
ITEM: NUMBER OF ATOMS
3
ITEM: BOX BOUNDS pp pp pp
0.0 10.0
0.0 10.0
0.0 10.0
ITEM: ATOMS id type x y z
1 1 1.234 2.345 3.456
2 2 4.567 5.678 6.789
3 1 7.890 8.901 9.012
```

***

#### Usage Example

```python
from io.reader import LAMMPSReader
from config.settings import Settings

# Path to the trajectory file which ends with the extension `.lammpstrj`
settings = Settings(
    file_location="trajectory.lammpstrj",
    range_of_frames=(0, 100),
    verbose=True
)

# Create reader and scan file
reader = LAMMPSReader(settings)
reader.scan()

print(f"Total frames: {reader.num_frames}")
print(f"Column layout: {reader.columns}")

# Parse specific frames
for frame_id in [0, 10, 50]:
    frame_gen = reader.parse(frame_id)
    frame = next(frame_gen)
    
    print(f"Frame {frame_id}:")
    print(f"  Atoms: {len(frame.nodes)}")
    print(f"  Lattice:\n{frame.lattice}")
```

***

#### Performance Characteristics

- **Indexing**: O(n) single-pass scan where n is the number of frames.
- **Frame Access**: O(1) direct seek using byte offsets.
- **Memory**: Minimal memory footprint—stores only frame metadata, not atomic data.
- **File I/O**: Buffered reading for efficient disk access.
- **Scalability**: Handles trajectories with thousands of frames efficiently.

***

#### Limitations

- **Box Geometry**: Currently supports only orthogonal (non-triclinic) simulation boxes. Triclinic boxes with xy, xz, yz tilts are not yet supported.
- **Column Detection**: Assumes standard LAMMPS dump format with `ITEM: ATOMS` header containing column names.
- **Atom Types**: Stores atom types as strings. Conversion to element symbols requires external mapping.

## `utils` Module

### aesthetics module

Utility functions for terminal output, color gradients, and file deduplication.


#### Functions

- `print_title(__version__: str) -> None`
  - Prints an ASCII art title banner with version information to the console.
  - Uses colorama for colored terminal output (light blue for banner).
  - **Parameters**:
    - `__version__` (str): Version string of the package.
  - **Returns**: None
  - Displays formatted banner and version in the terminal.

- `print_title_to_file(__version__: str, path: str) -> None`
  - Writes an ASCII art title banner with version information to a file.
  - Identical to `print_title` but outputs to a file instead of console.
  - **Parameters**:
    - `__version__` (str): Version string of the package.
    - `path` (str): File path where the title should be written.
  - **Returns**: None
  - Creates or overwrites the specified file with the banner and version.

- `generate_color_gradient(num_iterations: int) -> List[Tuple[int, int, int]]`
  - Generates a gradient of RGB colors interpolating from red to blue.
  - Used for dynamic color updates in progress bars (tqdm).
  - **Parameters**:
    - `num_iterations` (int): Number of color steps in the gradient.
  - **Returns**: List of RGB tuples, each containing three integers (0-255).
  - **Behavior**:
    - If `num_iterations == 0`: Returns single start color (red).
    - If `num_iterations == 1`: Returns both start and end colors (red and blue).
    - Otherwise: Returns smooth HSV-interpolated gradient from red to blue.

- `remove_duplicate_lines(filepath: str) -> None`
  - Removes duplicate lines from a file while preserving line order.
  - Reads file, identifies unique lines using `OrderedDict`, and rewrites the file.
  - **Parameters**:
    - `filepath` (str): Path to the file to process.
  - **Returns**: None
  - **Use Case**: Prevents redundant data entries in output files when appending analysis results across multiple runs.
  - Preserves first occurrence of each unique line and maintains original ordering.


#### Module Exports

```python
__all__ = [
    'print_title',
    'generate_color_gradient',
    'remove_duplicate_lines',
    'print_title_to_file'
]
```

### geometry module

Numba-accelerated geometric calculations for periodic systems: coordinate transformations, distance and angle computations, and gyration radius.


#### Functions

- `wrap_position(position: np.ndarray, lattice: np.ndarray) -> np.ndarray`
  - Wraps a single position vector into the primary periodic cell defined by the lattice matrix.
  - **JIT-compiled** with Numba for high performance (nopython mode, caching, fast math).
  - **Parameters**:
    - `position` (np.ndarray): 3D Cartesian position vector to wrap.
    - `lattice` (np.ndarray): 3×3 lattice matrix defining the periodic cell.
  - **Returns**: Wrapped position vector inside the primary cell.
  - **Method**: Converts to fractional coordinates, applies modulo operation via floor subtraction, converts back to Cartesian.
  - **Reference**: [Fractional coordinates on Wikipedia](https://en.wikipedia.org/wiki/Fractional_coordinates#Relationship_between_fractional_and_Cartesian_coordinates)

- `wrap_positions(positions: np.ndarray, lattice: np.ndarray) -> np.ndarray`
  - Wraps multiple position vectors into the primary periodic cell.
  - **JIT-compiled** with Numba for vectorized high-performance computation.
  - **Parameters**:
    - `positions` (np.ndarray): N×3 array of Cartesian position vectors.
    - `lattice` (np.ndarray): 3×3 lattice matrix defining the periodic cell.
  - **Returns**: N×3 array of wrapped position vectors.
  - **Method**: Iteratively applies `wrap_position` logic to each position in the array.

- `calculate_direct_distance(position1: np.ndarray, position2: np.ndarray) -> float`
  - Computes the Euclidean distance between two positions in direct space (no periodic boundaries).
  - **JIT-compiled** with Numba.
  - **Parameters**:
    - `position1` (np.ndarray): First 3D position vector.
    - `position2` (np.ndarray): Second 3D position vector.
  - **Returns**: Distance as a float.
  - **Formula**: $d = \|\mathbf{r}_1 - \mathbf{r}_2\|$

- `calculate_pbc_distance(position1: np.ndarray, position2: np.ndarray, lattice: np.ndarray) -> float`
  - Computes the minimum image distance between two positions under periodic boundary conditions.
  - **JIT-compiled** with Numba.
  - **Parameters**:
    - `position1` (np.ndarray): First 3D position vector.
    - `position2` (np.ndarray): Second 3D position vector.
    - `lattice` (np.ndarray): 3×3 lattice matrix defining the periodic cell.
  - **Returns**: Minimum distance as a float.
  - **Method**: 
    - Converts displacement vector to fractional coordinates.
    - Applies minimum image convention: $\mathbf{f} \leftarrow \mathbf{f} - \text{round}(\mathbf{f})$.
    - Converts back to Cartesian and computes norm.
  - **Physical Significance**: Ensures distances are computed considering the nearest periodic image, essential for accurate neighbor detection in periodic systems.
  - **Reference**: [Fractional coordinates on Wikipedia](https://en.wikipedia.org/wiki/Fractional_coordinates#Relationship_between_fractional_and_Cartesian_coordinates)

- `calculate_direct_angle(position1: np.ndarray, position2: np.ndarray, position3: np.ndarray) -> float`
  - Calculates the angle formed by three positions in direct space, with `position2` as the vertex.
  - **JIT-compiled** with Numba.
  - **Parameters**:
    - `position1` (np.ndarray): First position (end of vector 1).
    - `position2` (np.ndarray): Vertex position (center).
    - `position3` (np.ndarray): Third position (end of vector 2).
  - **Returns**: Angle in degrees.
  - **Formula**: 
    $$
    \theta = \arccos\left(\frac{(\mathbf{r}_1 - \mathbf{r}_2) \cdot (\mathbf{r}_3 - \mathbf{r}_2)}{\|\mathbf{r}_1 - \mathbf{r}_2\| \|\mathbf{r}_3 - \mathbf{r}_2\|}\right)
    $$

- `calculate_pbc_angle(position1: np.ndarray, position2: np.ndarray, position3: np.ndarray, lattice: np.ndarray) -> float`
  - Calculates the angle formed by three positions under periodic boundary conditions.
  - **JIT-compiled** with Numba.
  - **Parameters**:
    - `position1` (np.ndarray): First position (end of vector 1).
    - `position2` (np.ndarray): Vertex position (center).
    - `position3` (np.ndarray): Third position (end of vector 2).
    - `lattice` (np.ndarray): 3×3 lattice matrix.
  - **Returns**: Angle in degrees.
  - **Method**: 
    - Converts displacement vectors to fractional coordinates.
    - Applies minimum image convention to both vectors.
    - Converts back to Cartesian and computes angle using dot product.
    - Clips cosine value to [-1, 1] to handle numerical precision issues.
  - **Physical Significance**: Ensures angle calculations account for nearest periodic images, critical for bond angle analysis in periodic systems.

- `cartesian_to_fractional(position: np.ndarray, lattice: np.ndarray) -> np.ndarray`
  - Converts Cartesian coordinates to fractional coordinates in the lattice basis.
  - **Parameters**:
    - `position` (np.ndarray): Cartesian position vector or array of vectors.
    - `lattice` (np.ndarray): 3×3 lattice matrix.
  - **Returns**: Fractional coordinates.
  - **Method**: Solves $\mathbf{L}^T \mathbf{f} = \mathbf{r}$ where $\mathbf{L}$ is the lattice, $\mathbf{r}$ is Cartesian, and $\mathbf{f}$ is fractional.
  - **Reference**: [Fractional coordinates on Wikipedia](https://en.wikipedia.org/wiki/Fractional_coordinates#Relationship_between_fractional_and_Cartesian_coordinates)

- `fractional_to_cartesian(position: np.ndarray, lattice: np.ndarray) -> np.ndarray`
  - Converts fractional coordinates to Cartesian coordinates.
  - **Parameters**:
    - `position` (np.ndarray): Fractional position vector or array of vectors.
    - `lattice` (np.ndarray): 3×3 lattice matrix.
  - **Returns**: Cartesian coordinates.
  - **Formula**: $\mathbf{r} = \mathbf{f} \cdot \mathbf{L}$ where $\mathbf{L}$ is the lattice matrix.
  - **Reference**: [Fractional coordinates on Wikipedia](https://en.wikipedia.org/wiki/Fractional_coordinates#Relationship_between_fractional_and_Cartesian_coordinates)

- `calculate_gyration_radius(positions: np.ndarray, center_of_mass: np.ndarray) -> float`
  - Calculates the radius of gyration for a set of positions about their center of mass.
  - **JIT-compiled** with Numba.
  - **Parameters**:
    - `positions` (np.ndarray): N×3 array of positions.
    - `center_of_mass` (np.ndarray): 3D center of mass position vector.
  - **Returns**: Radius of gyration as a float.
  - **Formula**: 
    $$
    R_g = \sqrt{\frac{1}{N} \sum_{i=1}^{N} |\mathbf{r}_i - \mathbf{r}_{\text{COM}}|^2}
    $$
  - **Physical Significance**: Quantifies the spatial extent and compactness of a cluster. Returns 0.0 for empty position arrays.


#### Module Exports

```python
__all__ = [
    'wrap_position',
    'wrap_positions',
    'calculate_direct_distance',
    'calculate_pbc_distance',
    'calculate_direct_angle',
    'calculate_pbc_angle',
    'calculate_gyration_radius',
    'cartesian_to_fractional',
    'fractional_to_cartesian'
]
```

#### Performance Optimization

Most functions use Numba's JIT compilation with:
- **nopython=True**: Forces compilation to pure machine code without Python interpreter overhead.
- **cache=True**: Caches compiled functions to disk for faster subsequent imports.
- **fastmath=True**: Enables aggressive floating-point optimizations for improved performance.
- **np.ascontiguousarray()**: Ensures memory layout compatibility for optimized operations.

### Performance

Stores and tracks performance metrics for pipeline operations. Records execution time, memory usage, and CPU usage for a named operation. Supports custom metrics and historical snapshots for trend analysis.


#### Attributes

**Required Fields:**
- `id` (str): Unique identifier for the performance record.
- `name` (str): Descriptive name of the operation or component being measured.
- `timestamp` (datetime): Timestamp of record creation. Defaults to current time via `datetime.now()`.

**Performance Metrics:**
- `execution_time_ms` (Optional[float]): Execution time in milliseconds. Default: `None`.
- `memory_usage_mb` (Optional[float]): Memory consumption in megabytes. Default: `None`.
- `cpu_usage_percent` (Optional[float]): CPU utilization as a percentage. Default: `None`.

**Extensibility:**
- `metrics` (Dict[str, Any]): Dictionary for storing additional custom metrics. Default: empty dictionary.
- `history` (List[Dict[str, Any]]): List of historical performance snapshots. Default: empty list.

#### Methods

- `add_metric(name: str, value: Any) -> None`
  - Adds a custom metric to the performance record.
  - **Parameters**:
    - `name` (str): Name of the metric to store.
    - `value` (Any): Value associated with the metric.
  - **Returns**: None
  - Stores the metric in the `metrics` dictionary, allowing flexible tracking of domain-specific performance indicators.

- `record_history() -> None`
  - Captures a snapshot of the current performance state and appends it to the history list.
  - **Returns**: None
  - **Snapshot Contents**:
    - Current timestamp
    - Execution time (if set)
    - Memory usage (if set)
    - CPU usage (if set)
    - Copy of all custom metrics
  - Useful for tracking performance evolution over multiple iterations or frames.

- `get_average_execution_time() -> Optional[float]`
  - Calculates the mean execution time from all historical records.
  - **Returns**: Average execution time in milliseconds, or `None` if no valid entries exist in history.
  - Filters out `None` values before averaging.

- `__str__() -> str`
  - Provides a human-readable string representation of the performance record.
  - **Returns**: Formatted multi-line string showing ID, name, timestamp, core metrics, and additional metrics.


#### Usage Example

```python
import nexus.config.settings as c

general_settings = c.GeneralSettings(
  ...
  save_performance=True,  # Save performance (save performance data to export_directory/performance.json)
  ...
)

...
```
