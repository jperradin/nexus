API Reference
=============

This section provides some detailed reference for the Nexus-CAT API.
For a complete list of classes and methods, refer to the source code.

## main.py

### `main` function

Analyzes a trajectory based on the provided settings.

**Methods:**
- `main(settings) -> None`: The main function that initializes and runs the application.

## Modules

### `core` module

#### `Atom` class

Represents an atom within a system.

**Attributes:**
- `element (str)`: Atomic element.
- `id (int)`: Identifier of the atom in the system.
- `position (np.array)`: XYZ coordinates.
- `frame (int)`: Frame index in the trajectory.
- `cutoffs (dict)`: Cutoff distances dictionary (Cutoff object).
- `extension (str)`: Extension used for method determination.
- `atomic_mass (float)`: Atomic mass of the atom.
- `neighbours (list)`: List of first neighbours (PBC applied).
- `coordination (int)`: Number of neighbours around the atom (PBC applied).
- `parent (Atom)`: Root of the cluster.
- `cluster_id (int)`: Cluster index to which the atom belongs.
- `cluster_type (str)`: Cluster type to which the atom belongs.

**Methods:**
- `__init__(self, element, id, position, frame, cutoffs, extension="SiOz")`: Initializes an Atom object.
- `get_element(self) -> str`: Returns the element of the Atom.
- `get_id(self) -> int`: Returns the unique identifier of the Atom.
- `get_position(self) -> np.array`: Returns the spatial coordinates of the Atom.
- `get_frame(self) -> int`: Returns the frame index associated with the Atom.
- `get_neighbours(self) -> list`: Returns the list of neighbours of the Atom.
- `get_atomic_mass(self) -> float`: Returns the atomic mass of the Atom.
- `get_coordination(self) -> int`: Returns the coordination number of the Atom.
- `get_cluster_id(self) -> int`: Returns the cluster id that the atom belongs to.
- `add_neighbour(self, neighbour) -> None`: Adds a neighbour to the list of neighbours of the Atom.
- `filter_neighbours(self, distances) -> None`: Removes neighbours not within cutoff distances (depending on pair of atoms).
- `set_cluster(self, cluster_id, cluster_type) -> None`: Sets the cluster and type that the atom belongs to.
- `reset_cluster(self) -> None`: Resets the cluster(s) of the atom.

#### `Box` class

Represents a simulation box in three-dimensional space at each frame of the trajectory.

**Attributes:**
- `length_x (list)`: Length of the box in the x-direction.
- `length_y (list)`: Length of the box in the y-direction.
- `length_z (list)`: Length of the box in the z-direction.
- `volume (list)`: Volume of the box.

**Methods:**
- `__init__(self) -> None`: Initializes a Box object.
- `add_box(self, length_x, length_y, length_z) -> None`: Adds a box to the list of boxes.
- `get_volume(self, configuration) -> float`: Calculates and returns the volume of the box.
- `get_box_dimensions(self, configuration) -> list`: Returns the dimensions of the box.

#### `Cluster` class

Represents a cluster of atoms within a system.

**Attributes:**
- `atoms (list)`: List of Atom objects belonging to the cluster.
- `box (Box)`: Box object representing the simulation box.
- `connectivity (str)`: Connectivity of the cluster.
- `root_id (int)`: Atom id that is the root of the cluster.
- `frame (int)`: Frame of the trajectory.
- `size (int)`: Size of the cluster (number of atoms).
- `center_of_mass (list)`: Center of mass of the cluster.
- `indices (list)`: List of indices of the atoms.
- `unwrapped_positions (list)`: List of unwrapped positions of the cluster.
- `percolation_probability (str)`: Percolation probability.
- `gyration_radius (float)`: Gyration radius of the cluster.

**Methods:**
- `__init__(self, atoms=None, box=None, connectivity="", root_id=None, frame=None, size=None) -> None`: Initializes a Cluster object.
- `add_atom(self, atom) -> None`: Adds an atom to the cluster.
- `get_atoms(self) -> list`: Returns the list of Atom objects belonging to the cluster.
- `set_indices_and_positions(self, positions_dict) -> None`: Sets the array of unique indices and positions of atoms in the cluster.
- `calculate_center_of_mass(self) -> None`: Calculates the center of mass of the cluster.
- `write_coordinates(self, path_to_directory) -> None`: Writes the cluster coordinates to an XYZ file.
- `calculate_gyration_radius(self) -> None`: Calculates the gyration radius of the cluster.
- `calculate_percolation_probability(self) -> None`: Calculates the percolation probability of the cluster.
- `calculate_unwrapped_positions(self, criteria, chain, quiet=False) -> None`: Calculates the unwrapped positions of atoms in the cluster.
- `unwrap_position(self, vector, box_size)`: Unwraps position considering periodic boundary conditions.

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

Represents a system of atoms and provides methods for analyzing and manipulating the system.

**Attributes:**
- `settings (Settings)`: Settings object containing the list of all the parameters.
- `atoms (list)`: List of all the atoms in the system.
- `box (Box)`: The Box object containing the lattice information at each frame.
- `clusters (list)`: List of all the clusters of the system.
- `counter_c (int)`: Counter of Cluster object.
- `frame (int)`: Frame of the system in the trajectory.
- `cutoffs (Cutoff)`: Cutoff object managing cutoff distances for pairs of elements.

**Methods:**
- `__init__(self, settings) -> None`: Initializes a System object.
- `add_atom(self, atom) -> None`: Adds an Atom object to the list of atoms.
- `add_cluster(self, cluster:object) -> None`: Adds a Cluster object to the list of clusters.
- `get_atoms(self) -> list`: Returns the list of atoms.
- `get_positions(self) -> tuple`: Returns the list of positions and elements of all Atom objects.
- `get_positions_by_element(self, element) -> np.array`: Returns the list of positions of all Atom objects of the same element.
- `get_atoms_by_element(self, element) -> list`: Returns the list of Atom objects belonging to the same species.
- `get_unique_element(self) -> np.array`: Returns the unique elements present in the system along with their counts.
- `reset_cluster_indexes(self) -> None`: Resets the cluster indexes for all Atom objects in the system.
- `wrap_atomic_positions(self) -> None`: Wraps atomic positions inside the simulation box using periodic boundary conditions.
- `compute_mass(self) -> float`: Returns the mass of the system in atomic unit.
- `calculate_neighbours(self) -> None`: Calculates the nearest neighbours of all atoms in the system.
- `calculate_concentrations(self, extension) -> None`: Determines the structural units and other structural properties.
- `set_concentrations(self, connectivity: str) -> None`: Sets the concentrations of the structural units.
- `get_concentration(self, connectivity: str) -> float`: Returns the concentration of the sites for a given connectivity.
- `get_all_clusters(self, connectivity: str) -> list`: Returns the list of all Cluster objects associated with the given connectivity.
- `get_filtered_clusters(self, connectivity:str) -> list`: Returns the list of Cluster objects associated with the given connectivity.
- `get_all_cluster_sizes(self, connectivity:str) -> list`: Returns the list of cluster sizes associated with the given connectivity.
- `get_filtered_cluster_sizes(self, connectivity:str) -> list`: Returns the list of cluster sizes associated with the given connectivity.
- `get_cluster_sizes_distribution(self, connectivity:str) -> dict`: Returns the cluster size distribution of the clusters associated with the given connectivity.
- `get_gyration_radius_distribution(self, connectivity:str, list_sizes:list) -> dict`: Returns the gyration radius distribution of the clusters associated with the given connectivity.
- `calculate_order_parameter(self, connectivity:str) -> list`: Calculates the order parameter of the percolating cluster.
- `calculate_percolation_probability(self, connectivity:str) -> list`: Calculates the percolation probability of the percolating cluster.
- `write_coordinates_all_in_one(self, connectivity, path_to_directory) -> None`: Writes the cluster coordinates to an XYZ file.
- `decrypt_connectivity(self, connectivity)`: Decrypts the connectivity string to get the atomic species and the number of neighbors.
- `find(self, atom) -> object`: Finds the root of the cluster to which the given atom belongs.
- `union(self, atom_1, atom_2) -> None`: Unions the two clusters to which the given atoms belong.
- `find_clusters(self, connectivity) -> None`: Finds clusters based on specified criteria.
- `find_extra_clusters(self) -> None`: Finds extra clusters based on the extension.

### `data` module

#### `chemical_symbols` array

An array containing the chemical symbols of elements in the periodic table.

#### `atomic_masses` array

An array containing the atomic masses of elements in the periodic table.

#### `load_data` function

Loads the chemical symbols and atomic masses data.

**Methods:**
- `load_data() -> tuple`: Loads the chemical symbols and atomic masses data.

#### `get_chemical_symbol` function

Returns the chemical symbol for a given atomic number.

**Methods:**
- `get_chemical_symbol(atomic_number: int) -> str`: Returns the chemical symbol for a given atomic number.

#### `get_atomic_mass` function

Returns the atomic mass for a given atomic number.

**Methods:**
- `get_atomic_mass(atomic_number: int) -> float`: Returns the atomic mass for a given atomic number.

### `extensions` module

#### `SiSi` module

This module contains all the methods/functions that are specific to Si-Si clusters.

**Classes:**
- `Silicon`: Represents a silicon atom within a Si-Si cluster.

**Functions:**
- `transform_into_subclass(atom: Atom) -> object`: Returns a Silicon object from the subclass Silicon.
- `get_connectivity(cluster_settings) -> list`: Returns the list of connectivity for the given cluster settings.
- `get_extra_connectivity(cluster_settings) -> list`: Returns the list of extra connectivity for the given cluster settings.
- `get_default_settings(criteria="distance") -> dict`: Loads the default parameters for the SiSi extension.
- `calculate_concentrations(atoms: list, criteria: str, quiet: bool) -> dict`: Calculates the concentrations for each cluster connectivity.
- `find_extra_clusters(atoms: list, box: Box, counter_c: int, settings: object) -> None`: Finds LD, HD, VHD, and HV clusters in the system.

#### `SiOz` module

This module contains all the methods/functions that are specific to SiOz-SiOz clusters.

**Classes:**
- `Silicon`: Represents a silicon atom within a SiOz cluster.
- `Oxygen`: Represents an oxygen atom within a SiOz cluster.

**Functions:**
- `transform_into_subclass(atom: Atom) -> object`: Returns a Silicon or Oxygen object from the subclass Silicon or Oxygen.
- `get_connectivity(cluster_settings) -> list`: Returns the list of connectivity for the given cluster settings.
- `get_extra_connectivity(cluster_settings) -> list`: Returns the list of extra connectivity for the given cluster settings.
- `get_default_settings(criteria="bond") -> dict`: Loads the default parameters for the SiOz extension.
- `calculate_concentrations(atoms: list, criteria: str, quiet: bool) -> dict`: Calculates the concentrations for each cluster connectivity.
- `find_extra_clusters(atoms: list, box: Box, counter_c: int, settings: object) -> None`: Finds stishovite clusters in the system.

#### `OO` module

This module contains all the methods/functions that are specific to O-O clusters.

**Classes:**
- `Oxygen`: Represents an oxygen atom within an O-O cluster.

**Functions:**
- `transform_into_subclass(atom: Atom) -> object`: Returns an Oxygen object from the subclass Oxygen.
- `get_connectivity(cluster_settings) -> list`: Returns the list of connectivity for the given cluster settings.
- `get_extra_connectivity(cluster_settings) -> list`: Returns the list of extra connectivity for the given cluster settings.
- `get_default_settings(criteria="distance") -> dict`: Loads the default parameters for the OO extension.
- `calculate_concentrations(atoms: list, criteria: str, quiet: bool) -> dict`: Calculates the concentrations for each cluster connectivity.
- `find_extra_clusters(atoms: list, box: Box, counter_c: int, settings: object) -> None`: Finds LD, HD, VHD, and HV clusters in the system.

#### `NaO` module

This module contains all the methods/functions that are specific to Na-Na clusters.

**Classes:**
- `Sodium`: Represents a sodium atom within a Na-Na cluster.
- `Oxygen`: Represents an oxygen atom within a Na-Na cluster.

**Functions:**
- `transform_into_subclass(atom: Atom) -> object`: Returns a Sodium or Oxygen object from the subclass Sodium or Oxygen.
- `get_connectivity(cluster_settings) -> list`: Returns the list of connectivity for the given cluster settings.
- `get_extra_connectivity(cluster_settings) -> list`: Returns the list of extra connectivity for the given cluster settings.
- `get_default_settings(criteria="distance") -> dict`: Loads the default parameters for the NaO extension.
- `calculate_concentrations(atoms: list, criteria: str, quiet: bool) -> dict`: Calculates the concentrations for each cluster connectivity.
- `find_extra_clusters(atoms: list, box: Box, counter_c: int, settings: object) -> None`: Finds LD, HD, VHD, and HV clusters in the system.

### `io` module

#### `read_lattices_properties` function

Reads lattice properties from a trajectory file and stores them in a Box object.

**Methods:**
- `read_lattices_properties(box, file_path, keyword="Lattice") -> None`: Creates the Box object for each frame in the trajectory file.

#### `count_configurations` function

Counts the number of configurations in a trajectory file.

**Methods:**
- `count_configurations(file_path, keyword="Lattice") -> int`: Counts the number of configurations in the trajectory file.

#### `read_and_create_system` function

Reads an XYZ file and creates a System object from the data.

**Methods:**
- `read_and_create_system(file_path, frame, frame_size, settings, cutoffs, start, end) -> System`: Reads the XYZ file and returns the frame as a System object.

#### `write_list_of_files` function

Writes a list of all-in-one unwrapped clusters files to a text file.

**Methods:**
- `write_list_of_files(dirpath: str) -> None`: Writes a list of the all-in-one unwrapped clusters files to a text file.

#### `make_lines_unique` function

Reads a file, removes duplicate lines, and rewrites the file with unique lines.

**Methods:**
- `make_lines_unique(filepath: str) -> None`: Reads a file, removes duplicate lines, and rewrites the file with unique lines.

#### `Result` class

Represents a generic result.

**Attributes:**
- `property (str)`: The property name.
- `info (str)`: Additional information about the result.
- `init_frame (int)`: The initial frame.
- `timeline (list)`: List of data points over time.
- `result (float)`: The computed result.
- `error (float)`: The error associated with the result.

#### `AverageClusterSize` class

Represents a result for average cluster sizes.

**Attributes:**
- `average_size (float)`: The computed average cluster size.
- `filepath (str)`: The path to the output file.

**Methods:**
- `add_to_timeline(self, value, concentration) -> None`: Appends a data point to the timeline.
- `calculate_average_cluster_size(self) -> None`: Calculates the average cluster size based on the timeline data.
- `write_file_header(self, overwrite, path_to_directory, number_of_frames) -> None`: Initializes the output file with a header.
- `append_results_to_file(self) -> None`: Appends the result to the output file.

#### `BiggestClusterSize` class

Represents a result for the biggest cluster size.

**Attributes:**
- `biggest_size (float)`: The computed biggest cluster size.
- `filepath (str)`: The path to the output file.

**Methods:**
- `add_to_timeline(self, value, concentration) -> None`: Appends a data point to the timeline.
- `calculate_biggest_cluster_size(self) -> None`: Calculates the biggest cluster size based on the timeline data.
- `write_file_header(self, overwrite, path_to_directory, number_of_frames) -> None`: Initializes the output file with a header.
- `append_results_to_file(self) -> None`: Appends the result to the output file.

#### `SpanningClusterSize` class

Represents a result for the spanning cluster size.

**Attributes:**
- `spanning_size (float)`: The computed spanning cluster size.
- `filepath (str)`: The path to the output file.

**Methods:**
- `add_to_timeline(self, value, concentration) -> None`: Appends a data point to the timeline.
- `calculate_spanning_cluster_size(self) -> None`: Calculates the spanning cluster size based on the timeline data.
- `write_file_header(self, overwrite, path_to_directory, number_of_frames) -> None`: Initializes the output file with a header.
- `append_results_to_file(self) -> None`: Appends the result to the output file.

#### `ClusterSizeDistribution` class

Represents a result for the cluster size distribution.

**Attributes:**
- `distribution (dict)`: The computed cluster size distribution.
- `filepath (str)`: The path to the output file.

**Methods:**
- `add_to_timeline(self, frame: int, value: list, concentration: float) -> None`: Appends a data point to the timeline.
- `calculate_cluster_size_distribution(self) -> None`: Calculates the cluster size distribution based on the timeline data.
- `write_file_header(self, path_to_directory, number_of_frames) -> None`: Initializes the output file with a header.
- `append_results_to_file(self) -> None`: Appends the result to the output file.

#### `GyrationRadiusDistribution` class

Represents a result for the gyration radius distribution.

**Attributes:**
- `distribution (dict)`: The computed gyration radius distribution.
- `filepath (str)`: The path to the output file.

**Methods:**
- `add_to_timeline(self, frame: int, value: list, concentration: float) -> None`: Appends a data point to the timeline.
- `calculate_gyration_radius_distribution(self) -> None`: Calculates the gyration radius distribution based on the timeline data.
- `write_file_header(self, path_to_directory, number_of_frames) -> None`: Initializes the output file with a header.
- `append_results_to_file(self) -> None`: Appends the result to the output file.

#### `CorrelationLength` class

Represents a result for the correlation length.

**Attributes:**
- `corre_length (float)`: The computed correlation length.
- `filepath (str)`: The path to the output file.

**Methods:**
- `add_to_timeline(self, frame, value, concentration) -> None`: Appends a data point to the timeline.
- `calculate_correlation_length(self) -> None`: Calculates the correlation length based on the timeline data.
- `write_file_header(self, overwrite, path_to_directory, number_of_frames) -> None`: Initializes the output file with a header.
- `append_results_to_file(self) -> None`: Appends the result to the output file.

#### `OrderParameter` class

Represents a result for the order parameter.

**Attributes:**
- `order_parameter (float)`: The computed order parameter.
- `filepath (str)`: The path to the output file.

**Methods:**
- `add_to_timeline(self, frame: int, value: list, concentration: float) -> None`: Appends a data point to the timeline.
- `calculate_order_parameter(self) -> None`: Calculates the order parameter based on the timeline data.
- `write_file_header(self, overwrite, path_to_directory, number_of_frames) -> None`: Initializes the output file with a header.
- `append_results_to_file(self) -> None`: Appends the result to the output file.

#### `PercolationProbability` class

Represents a result for the percolation probability.

**Attributes:**
- `percolation_probability (float)`: The computed percolation probability.
- `filepath (str)`: The path to the output file.

**Methods:**
- `add_to_timeline(self, frame: int, value:list, concentration: float) -> None`: Appends a data point to the timeline.
- `calculate_percolation_probability(self) -> None`: Calculates the percolation probability based on the timeline data.
- `write_file_header(self, overwrite, path_to_directory, number_of_frames) -> None`: Initializes the output file with a header.
- `append_results_to_file(self) -> None`: Appends the result to the output file.

### `settings` module

#### `Settings` class

Handles the configuration settings for the application.

**Attributes:**
- `project_name (Parameter)`: Name of the project used for the output directory.
- `export_directory (Parameter)`: Parent directory where the output files will be saved.
- `build_fancy_recaps (Parameter)`: Whether to build fancy recaps of the results into a single file.
- `build_fancy_plots (Parameter)`: Whether to build fancy plots of the results into a single file.
- `path_to_xyz_file (Parameter)`: Path to the XYZ file containing the atomic coordinates.
- `number_of_atoms (Parameter)`: Number of atoms in the system.
- `number_of_frames (Parameter)`: Number of frames in the XYZ file.
- `header (Parameter)`: Number of lines in the header of the XYZ file.
- `range_of_frames (Parameter)`: Range of frames to be analysed.
- `frames_to_analyse (Parameter)`: Frames to be analysed.
- `timestep (Parameter)`: Timestep of the simulation.
- `lbox (Parameter)`: Box length.
- `temperature (Parameter)`: Temperature of the system.
- `pressure (Parameter)`: Pressure of the system.
- `version (Parameter)`: Version of the software.
- `quiet (Parameter)`: Whether to print settings or not.
- `overwrite_results (Parameter)`: Whether to overwrite files by default.
- `print_clusters_positions (Parameter)`: Whether to print the positions of the clusters.
- `max_size (Parameter)`: Maximum size of the clusters for the cluster size distribution.
- `supported_extensions (Parameter)`: List of supported extensions.

**Methods:**
- `__init__(self, extension="SiOz") -> None`: Initializes a Settings object with default settings.
- `load_default_settings(self, extension) -> None`: Loads default settings based on the extension.
- `print_settings(self) -> None`: Prints the current settings.
- `print_all_settings(self) -> None`: Prints all settings, including those not recommended for printing.

#### `Parameter` class

Represents a generic parameter with a name and value.

**Attributes:**
- `name (str)`: Name of the parameter.
- `value`: Value associated with the parameter.

**Methods:**
- `__init__(self, name, value) -> None`: Initializes a Parameter object with a name and value.
- `get_name(self) -> str`: Returns the name of the parameter.
- `get_value(self)`: Returns the value associated with the parameter.
- `set_value(self, new_value) -> None`: Sets a new value for the parameter.

#### `ClusterParameter` class

Represents a parameter specific to cluster settings.

**Attributes:**
- `name (str)`: The name of the cluster parameter.
- `value (dict)`: The value of the cluster parameter, stored as a dictionary.

**Methods:**
- `__init__(self, name, value) -> None`: Initializes a ClusterParameter object.
- `set_cluster_parameter(self, key, new_value) -> None`: Replaces a value of the settings.

### `utils` module

#### `print_title` function

Prints the title and the version of the package.

**Methods:**
- `print_title(__version__) -> None`: Prints the title and the version of the package.

#### `generate_color_gradient` function

Generates a color gradient between two colors.

**Methods:**
- `generate_color_gradient(num_iterations) -> list`: Generates a color gradient between two colors.
