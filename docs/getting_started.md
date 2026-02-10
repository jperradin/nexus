.. _getting_started:

Getting Started
===============

### Quick Start Example

Run the following commands to start a pre-made launch script `quickstart-SiO2.py`.

```bash
# get the github repo with the example files and scripts.
git clone https://github.com/jperradin/nexus.git
cd nexus/
# install if it is not installed yet.
pip install -e .

# run the quickstart example script.
python examples/inputs/quickstart-SiO2.py

# go to the export directory to see the results.
cd examples/outputs/quickstart-SiO2
```

### Workflow Overview

This section will guide you through the initial setup and usage of Nexus-CAT with the example file: `examples/inputs/example-SiO2-1008at.xyz`.

![Nexus Workflow Diagram](https://github.com/jperradin/nexus/blob/main/docs/_image/workflow-Nexus.jpg)

Nexus follows a modular, configuration-driven workflow that separates user input from the core analysis pipeline. The framework processes trajectory files through distinct stages: configuration, data reading, frame analysis, clustering, and output generation.

The analysis pipeline consists of five main stages:

#### 0. Module importation

Import `main`, `SettingsBuilder`, `settings` from Nexus-CAT before starting.

```python
from nexus import SettingsBuilder, main
import nexus.config.settings as c 
```
```
```

#### 1. Configuration Setup and Execution

The workflow begins with the user defining analysis parameters through configuration objects:

- **GeneralSettings**: Project metadata, file locations, frame ranges, and verbosity

```python
# General settings
config_general = c.GeneralSettings(
    project_name="my_first_project",                          # Project name
    export_directory="./examples/outputs",                    # Export directory
    file_location="examples/inputs/example-SiO2-1008at.xyz",  # File location
    range_of_frames=(0, -1),                                  # Range of frames : (0, -1) analyze all the frames.
    apply_pbc=True,                                           # Apply periodic boundary conditions
    verbose=True,                                             # Verbose mode (if True, print title, progress bars, etc.)
    save_logs=True,                                           # Save logs (save logs to export_directory/logs.txt)
    save_performance=False,                                   # Save performance (save performance data to export_directory/performance...json)
)
```

- **LatticeSettings**: Simulation box parameters and periodic boundary conditions

```python
# Lattice settings
config_lattice = c.LatticeSettings(
    apply_custom_lattice=False,  # For this first project we will read lattice from trajectory file
)
```

- **ClusteringSettings**: Clustering criterion, connectivity patterns, cutoff distances, and coordination constraints

```python

# Clustering settings calling the CoordinationStrategy
config_clustering = c.ClusteringSettings(
    criterion="bond",                 # Choose a three node connectivity

    node_types=["Si", "O"],           # Node types that participate in the connectivity

    node_masses=[28.0855, 15.9994],   # Nodes' masses in reduced units
  
    connectivity=["Si", "O", "Si"],   # Cluster connectivity, the node order is important. 
                                      # In this case Si nodes are the networking nodes and O are the bridging nodes.
 
    cutoffs=[
        c.Cutoff(type1="Si", type2="Si", distance=3.50), # Si-Si distance cutoff
        c.Cutoff(type1="Si", type2="O", distance=2.30),  # Si-O  distance cutoff
        c.Cutoff(type1="O", type2="O", distance=3.05),   # O-O   distance cutoff
    ],
  
    with_coordination_number=True,         # Calling the CoordinationStrategy

    coordination_mode="O",                 # Choose to calculated the coordination number of networking nodes (i.e. Si nodes) 
                                           # as the number of nearest bridging nodes (i.e. O nodes) within the specified cutoff above
                                           # other modes : "all_types" or "same_type" or "different_type" or "<node_type>"

    coordination_range=[4, 6],             # The range of Si nodes' coordination number to consider

    # if all below are False, the program enters `default` mode and find A_z-B_z cluster connectivity with z = coordination range
    with_pairwise=True,       # if with_coordination_number is True, calculate pairwise coordination number ie 4-4, 5-5, 6-6 ...
    with_mixing=False,        # if with_coordination_number is True, calculate mixing coordination number ie 4-5, 5-6, 4-6 ...
    with_alternating=False,   # if with_coordination_number is True, calculate alternating coordination number ie 4-5, 5-6 ...

    with_printed_unwrapped_clusters=True, # Choose to print the unwrapped coordinates of the clusters.
    print_mode="connectivity",            # Choose to print the clusters in one file per connectivity (and per frame). 
                                          # other modes: "all", "connectivity", "individual", "none"
)
```

- **AnalysisSettings**: Selection of metrics to compute (cluster sizes, gyration radii, percolation probabilities, etc.)

```python
# Analysis settings
config_analysis = c.AnalysisSettings(
    with_all=True, # compute all the cluster properties
)
config_analysis.overwrite = True  # Overwrite previous results in the exported files. 
                                  # It can be disabled to append the new results to the existing files.
```

These settings are assembled using the `SettingsBuilder` pattern:

```python
settings = (
    SettingsBuilder()
    .with_general(general_settings)
    .with_lattice(lattice_settings)
    .with_clustering(clustering_settings)
    .with_analysis(analysis_settings)
    .build()
)
```

The `settings` object created with `SettingsBuilder` serves as a parameter for the `main` function.
The function will dispatched our settings to the various modules and start the next steps below.

Execute the `main` function to launch the program. 
```python
main(settings)
```

#### 2. File Reading and Frame Indexing

The `Reader` component (automatically selected based on file format) performs two operations:

- **Scanning**: The `scan_file()` method indexes the trajectory file, recording byte offsets and metadata (number of atoms, lattice parameters) for each frame
- **Parsing**: The `parse_data()` method extracts atomic positions and element symbols for specific frames on-demand

Supported formats include extended XYZ and LAMMPS dump files. The reader implements lazy loading—frames are parsed only when needed, minimizing memory usage for large trajectories.

#### 3. Frame Initialization

Each parsed frame creates a `Frame` object containing:

- **nodes**: List of `Node` objects representing atoms with positions, symbols, and connectivity
- **lattice**: 3×3 lattice matrix defining the simulation box
- **clusters**: Initially empty; populated during clustering

The `init_nodes()` method converts raw atomic data into `Node` objects, assigning unique IDs and initializing parent pointers for the union-find algorithm.

#### 4. Cluster Identification

The `ClusterStrategy` component (selected based on `criterion` in settings) identifies clusters through three steps:

**a. Neighbor Finding**
- Computes pairwise distances using periodic boundary conditions
- Identifies neighbors within specified cutoff distances
- Builds neighbor lists for each node

**b. Cluster Building**
- Applies the union-find algorithm to connect nodes meeting connectivity criteria
- For **distance criterion**: Connects nodes of specified types within cutoff distance
- For **bond criterion**: Connects nodes sharing a bridging atom (e.g., Si-O-Si linkages)
- For **coordination criterion**: Adds coordination number constraints
- For **shared criterion**: Requires minimum number of shared neighbors

**c. Cluster Property Calculation**
- Unwraps node positions across periodic boundaries
- Computes center of mass, radius of gyration, and percolation status
- Assigns cluster IDs and connectivity descriptors

The `get_clusters()` method returns the complete list of identified clusters for the frame.

#### 5. Analysis and Output

The `Analyzer` component processes clusters through multiple analyzers in parallel:

- Each analyzer (e.g., `AverageClusterSizeAnalyzer`, `PercolationProbabilityAnalyzer`) implements the `analyze(frame)` method
- Analyzers accumulate statistics across frames
- The `finalize()` method computes ensemble averages and standard deviations
- The `write()` method outputs results to CSV files in the export directory

### Frame Processing Loop

The main execution loop iterates over the specified frame range:

```python
for frame_id in range(start_frame, end_frame):
    # Parse frame data
    frame = reader.parse_data(frame_id)
    
    # Initialize nodes
    frame.init_nodes()
    
    # Find clusters
    strategy.find_neighbors()
    clusters = strategy.build_clusters()
    frame.clusters = clusters
    
    # Analyze frame
    for analyzer in analyzers:
        analyzer.analyze(frame, connectivities)
```

After all frames are processed, analyzers finalize results and write output files to the export directory.

### Output Files

Nexus generates multiple output files in the export directory depending on the analysis settings:

- **Cluster configurations**: Extended XYZ files with unwrapped atomic positions (`<export_directory>/unwrapped_clusters/`)
- **Bond files**: Connectivity information for visualization (`.bonds` format)
- **Analysis data**: CSV files with statistics (e.g., `average_cluster_size.dat`, `percolation_probability.dat`)
- **Logs**: Execution record with version and settings (`log.txt`)
- **Performance metrics**: Timing and resource usage data (`performance_*.json`)

### Key Design Features

**Modularity**: Each component (reader, strategy, analyzer) operates independently, allowing easy extension with new file formats, clustering algorithms, or analysis metrics.

**Lazy Loading**: Frames are parsed on-demand rather than loading the entire trajectory into memory, enabling analysis of multi-gigabyte files.

**Factory Patterns**: `ReaderFactory`, `StrategyFactory`, and `AnalyzerFactory` automatically select appropriate implementations based on configuration, hiding complexity from the user.

**Generator-Based Parsing**: The `parse()` method yields frames one at a time, supporting streaming workflows for arbitrarily long trajectories.

**Reproducibility**: Complete settings and version information are logged with each analysis run, ensuring results can be regenerated.

***

This workflow design prioritizes **flexibility** (multiple clustering criteria and analysis metrics), **performance** (JIT-compiled geometry functions, efficient file indexing), and **scalability** (memory-efficient frame processing).

[1](https://github.com/jperradin/nexus/blob/main/assets/workflow-Nexus.jpg../assets/workflow-Nexus.jpg)

