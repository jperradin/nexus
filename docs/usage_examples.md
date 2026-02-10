# Usage Examples

This section provides practical examples of how to use Nexus-CAT for cluster analysis with the new API that uses the builder pattern and configuration dataclasses.

### Example 1: Silica Analysis with Coordination Strategy

This example demonstrates how to analyze Si-O-Si bonds in a silica trajectory, focusing on coordination numbers between 4 and 6.

```python
# Import necessary modules
from nexus import SettingsBuilder, main
import nexus.config.settings as c

# Path to the trajectory file
path = 'examples/inputs/example-SiO2-27216at.xyz'

# General settings
config_general = c.GeneralSettings(
    project_name="SiO2",                    # Project name
    export_directory="examples/exports",    # Export directory
    file_location=path,                     # File location
    range_of_frames=(0, 10),                # Range of frames
    apply_pbc=True,                         # Apply periodic boundary conditions
    verbose=True,                           # Verbose mode (print title, progress bars, etc.)
    save_logs=True,                         # Save logs
    save_performance=True                   # Save performance metrics
)

# Lattice settings
config_lattice = c.LatticeSettings(
    apply_custom_lattice=False,             # Read lattice from trajectory file
)

# Clustering settings
config_clustering = c.ClusteringSettings(
    criterion="bond",                        # Use bond criterion for clustering
    node_types=["Si", "O"],                 # Types of nodes in the system
    node_masses=[28.0855, 15.9994],         # Masses of nodes in reduced units
    connectivity=["Si", "O", "Si"],         # Connectivity pattern to analyze
    cutoffs=[c.Cutoff(type1="Si", type2="Si", distance=3.50), # Cutoff distances
             c.Cutoff(type1="Si", type2="O", distance=2.30),
             c.Cutoff(type1="O", type2="O", distance=3.05)],

    with_coordination_number=True,          # Calculate coordination numbers
    coordination_mode="O",                  # Mode for coordination calculations
    coordination_range=[4, 6],              # Range of coordination numbers

    with_alternating=True,                  # Calculate alternating coordination

    with_number_of_shared=False,            # Don't calculate shared neighbors

    with_printed_unwrapped_clusters=False,  # Don't print cluster coordinates
    print_mode="connectivity"               # Mode for printing clusters
)

# Analysis settings
config_analysis = c.AnalysisSettings(
    with_all=True,                          # Enable all analysis methods
)

# Build Settings object
settings = (SettingsBuilder()
    .with_general(config_general)           # General settings
    .with_lattice(config_lattice)           # Lattice settings
    .with_clustering(config_clustering)     # Clustering settings
    .with_analysis(config_analysis)         # Analysis settings
    .build()                                # Build the settings object
)

# Run the main function to process the trajectory
main(settings)
```

### Example 2: Reconfiguring and Rerunning with Shared Neighbors

This example shows how to change settings and rerun the analysis with different parameters:

```python
# After running an initial analysis, you can modify settings and run again:

# Update clustering settings to focus on specific coordination
config_clustering.with_number_of_shared = True  # Now calculate shared neighbors
config_clustering.coordination_range = [6, 6]   # Look only for SiO6-SiO6 connections
config_clustering.shared_mode="O",              # O atoms as shared neighbors (bridges between Si)
config_clustering.shared_threshold=2,           # threshold for the number of shared O atoms between Si
config_clustering.shared_threshold_mode="exact" # exact match for shared neighbors
config_analysis.overwrite = False               # Preserve previous results

# Rebuild settings with updated configuration
settings = (SettingsBuilder()
    .with_general(config_general)
    .with_lattice(config_lattice)
    .with_clustering(config_clustering)
    .with_analysis(config_analysis)
    .build()
)

# Run the analysis again with the new settings
main(settings)
```

### Example 3: Distance-Based Clustering

This example demonstrates how to use the simpler distance criterion instead of the bond criterion, analyzing direct Si-Si connectivity in silica:

```python
# Import necessary modules
from nexus import SettingsBuilder, main
import nexus.config.settings as c

# Path to the trajectory file
path = './examples/inputs/example-SiO2-1008at.xyz'

# General settings
config_general = c.GeneralSettings(
    project_name="SiO2-distance",           # Project name
    export_directory="examples/exports",    # Export directory
    file_location=path,                     # File location
    range_of_frames=(0, -1),                # All frames (-1 means to the end)
    apply_pbc=True,                         # Apply periodic boundary conditions
    verbose=True,                           # Verbose mode
    save_logs=True,                         # Save logs
    save_performance=False                  # Don't save performance metrics
)

# Lattice settings
config_lattice = c.LatticeSettings(
    apply_custom_lattice=False,             # Read lattice from trajectory file
)

# Clustering settings using the distance criterion
config_clustering = c.ClusteringSettings(
    criterion="distance",                    # Use distance criterion for clustering
    node_types=["Si", "O"],                 # Types of nodes in the system
    node_masses=[28.0855, 15.9994],         # Masses of nodes in reduced units
    connectivity=["Si", "Si"],              # Study direct Si-Si connectivity
    cutoffs=[c.Cutoff(type1="Si", type2="Si", distance=3.50), # Cutoff distances
             c.Cutoff(type1="Si", type2="O", distance=2.30),
             c.Cutoff(type1="O", type2="O", distance=3.05)]
)

# Analysis settings
config_analysis = c.AnalysisSettings(
    with_all=True,                          # Enable all analysis methods
)

# Build Settings object
settings = (SettingsBuilder()
    .with_general(config_general)           # General settings
    .with_lattice(config_lattice)           # Lattice settings
    .with_clustering(config_clustering)     # Clustering settings
    .with_analysis(config_analysis)         # Analysis settings
    .build()                                # Build the settings object
)

# Run the main function to process the trajectory
main(settings)
```

### Example 4: Custom Lattice Configuration

This example shows how to specify a custom lattice rather than reading from the trajectory file:

```python
import numpy as np
from nexus import SettingsBuilder, main
import nexus.config.settings as c

# Define a custom lattice (3x3 matrix)
custom_lattice = np.array([
    [10.0, 0.0, 0.0],
    [0.0, 10.0, 0.0],
    [0.0, 0.0, 10.0]
])

# Lattice settings with custom lattice
config_lattice = c.LatticeSettings(
    apply_custom_lattice=True,              # Use the custom lattice
    custom_lattice=custom_lattice,          # Custom lattice matrix
)

# Other settings as needed...
config_general = c.GeneralSettings(
    project_name="custom_lattice",
    export_directory="examples/exports",
    file_location="examples/inputs/example-SiO2-1008at.xyz",
    range_of_frames=(0, 5),
    apply_pbc=True,
    verbose=True
)

# Define clustering settings
config_clustering = c.ClusteringSettings(
    criterion="bond",
    node_types=["Si", "O"],
    node_masses=[28.0855, 15.9994],
    connectivity=["Si", "O", "Si"],
    cutoffs=[c.Cutoff(type1="Si", type2="Si", distance=3.50),
             c.Cutoff(type1="Si", type2="O", distance=2.30),
             c.Cutoff(type1="O", type2="O", distance=3.05)],
)

config_analysis = c.AnalysisSettings(
    with_all=True
)

# Build and run as usual
settings = (SettingsBuilder()
    .with_general(config_general)
    .with_lattice(config_lattice)
    .with_clustering(config_clustering)
    .with_analysis(config_analysis)
    .build()
)

main(settings)
```

These examples demonstrate various ways to use the Nexus-CAT package with its modern API. For more specific scenarios or advanced features, refer to the API Reference section.
