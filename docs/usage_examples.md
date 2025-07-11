Usage Examples
==============

This section provides practical examples of how to use Nexus-CAT for cluster analysis with the new API that uses the builder pattern and configuration dataclasses.

### Example 1: Silica Analysis

This example demonstrates how to analyze Si-O-Si bonds in a silica trajectory, focusing on coordination numbers between 4 and 6.

```python
# Import necessary modules
from nexus import SettingsBuilder, main
import nexus.config.settings as c

# Path to the trajectory file
path = 'examples/inputs/SiO2-27216at-pos67B.xyz'

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
    criteria="bond",                        # Use bond criteria for clustering
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
    shared_mode="O",                        
    shared_threshold=2,                     

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

### Example 2: Reconfiguring and Rerunning

This example shows how to change settings and rerun the analysis with different parameters:

```python
# After running an initial analysis, you can modify settings and run again:

# Update clustering settings to focus on specific coordination
config_clustering.with_number_of_shared = True  # Now calculate shared neighbors
config_clustering.coordination_range = [6, 6]   # Look only for SiO6-SiO6 connections
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

### Example 3: Water Analysis

This example demonstrates how to analyze water/ice systems:

```python
# Import necessary modules
from nexus import SettingsBuilder, main
import nexus.config.settings as c

# Path to the trajectory file
path = './examples/inputs/waterM2825-5kbar.xyz'

# General settings
config_general = c.GeneralSettings(
    project_name="H2O",                     # Project name
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

# Define H2O clustering settings
config_clustering = c.ClusteringSettings(
    criteria="distance",                    # Use distance criteria for clustering 
    node_types=["O", "H"],                  # Types of nodes in the system
    node_masses=[15.9994, 1.00794],         # Masses of nodes in reduced units
    connectivity=["O", "O"],                # Study O-O connectivity
    cutoffs=[c.Cutoff(type1="O", type2="O", distance=3.5), # Cutoff distances
             c.Cutoff(type1="O", type2="H", distance=1.2),
             c.Cutoff(type1="H", type2="H", distance=2.0)]
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

### Example 4: Processing Multiple Files with the Parser

This example shows how to use the Parser to process multiple related files:

```python
# Import necessary modules
from nexus import SettingsBuilder, main
import nexus.config.settings as c
import nexus.io.parser.parser as p

# Lattice settings (same for all files)
config_lattice = c.LatticeSettings(
    apply_custom_lattice=False,             # Read lattice from trajectory file
)

# Define clustering settings for percolation analysis
config_clustering = c.ClusteringSettings(
    criteria="distance",                    # Use distance criteria
    node_types=["1"],                       # Simple node type
    node_masses=[1.0],                      # Node mass
    connectivity=["1", "1"],                # Analyze connectivity between same type
    cutoffs=[c.Cutoff(type1="1", type2="1", distance=1.1)],  # Cutoff distance
)

# Analysis settings
config_analysis = c.AnalysisSettings(
    with_all=True,                          # Enable all analysis methods
)

# Use the parser to get multiple files from a directory
rootdir = './examples/inputs/ordinary_percolation'
parser = p.Parser(file_location=rootdir, format='xyz')
files = parser.get_files()
infos = parser.get_infos()

# Process each file
for i, file in enumerate(files):
    path = file
    project_name = infos['project_name'][i]
    
    # Configure general settings for this file
    config_general = c.GeneralSettings(
        project_name=project_name,
        export_directory='examples/exports/ordinary_percolation',
        file_location=path,
        range_of_frames=(0, 0),             # Only first frame
        apply_pbc=True,
        verbose=True,
        save_logs=True
    )
    
    # Build settings for this file
    settings = (SettingsBuilder()
        .with_general(config_general)       # General settings
        .with_lattice(config_lattice)       # Lattice settings
        .with_clustering(config_clustering) # Clustering settings
        .with_analysis(config_analysis)     # Analysis settings
        .build()
    )

    # Run the analysis for this file
    main(settings)

# Generate a summary of all processed files
from nexus.io.writer.writer_factory import WriterFactory
writer_factory = WriterFactory(settings)
writer = writer_factory.get_writer("MultipleFilesSummaryWriter", mode="connectivity")
writer.write()
```

### Example 5: Custom Lattice Configuration

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
    file_location="path/to/trajectory.xyz",
    range_of_frames=(0, 5),
    apply_pbc=True,
    verbose=True
)

# Define clustering settings
config_clustering = c.ClusteringSettings(
    # Your clustering settings...
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
