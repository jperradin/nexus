# Import necessary modules
from nexus import SettingsBuilder, main
import nexus.config.settings as c

# Path to the trajectory file
path = "./examples/inputs/example-SiO2-1008at.xyz"
# path = "./examples/inputs/example-SiO2-1008at-high_pressure.xyz"
# path = "./examples/inputs/example-SiO2-27216at-low_pressure.xyz"
# path = "./examples/inputs/example-SiO2-27216at.xyz" # critical point for SiO5-SiO5
# path = "./examples/inputs/example-SiO2-96000at.xyz"
# path = "./examples/inputs/example-SiO2-1056000at.xyz"

# General settings
config_general = c.GeneralSettings(
    project_name="quickstart-SiO2",         # Project name
    export_directory="./examples/outputs",  # Export directory
    file_location=path,                     # File location
    range_of_frames=(0, -1),                # Range of frames : (0, -1) analyze all the frames.
    apply_pbc=True,                         # Apply periodic boundary conditions
    verbose=True,                           # Verbose mode (if True, print title, progress bars, etc.)
    save_logs=True,                         # Save logs (save logs to export_directory/logs.txt)
    save_performance=False,                 # Save performance (save performance data to export_directory/performance...json)
)

# Lattice settings
config_lattice = c.LatticeSettings(
    apply_custom_lattice=False,  # If False, read lattice from trajectory file
)

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

# Analysis settings
config_analysis = c.AnalysisSettings(
    with_all=True,
)
config_analysis.overwrite = False  # Overwrite previous results

# Build Settings object
settings = (
    SettingsBuilder()
    .with_general(config_general)  # General settings \
    .with_lattice(config_lattice)  # Lattice settings \
    .with_clustering(config_clustering)  # Clustering settings \
    .with_analysis(config_analysis)  # Analysis settings \
    .build()  # Don't forget to build the settings object
)

# Run the main function to process the trajectory
main(settings)
