# Import necessary modules
from nexus import SettingsBuilder, main
import nexus.config.settings as c

# Path to the trajectory file
# path = "./example-a_H2O-8192mol-0.01kbar.xyz"
path = "./example-a_H2O-8192mol-0.5kbar.xyz"
# path = "./example-a_H2O-8192mol-1kbar.xyz"

# General settings
config_general = c.GeneralSettings(
    project_name=path.replace('.xyz',''),         # Project name
    export_directory="../outputs",  # Export directory
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
    criterion="distance",                 # Choose a three node connectivity

    node_types=["O"],           # Node types that participate in the connectivity

    node_masses=[15.9994],   # Nodes' masses in reduced units
  
    connectivity=["O", "O"],   # Cluster connectivity, the node order is important. 
                                      # In this case Si nodes are the networking nodes and O are the bridging nodes.
 
    cutoffs=[
        c.Cutoff(type1="O", type2="O", distance=3.5),   # O-O   distance cutoff
    ],
  
    with_coordination_number=True,         # Calling the CoordinationStrategy

    coordination_mode="O",                 # Choose to calculated the coordination number of networking nodes (i.e. O nodes) 

    coordination_range=[4, 4],             # The range of O nodes' coordination number to consider : Z = 4

    # if all below are False, the program enters `default` mode and find A_z-B_z cluster connectivity with z = coordination range
    with_pairwise=False,       # if with_coordination_number is True, calculate pairwise coordination number ie 4-4, 5-5, 6-6 ...
    with_mixing=False,        # if with_coordination_number is True, calculate mixing coordination number ie 4-5, 5-6, 4-6 ...
    with_alternating=False,   # if with_coordination_number is True, calculate alternating coordination number ie 4-5, 5-6 ...
    with_connectivity_name='LD',

    with_printed_unwrapped_clusters=True, # Choose to print the unwrapped coordinates of the clusters.
    print_mode="connectivity",            # Choose to print the clusters in one file per connectivity (and per frame). 
                                          # other modes: "all", "connectivity", "individual", "none"
)

# Analysis settings
config_analysis = c.AnalysisSettings(
    with_all=True,
)
config_analysis.overwrite = True # Overwrite previous results

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


#### Reconfigure to compute the HD phase
# Clustering settings calling the CoordinationStrategy
config_clustering = c.ClusteringSettings(
    criterion="distance",                 # Choose a three node connectivity

    node_types=["O"],           # Node types that participate in the connectivity

    node_masses=[15.9994],   # Nodes' masses in reduced units
  
    connectivity=["O", "O"],   # Cluster connectivity, the node order is important. 
                                      # In this case Si nodes are the networking nodes and O are the bridging nodes.
 
    cutoffs=[
        c.Cutoff(type1="O", type2="O", distance=3.5),   # O-O   distance cutoff
    ],
  
    with_coordination_number=True,         # Calling the CoordinationStrategy

    coordination_mode="O",                 # Choose to calculated the coordination number of networking nodes (i.e. O nodes) 

    coordination_range=[5, 7],             # The range of O nodes' coordination number to consider : 4 < Z < 8

    # if all below are False, the program enters `default` mode and find A_z-B_z cluster connectivity with z = coordination range
    with_pairwise=False,       # if with_coordination_number is True, calculate pairwise coordination number ie 4-4, 5-5, 6-6 ...
    with_mixing=False,        # if with_coordination_number is True, calculate mixing coordination number ie 4-5, 5-6, 4-6 ...
    with_alternating=False,   # if with_coordination_number is True, calculate alternating coordination number ie 4-5, 5-6 ...
    with_connectivity_name='HD',

    with_printed_unwrapped_clusters=True, # Choose to print the unwrapped coordinates of the clusters.
    print_mode="connectivity",            # Choose to print the clusters in one file per connectivity (and per frame). 
                                          # other modes: "all", "connectivity", "individual", "none"
)

config_analysis.overwrite = False  # Do not overwrite previous results, simply append the results to the existing files.

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

#### Reconfigure to compute the VHD phase
# Clustering settings calling the CoordinationStrategy
config_clustering = c.ClusteringSettings(
    criterion="distance",                 # Choose a three node connectivity

    node_types=["O"],           # Node types that participate in the connectivity

    node_masses=[15.9994],   # Nodes' masses in reduced units
  
    connectivity=["O", "O"],   # Cluster connectivity, the node order is important. 
                                      # In this case Si nodes are the networking nodes and O are the bridging nodes.
 
    cutoffs=[
        c.Cutoff(type1="O", type2="O", distance=3.5),   # O-O   distance cutoff
    ],
  
    with_coordination_number=True,         # Calling the CoordinationStrategy

    coordination_mode="O",                 # Choose to calculated the coordination number of networking nodes (i.e. O nodes) 

    coordination_range=[8, 999],           # The range of O nodes' coordination number to consider : Z >= 8

    # if all below are False, the program enters `default` mode and find A_z-B_z cluster connectivity with z = coordination range
    with_pairwise=False,       # if with_coordination_number is True, calculate pairwise coordination number ie 4-4, 5-5, 6-6 ...
    with_mixing=False,        # if with_coordination_number is True, calculate mixing coordination number ie 4-5, 5-6, 4-6 ...
    with_alternating=False,   # if with_coordination_number is True, calculate alternating coordination number ie 4-5, 5-6 ...
    with_connectivity_name='VHD',

    with_printed_unwrapped_clusters=True, # Choose to print the unwrapped coordinates of the clusters.
    print_mode="connectivity",            # Choose to print the clusters in one file per connectivity (and per frame). 
                                          # other modes: "all", "connectivity", "individual", "none"
)

config_analysis.overwrite = False  # Do not overwrite previous results, simply append the results to the existing files.

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
