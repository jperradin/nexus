# Import necessary modules
from nexus import SettingsBuilder, main
import nexus.config.settings as c
import nexus.io.parser.parser as p
import numpy as np

# Lattice settings
config_lattice = c.LatticeSettings(
    apply_custom_lattice=False,  # If False, read lattice from trajectory file
)

# Clustering settings
config_clustering = c.ClusteringSettings(
    criterion="distance",
    node_types=["1"],
    node_masses=[1.0],
    connectivity=["1", "1"],
    cutoffs=[
        c.Cutoff(type1="1", type2="1", distance=1.1)
    ],  # cutoff distance in reduced units
    # with_printed_unwrapped_clusters=True,
    # print_mode="all", # "all", "connectivity", "individual", "none"
)

# Analysis settings
config_analysis = c.AnalysisSettings(
    with_all=True,
)
config_analysis.overwrite = False

# Path to the trajectory file
j = range(85, 100, 5)
for j in j:
    path = f"benchmarks/perco-{j}_0.3116.xyz"
    project_name = f"{j}"
    config_general = c.GeneralSettings(
        project_name=project_name,
        export_directory=f"benchmarks/outputs/",
        file_location=path,
        range_of_frames=(0, -1),
        apply_pbc=True,
        verbose=True,
        save_logs=True,
        save_performance=True,
    )
    # Settings builder
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
