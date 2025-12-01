# Import necessary modules
from nexus import SettingsBuilder, main
import nexus.config.settings as c
from multiprocessing import Pool
from tqdm import tqdm
import numpy as np


def run_distance_strategy_benchmark(task):
    path, name = task

    try:
        # Lattice settings
        config_lattice = c.LatticeSettings(
            apply_custom_lattice=False,  # If False, read lattice from trajectory file
        )

        # Clustering settings
        config_clustering = c.ClusteringSettings(
            criterion="distance",
            node_types=["1", "2"],
            node_masses=[1.0, 1.0],
            connectivity=["1", "2"],
            cutoffs=[
                c.Cutoff(type1="1", type2="2", distance=0.9)
            ],  # cutoff distance in reduced units
            # with_printed_unwrapped_clusters=True,
            # print_mode="all", # "all", "connectivity", "individual", "none"
        )

        # Analysis settings
        config_analysis = c.AnalysisSettings(
            with_all=True,
        )

        # Path to the trajectory file
        config_general = c.GeneralSettings(
            project_name=name,
            export_directory=f"benchmarks/output/distance_strategy/",
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

    except Exception as e:
        print(f"An error occurred while processing {path}: {e}")

def run_bond_strategy_benchmark(task):
    path, name = task

    try:
        # Lattice settings
        config_lattice = c.LatticeSettings(
            apply_custom_lattice=False,  # If False, read lattice from trajectory file
        )

        # Clustering settings
        config_clustering = c.ClusteringSettings(
            criterion="bond",
            node_types=["1", "2"],
            node_masses=[1.0, 1.0],
            connectivity=["1", "2", "1"],
            cutoffs=[
                c.Cutoff(type1="1", type2="2", distance=0.9)
            ],  # cutoff distance in reduced units
            # with_printed_unwrapped_clusters=True,
            # print_mode="all", # "all", "connectivity", "individual", "none"
        )

        # Analysis settings
        config_analysis = c.AnalysisSettings(
            with_all=True,
        )

        # Path to the trajectory file
        config_general = c.GeneralSettings(
            project_name=name,
            export_directory=f"benchmarks/output/bond_strategy/",
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

    except Exception as e:
        print(f"An error occurred while processing {path}: {e}")

def run_coordination_strategy_benchmark(task):
    path, name = task

    try:
        # Lattice settings
        config_lattice = c.LatticeSettings(
            apply_custom_lattice=False,  # If False, read lattice from trajectory file
        )

        # Clustering settings
        config_clustering = c.ClusteringSettings(
            criterion="bond",
            node_types=["1", "2"],
            node_masses=[1.0, 1.0],
            connectivity=["1", "2", "1"],
            cutoffs=[
                c.Cutoff(type1="1", type2="2", distance=0.9)
            ],  # cutoff distance in reduced units
            with_printed_unwrapped_clusters=True,
            print_mode="connectivity", # "all", "connectivity", "individual", "none"
            with_coordination_number=True,
            with_pairwise=True,
            coordination_mode="2",
            coordination_range=[1, 6],
            with_number_of_shared=False,
        )

        # Analysis settings
        config_analysis = c.AnalysisSettings(
            with_all=True,
        )

        # Path to the trajectory file
        config_general = c.GeneralSettings(
            project_name=name,
            export_directory=f"benchmarks/output/coordination_strategy/",
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

    except Exception as e:
        print(f"An error occurred while processing {path}: {e}")

def run_shared_strategy_benchmark(task):
    path, name = task

    try:
        # Lattice settings
        config_lattice = c.LatticeSettings(
            apply_custom_lattice=False,  # If False, read lattice from trajectory file
        )

        # Clustering settings
        config_clustering = c.ClusteringSettings(
            criterion="bond",
            node_types=["1", "2"],
            node_masses=[1.0, 1.0],
            connectivity=["1", "2", "1"],
            cutoffs=[
                c.Cutoff(type1="1", type2="1", distance=1.1),
                c.Cutoff(type1="2", type2="2", distance=1.1),
                c.Cutoff(type1="1", type2="2", distance=0.9)
            ],  # cutoff distance in reduced units
            # with_printed_unwrapped_clusters=True,
            # print_mode="all", # "all", "connectivity", "individual", "none"
            with_coordination_number=True,
            with_pairwise=True,
            coordination_mode="2",
            coordination_range=[1, 6],
            with_number_of_shared=True,
            shared_mode="2",
            shared_threshold=2,
        )

        # Analysis settings
        config_analysis = c.AnalysisSettings(
            with_all=True,
        )

        # Path to the trajectory file
        config_general = c.GeneralSettings(
            project_name=name,
            export_directory=f"benchmarks/output/shared_strategy/",
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

    except Exception as e:
        print(f"An error occurred while processing {path}: {e}")

# --- Main Execution Block ---
if __name__ == "__main__":
    benchmarks = [
        # ('benchmarks/data/universal_benchmark-50000-50.xyz', 'universal_benchmark-50000-50'),
        # ('benchmarks/data/universal_benchmark-75000-50.xyz', 'universal_benchmark-75000-50'),
        # ('benchmarks/data/universal_benchmark-86400-60.xyz', 'universal_benchmark-86400-60'),
        # ('benchmarks/data/universal_benchmark-129600-60.xyz', 'universal_benchmark-129600-60'),
        ('benchmarks/data/universal_benchmark-137200-70.xyz', 'universal_benchmark-137200-70'),
        ('benchmarks/data/universal_benchmark-205800-70.xyz', 'universal_benchmark-205800-70'),
        # ('benchmarks/data/universal_benchmark-204800-80.xyz', 'universal_benchmark-204800-80'),
        # ('benchmarks/data/universal_benchmark-307200-80.xyz', 'universal_benchmark-307200-80'),
        ('benchmarks/data/universal_benchmark-291600-90.xyz', 'universal_benchmark-291600-90'),
        ('benchmarks/data/universal_benchmark-437400-90.xyz', 'universal_benchmark-437400-90'),
        # ('benchmarks/data/universal_benchmark-400000-100.xyz', 'universal_benchmark-400000-100'),
        # ('benchmarks/data/universal_benchmark-600000-100.xyz', 'universal_benchmark-600000-100'),
    ]

    for benchmark in benchmarks:
        run_distance_strategy_benchmark(benchmark)
        run_bond_strategy_benchmark(benchmark)
        run_coordination_strategy_benchmark(benchmark)
        run_shared_strategy_benchmark(benchmark)