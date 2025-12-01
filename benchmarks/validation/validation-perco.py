# Import necessary modules
from nexus import SettingsBuilder, main
import nexus.config.settings as c
import nexus.io.parser.parser as p
from multiprocessing import Pool
from tqdm import tqdm
import numpy as np


def run_single_analysis(task):
    path, name, j = task

    try:
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

        # Path to the trajectory file
        config_general = c.GeneralSettings(
            project_name=name,
            export_directory=f"benchmarks/validation/output/{j}",
            file_location=path,
            range_of_frames=(0, -1),
            apply_pbc=True,
            verbose=False,
            save_logs=False,
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
    # Load file paths and output names from the input files
    try:
        files, outputs, j = np.loadtxt(
            "benchmarks/validation/nexus_inputs_critical_point",
            dtype="<U100",
            unpack=True,
            comments="#",
        )  # Increased size for longer paths
    except IOError as e:
        print(f"Error reading input files: {e}")
        print(
            "Please ensure 'inputs' and 'outputs' files exist and are correctly formatted."
        )
        exit()

    # Create a list of tasks to be processed
    if files.shape != outputs.shape:
        print(
            "Error: The 'inputs' and 'outputs' files have a different number of lines."
        )
        exit()

    tasks = list(zip(files, outputs, j))

    if not tasks:
        print("No tasks found in input files.")
    else:
        # Define the number of parallel workers
        num_workers = 4
        print(
            f"Starting parallel processing for {len(tasks)} files with {num_workers} workers..."
        )

        # Create a multiprocessing Pool
        with Pool(processes=num_workers) as pool:
            # Use tqdm to create a progress bar for the parallel execution
            # pool.imap_unordered processes tasks as they are submitted and returns results as they complete

            list(
                tqdm(
                    pool.imap_unordered(run_single_analysis, tasks),
                    total=len(tasks),
                    desc="Analyzing Files",
                )
            )

        print("\n--- All analyses complete ---")
