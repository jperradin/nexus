.. _getting_started:

Getting Started
===============

This section will guide you through the initial setup and usage of Nexus-CAT.

### Quick Start Example

1. Prepare a trajectory file in the `Extended XYZ` format, or use the example file provided at `examples/inputs/SiO2-27216at-pos67B.xyz`.
    
    The trajectory file must have the following format:

    ```plaintext
    1008 # Number of atoms
    Lattice="21.3850 0.0 0.0 0.0 21.3850 0.0 0.0 0.0 21.3850" # Lattice informations
    O 4.81534 7.79666 13.4134       # Element x y z
    O 15.3218 18.7985 16.3627
    ...
    Si 16.9122 10.503 3.05908
    Si 5.69312 8.20448 18.6908
    ...
    ```

    Nexus-CAT reads lattice properties directly from the trajectory file.

2. Create a python script and import the necessary modules

    ```python
    # Import necessary modules
    from nexus import SettingsBuilder, main
    import nexus.config.settings as c
    ```

3. Set up the general configuration settings

    ```python
    # General settings
    config_general = c.GeneralSettings(
        project_name="SiO2",                    # Project name
        export_directory="examples/exports",    # Export directory
        file_location="path/to/trajectory.xyz", # File location
        range_of_frames=(0, 10),                # Range of frames to analyze
        apply_pbc=True,                         # Apply periodic boundary conditions
        verbose=True,                           # Verbose mode (print progress bars, etc.)
        save_logs=True,                         # Save logs to export directory
        save_performance=True                   # Track & save performance metrics
    )
    ```

4. Configure lattice settings

    ```python
    # Lattice settings
    config_lattice = c.LatticeSettings(
        apply_custom_lattice=False,             # If False, read lattice from trajectory file
        # Custom lattice can be specified as a numpy array if needed
    )
    ```

5. Set up the clustering settings
   
    ```python
    # Clustering settings
    config_clustering = c.ClusteringSettings(
        criteria="bond",                        # "bond" or "distance" criteria
        node_types=["Si", "O"],                 # Types of nodes in the system
        node_masses=[28.0855, 15.9994],        # Masses of nodes in reduced units
        connectivity=["Si", "O", "Si"],         # Connectivity pattern to analyze
        cutoffs=[c.Cutoff(type1="Si", type2="Si", distance=3.50),  # Cutoff distances
                 c.Cutoff(type1="Si", type2="O", distance=2.30),   # between different
                 c.Cutoff(type1="O", type2="O", distance=3.05)],   # types of nodes
        
        with_coordination_number=True,          # Calculate coordination numbers
        coordination_mode="O",                  # Mode for coordination calculations
        coordination_range=[4, 6],              # Range of coordination numbers to consider
        
        with_alternating=True,                  # Calculate alternating coordinations
        
        with_printed_unwrapped_clusters=False,  # Whether to output cluster coordinates
        print_mode="connectivity"               # Mode for printing clusters
    )
    ```

6. Configure the analysis settings

    ```python
    # Analysis settings
    config_analysis = c.AnalysisSettings(
        with_all=True,                          # Enable all analysis methods
        # Specific analyses can be enabled individually as needed
    )
    ```

7. Build the settings object using the builder pattern and run the analysis

    ```python
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

8. Run the script and check the results in the export directory.

    ```bash
    python your_script.py
    ```

9. View the results in the specified export directory.

    ```bash
    ls examples/exports/SiO2/
    ```

### Advanced Usage: Reconfiguring and Rerunning

You can modify settings and rerun the analysis with new parameters:

```python
# Update clustering settings
config_clustering.with_number_of_shared = True
config_clustering.coordination_range = [6, 6]  # Looking for SiO6-SiO6 connections
config_analysis.overwrite = False              # Preserve previous results

# Rebuild settings with updated configuration
settings = (SettingsBuilder()
    .with_general(config_general)
    .with_lattice(config_lattice)
    .with_clustering(config_clustering)
    .with_analysis(config_analysis)
    .build()
)

# Run the analysis again with new settings
main(settings)
```

For different and more advanced examples, see the [usage_examples](usage_examples.rst) section.
