.. _getting_started:

Getting Started
===============

This section will guide you through the initial setup and usage of Nexus-CAT.

### Quick Start Example

1.  **Prepare a trajectory file** in the `Extended XYZ` format. The trajectory file must have the following format:

    ```plaintext
    1008
    Lattice="21.3850 0.0 0.0 0.0 21.3850 0.0 0.0 0.0 21.3850"
    O 4.81534 7.79666 13.4134
    O 15.3218 18.7985 16.3627
    ...
    Si 16.9122 10.503 3.05908
    Si 5.69312 8.20448 18.6908
    ...
    ```
    Nexus-CAT reads lattice properties directly from the trajectory file.

2.  **Create a python script** and import the necessary modules:

    ```python
    from nexus import SettingsBuilder, main
    import nexus.config.settings as c
    ```

3.  **Set up the general configuration settings**:

    ```python
    config_general = c.GeneralSettings(
        project_name="SiO2",
        export_directory="examples/exports",
        file_location="path/to/trajectory.xyz",
        range_of_frames=(0, 10),
        apply_pbc=True,
        verbose=True,
        save_logs=True,
        save_performance=True
    )
    ```

4.  **Configure lattice settings**:

    ```python
    config_lattice = c.LatticeSettings(
        apply_custom_lattice=False,
    )
    ```

5.  **Set up the clustering settings**:

    ```python
    config_clustering = c.ClusteringSettings(
        criterion="bond",
        node_types=["Si", "O"],
        node_masses=[28.0855, 15.9994],
        connectivity=["Si", "O", "Si"],
        cutoffs=[c.Cutoff(type1="Si", type2="Si", distance=3.50),
                 c.Cutoff(type1="Si", type2="O", distance=2.30),
                 c.Cutoff(type1="O", type2="O", distance=3.05)],
        with_coordination_number=True,
        coordination_mode="O",
        coordination_range=[4, 6],
        with_alternating=True,
        with_printed_unwrapped_clusters=False,
        print_mode="connectivity"
    )
    ```

6.  **Configure the analysis settings**:

    ```python
    config_analysis = c.AnalysisSettings(
        with_all=True,
    )
    ```

7.  **Build the settings object** using the builder pattern and run the analysis:

    ```python
    settings = (SettingsBuilder()
        .with_general(config_general)
        .with_lattice(config_lattice)
        .with_clustering(config_clustering)
        .with_analysis(config_analysis)
        .build()
    )

    main(settings)
    ```

8.  **Run the script** and check the results in the export directory.

    ```bash
    python your_script.py
    ```