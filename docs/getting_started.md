.. _getting_started:

Getting Started
===============

This section will guide you through the initial setup and usage of Nexus-CAT.

### Quick Start Example

1. Prepare a trajectory file in the `Extended XYZ` format.
    
    For now, the trajectory file must have the following format:

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

    (Next versions of Nexus-CAT will handle trajectory file format automatically)

2. Create a python script and import nexus 

    ```python
    import nexus
    ```

3. Initialize the settings using the extention `SiOz`.

    The core of the program of Nexus-CAT works alongside extensions. 
    They contain the necessary settings and functions to analyze the trajectory file.
    The `SiOz` extension is used to analyze SiO2 systems with z being the coordination number of Si atoms.
    (Refer to the [extensions](extensions.rst) section for more information).

    ```python
    # Initialize the settings object
    extension = "SiOz"
    settings = nexus.settings.Settings(extension)
    ```

4. Set up the input/output and structure settings
    ```python
    # Set the trajectory file
    settings.project_name.set_value("quick_start") 
    settings.export_directory.set_value("./export")

    # Load the trajectory file
    trajectory = "./path/to/trajectory.xyz"         # path of trajectory file
    settings.path_to_xyz_file.set_value(trajectory) # set the path to the trajectory file
    settings.number_of_atoms.set_value(100)         # set the number of atoms in the trajectory file
    settings.header.set_value(2)                    # set the number of header lines in the trajectory file
    settings.range_of_frames.set_value([0, 1])      # set the range of frames to be analyzed \[start, end\]

    # Set structure informations
    settings.structure.set_value([
        {"element": "Si", "number": 336},
        {"element": "O", "number": 672},     # each element can be added separately
    ])
    ```

5. Set up the cluster settings

    In this example, we will set the cluster analysis criteria to `bond`, the connectivities to look for to `['Si', 'O', 'Si']`, and the polyhedra to look for to `[[4, 4], [5, 5], [6, 6]]`.
    Which means we are looking for SiO$_4$-SiO$_4$, SiO$_5$-SiO$_5$, and SiO$_6$-SiO$_6$ clusters where the connectivities between SiO$_z$ polyhedra are bridges formed by the Si-O-Si bonds.

    ```python
    # Set cluster analysis criteria (bond or distance)
    settings.cluster_settings.set_cluster_parameter(
        'criteria', 'bond'
    )
    # Set cluster connectivities to look for
    settings.cluster_settings.set_cluster_parameter(
        'connectivity', ['Si', 'O', 'Si']
    )
    # Set polyhedra to look for
    settings.cluster_settings.set_cluster_parameter(
        'polyhedra', [[4, 4], [5, 5], [6, 6]]
    )
    ```

6. Run the analysis through the main function using the provided settings object

    ```python
    nexus.main(settings)

    # _output_directory is the combined path of export_directory and project_name
    print(f"Results are saved here \u279c {settings._output_directory}\n\n") 
    ```

7. Run the script and check the results in the export directory.

    ```bash
    python quick_start.py
    ```

8. View the results in the specified export directory.

    ```bash
    ls export/quick_start/
    ```

For different and more advanced examples, see the [usage_examples](usage_examples.rst) section.
