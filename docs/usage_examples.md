Usage Examples
==============

This section provides practical examples of how to use Nexus-CAT for cluster analysis.

### Example 1: Water Ice Analysis

```python
# Import the package
import nexus

# Load trajectory data
# /!\ Only extended XYZ files are supported
#     No need to provide the lattice vectors, the code will automatically detect them
#     with the keyword 'Lattice' in the XYZ file.
#     The code will crash if the file is not an extended XYZ file
trajectory = "./examples/inputs/waterM2825-5kbar.xyz"

# Initialize settings
settings = nexus.settings.Settings(extension="OO")
settings.quiet.set_value(False)

# Set project name, this will be used to name the output directory in the export directory
settings.project_name.set_value("water_ice")

# Set extension
settings.extension.set_value("OO")

# Set export directory
settings.export_directory.set_value(f"./examples/export/")

# Set path to XYZ file
settings.path_to_xyz_file.set_value(trajectory)

# Set number of atoms
# /!\ This value must be the same as the number of atoms in the XYZ file.
#     If the number of atoms is not provided, or the value provided is wrong, the code will crash.
settings.number_of_atoms.set_value(2825)

# Set range of frames (optional)
# settings.range_of_frames.set_value([2, 5]) # Only frames 2 to 5 will be processed

# Set header of the XYZ file
#   (ie, number of atoms in the first line, lattice properties in the second line)
settings.header.set_value(2)

# Set structure
# /!\ This value must be the same as the number of atoms in the XYZ file.
#     If the number of atoms is not provided, or the value provided is wrong, the code will crash.
settings.structure.set_value(
    [
        {"element": "Si", "number": 2825},
    ]
)

# Set temperature in Kelvin (optional but recommended for the recap. of the results)
settings.temperature.set_value(0)

# Set pressure in GPa (optional but recommended for the recap. of the results)
settings.pressure.set_value(0.05)

# Set to print cluster positions (optional, default is False, if set to True, the user will be prompted to confirm the action)
# uncomment the following to remove the warning
settings.print_clusters_positions.disable_warnings = True
settings.print_clusters_positions.set_value(True)

# Set to not overwrite results to compare with the previous results (optional, default is True)
settings.overwrite_results.set_value(True)

# Set cluster analysis criteria (bond or distance)
settings.cluster_settings.set_cluster_parameter("criteria", "distance")
# Set cluster connectivities to look for
settings.cluster_settings.set_cluster_parameter("connectivity", ["O", "O"])
# Set polyhedra to look for
# settings.cluster_settings.set_cluster_parameter(
#     'polyhedra', [[4, 4], [5, 5], [6, 6]]
# )

# Run the main function with the provided settings
print("Processing the trajectory with 'distance' criteria ...")
nexus.main(settings)

# Print the path to the results
print("\n\n\t\tAll trajectories have been processed successfully.")
print(f"\n\t\tResults are saved here \u279c {settings._output_directory}\n\n")
```

### Example 2: Silica Sequential Processing

```python
# Import the package
import nexus

# Load trajectory data
# /!\ Only extended XYZ files are supported
#     No need to provide the lattice vectors, the code will automatically detect them
#     with the keyword 'Lattice' in the XYZ file.
#     The code will crash if the file is not an extended XYZ file
trajectory = "./examples/inputs/SiO2-27216at-pos67B.xyz"

# Initialize settings
settings = nexus.settings.Settings(extension="SiOz")
settings.quiet.set_value(True)

# Set project name, this will be used to name the output directory in the export directory
settings.project_name.set_value("silica-sequential_processing")

# Set extension
settings.extension.set_value("SiOz")

# Set export directory
settings.export_directory.set_value(f"./examples/export/")

# Set path to XYZ file
settings.path_to_xyz_file.set_value(trajectory)

# Set number of atoms
# /!\ This value must be the same as the number of atoms in the XYZ file.
#     If the number of atoms is not provided, or the value provided is wrong, the code will crash.
settings.number_of_atoms.set_value(27216)

# Set range of frames (optional)
# settings.range_of_frames.set_value([2, 5]) # Only frames 2 to 5 will be processed

# Set header of the XYZ file
#   (ie, number of atoms in the first line, lattice properties in the second line)
settings.header.set_value(2)
settings.range_of_frames.set_value([0, 1])

nSi = int(27216 / 3)
nO = int(nSi * 2)

# Set structure
# /!\ This value must be the same as the number of atoms in the XYZ file.
#     If the number of atoms is not provided, or the value provided is wrong, the code will crash.
settings.structure.set_value(
    [
        {"element": "Si", "number": nSi},
        {"element": "O", "number": nO},
    ]
)

# Set temperature in Kelvin (optional but recommended for the recap. of the results)
settings.temperature.set_value(300)

# Set pressure in GPa (optional but recommended for the recap. of the results)
settings.pressure.set_value(10.0)

# Set to print cluster positions (optional, default is False, if set to True, the user will be prompted to confirm the action)
# uncomment the following to remove the warning
settings.print_clusters_positions.disable_warnings = True
settings.print_clusters_positions.set_value(True)

# Set to not overwrite results to compare with the previous results (optional, default is True)
settings.overwrite_results.set_value(True)

# Set cluster analysis criteria (bond or distance)
settings.cluster_settings.set_cluster_parameter("criteria", "bond")
# Set cluster connectivities to look for
settings.cluster_settings.set_cluster_parameter("connectivity", ["Si", "O", "Si"])
# Set polyhedra to look for
settings.cluster_settings.set_cluster_parameter("polyhedra", [[4, 4], [5, 5], [6, 6]])

# Run the main function with the provided settings
print("Processing the trajectory with 'bond' criteria ...")
nexus.main(settings)

# Print the path to the results
print("\n\n\t\tAll trajectories have been processed successfully.")
print(f"\n\t\tResults are saved here \u279c {settings._output_directory}\n\n")
```

### Example 3: Silica Parallel Processing

```python
import nexus
import os
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor

def process_trajectory(trajectory, output, pressure):
    settings = nexus.settings.Settings(extension="SiOz")
    
    settings.quiet.set_value(True)
    
    # Set project name, this will be used to name the output directory in the export directory
    settings.project_name.set_value(output)
    
    # Set extension
    settings.extension.set_value("SiOz")
    
    # Set export directory
    settings.export_directory.set_value(f"./examples/export/silica-parallel_processing")
    
    # Set path to XYZ file
    settings.path_to_xyz_file.set_value(trajectory)
    
    # Set number of atoms
    # /!\ This value must be the same as the number of atoms in the XYZ file.
    #     If the number of atoms is not provided, or the value provided is wrong, the code will crash.
    settings.number_of_atoms.set_value(1008)
    
    # Set header of the XYZ file
    #   (ie, number of atoms in the first line, lattice properties in the second line)
    settings.header.set_value(2)
    
    # Set structure
    # /!\ This value must be the same as the number of atoms in the XYZ file.
    #     If the number of atoms is not provided, or the value provided is wrong, the code will crash.
    settings.structure.set_value([
                    {"element": "Si", "alias": 2, "number": 336},
                    {"element": "O" , "alias": 1, "number": 672},
                ])
    
    # Set temperature in Kelvin (optional but recommended for the recap. of the results)
    settings.temperature.set_value(300) 
    
    # Set pressure in GPa (optional but recommended for the recap. of the results)
    settings.pressure.set_value(pressure)
    
    # Set to not overwrite results to compare with the previous results (optional, default is True)
    settings.overwrite_results.set_value(True)
    
    # Set cluster parameter SiOz polyhedra + 'bond' criteria
    settings.cluster_settings.set_cluster_parameter("connectivity", ["Si", "O", "Si"])
    settings.cluster_settings.set_cluster_parameter("criteria", "bond")
    settings.cluster_settings.set_cluster_parameter("polyhedra", [[4, 4], [4, 5], [5, 5], [5, 6], [6, 6]])
    
    # Run the main function for SiOz polyhedra clusters
    nexus.main(settings)
    
    # Set cluster parameter OSiz units + 'bond' criteria
    settings.cluster_settings.set_cluster_parameter("connectivity", ["O", "Si", "O"])
    settings.cluster_settings.set_cluster_parameter("criteria", "bond")
    settings.cluster_settings.set_cluster_parameter("polyhedra", [[2, 2], [2, 3], [3, 3]])
    
    # Do not overwrite results if exported files already exist
    settings.overwrite_results.set_value(False)

    # Run the main function
    nexus.main(settings)
    
    # Set cluster parameter OO units + 'distance' criteria
    settings.cluster_settings.set_cluster_parameter("connectivity", ["O", "O"])
    settings.cluster_settings.set_cluster_parameter("criteria", "distance")
    settings.cluster_settings.set_cluster_parameter("polyhedra", [[2, 2], [2, 3], [3, 3]])
    
    # Run the main function
    nexus.main(settings)
    
    # Set cluster parameter SiSi units + 'distance' criteria
    settings.cluster_settings.set_cluster_parameter("connectivity", ["Si", "Si"])
    settings.cluster_settings.set_cluster_parameter("criteria", "distance")
    settings.cluster_settings.set_cluster_parameter("polyhedra", [[4, 4], [4, 5], [5, 5], [5, 6], [6, 6]])

    # Run the main function
    nexus.main(settings)

# Load trajectories, output names, and pressures from files
directory = "path/to/directory"
trajectories = []
pressures = []
outputs = []

with open(os.path.join(directory, "inputs")) as f:
    data = f.readlines()
    for i, l in enumerate(data):
        trajectories.append(l.strip())
f.close()

with open(os.path.join(directory, "outputs")) as f:
    data = f.readlines()
    for i, l in enumerate(data):
        outputs.append(l.strip())
f.close()

with open(os.path.join(directory, "pressure")) as f:
    data = f.readlines()
    for i, l in enumerate(data):
        pressures.append(float(l.strip()))
f.close()

# Initialize a progress bar
progress_bar = tqdm(enumerate(trajectories), total=len(trajectories), desc="", colour="#510e4c", unit='file', leave=False)

# Fancy color bar
color_gradient = nexus.utils.generate_color_gradient(len(trajectories))

# Initialize settings
settings = nexus.settings.Settings(extension="SiOz")

# Enable print clusters positions
settings.print_clusters_positions.disable_warnings = True # Disable warnings if false, it ask to user to confirm the action (every loop)
settings.print_clusters_positions.set_value(False) # Set to True to print cluster coordinates to XYZ files (/!\ careful for storage consumption)

with ProcessPoolExecutor(max_workers=8) as executor:
    # Parallel processing of the trajectories
    futures = []
    for i, trajectory in progress_bar:
        output = outputs[i]
        pressure = pressures[i]
        progress_bar.set_description(f"Processing ... \u279c {str(trajectory).split('/')[-1]}")
        progress_bar.colour = "#%02x%02x%02x" % color_gradient[i]
        future = executor.submit(process_trajectory, trajectory, output, pressure)
        futures.append(future)
    
    for future in tqdm(futures, total=len(futures), desc="Waiting for the results", colour="#510e4c", unit='file', leave=False):
        future.result()
        
# Print success message (if it appears, everything went well)
print("\n\n\t\tAll trajectories have been processed successfully.")
print(f"\n\t\tResults are saved here \u279c {settings.export_directory.get_value()}\n\n")
```

### Example 4: Amorphous Silicon Analysis

```python
# Import the package
import nexus

# Load trajectory data
# /!\ Only extended XYZ files are supported
#     No need to provide the lattice vectors, the code will automatically detect them
#     with the keyword 'Lattice' in the XYZ file.
#     The code will crash if the file is not an extended XYZ file
trajectory = "./examples/inputs/aSi-compress_GAP-18_10GPa.xyz"

# Initialize settings
settings = nexus.settings.Settings(extension="SiSi")
settings.quiet.set_value(False)

# Set project name, this will be used to name the output directory in the export directory
settings.project_name.set_value("amorphous_silicon")

# Set extension
settings.extension.set_value("SiSi")

# Set export directory
settings.export_directory.set_value(f"./examples/export/")

# Set path to XYZ file
settings.path_to_xyz_file.set_value(trajectory)

# Set number of atoms
# /!\ This value must be the same as the number of atoms in the XYZ file.
#     If the number of atoms is not provided, or the value provided is wrong, the code will crash.
settings.number_of_atoms.set_value(100000)

# Set range of frames (optional)
# settings.range_of_frames.set_value([2, 5]) # Only frames 2 to 5 will be processed

# Set header of the XYZ file
#   (ie, number of atoms in the first line, lattice properties in the second line)
settings.header.set_value(2)

# Set structure
# /!\ This value must be the same as the number of atoms in the XYZ file.
#     If the number of atoms is not provided, or the value provided is wrong, the code will crash.
settings.structure.set_value(
    [
        {"element": "Si", "number": 100000},
    ]
)

# Set temperature in Kelvin (optional but recommended for the recap. of the results)
settings.temperature.set_value(500)

# Set pressure in GPa (optional but recommended for the recap. of the results)
settings.pressure.set_value(10.0)

# Set to print cluster positions (optional, default is False, if set to True, the user will be prompted to confirm the action)
# uncomment the following to remove the warning
settings.print_clusters_positions.disable_warnings = True
settings.print_clusters_positions.set_value(True)

# Set to not overwrite results to compare with the previous results (optional, default is True)
settings.overwrite_results.set_value(True)

# Set cluster analysis criteria (bond or distance)
settings.cluster_settings.set_cluster_parameter("criteria", "distance")
# Set cluster connectivities to look for
settings.cluster_settings.set_cluster_parameter("connectivity", ["Si", "Si"])
# Set polyhedra to look for
settings.cluster_settings.set_cluster_parameter(
    'polyhedra', [[4, 4], [5, 5], [6, 6]]
)

# Run the main function with the provided settings
print("Processing the trajectory with 'distance' criteria ...")
nexus.main(settings)

# Print the path to the results
print("\n\n\t\tAll trajectories have been processed successfully.")
print(f"\n\t\tResults are saved here \u279c {settings._output_directory}\n\n")
```

For more examples, please refer to the `examples` folder in the repository.

