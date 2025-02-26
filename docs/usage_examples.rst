Usage Examples
==============

This section provides practical examples of how to use Nexus-CAT for cluster analysis.

### Example 1: Basic Run with SiOz / SiSi / OO Extension 

reference: [Getting Started](getting_started.rst)
```python
import nexus


# Initialize the settings object
extension = "SiOz" # or "SiSi" or "OO"
settings = nexus.settings.Settings(extension)

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

# Run the analysis through the main function using the provided settings object
nexus.main(settings)

print(f"Results are saved here \u279c {settings._output_directory}\n\n") # _output_directory is the combined path of export_directory and project_name

```

### Example 2: Parallel Processing of Multiple Files

reference: [Parallel Processing](https://github.com/JulienPerradin/nexus/blob/f630a7e87da3ea4a1e2028fee54de1002c6fd4d3/scripts/launch-nexus-parallel-multiple-files.py)

```python
import nexus
import os
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor

def process_trajectory(trajectory, output, pressure):
    settings = nexus.settings.Settings(extension="SiOz")
    
    settings.quiet.set_value(True)
    
    # Set output file name
    settings.project_name.set_value(output)
    
    # Set various parameters
    settings.extension.set_value("SiOz")
    settings.export_directory.set_value(f"./export/")
    settings.path_to_xyz_file.set_value(trajectory)
    
    settings.number_of_atoms.set_value(1008)
    settings.header.set_value(2)
    settings.structure.set_value([
                    {"element": "Si", "alias": 2, "number": 336},
                    {"element": "O" , "alias": 1, "number": 672},
                ])
    
    # Set temperature and pressure values (optional)
    settings.temperature.set_value(300) 
    settings.pressure.set_value(pressure)
    
    # Overwrite results if exported files already exist
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

For more examples, please refer to the `scripts` folder in the repository.

