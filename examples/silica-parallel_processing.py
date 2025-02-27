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

