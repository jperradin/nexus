import nexus
import os
from tqdm import tqdm

# Load trajectories, output names, and pressures from files
directory = "tests/inputs/SiO2/1008/sio2-1008at-11frames/"
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
settings.print_clusters_positions.set_value(False)

# Loop over each trajectory
for i, trajectory in progress_bar:
    
    # Update progress bar description
    progress_bar.set_description(f"Processing ... \u279c {str(trajectory).split('/')[-1]}")
    progress_bar.colour = "#%02x%02x%02x" % color_gradient[i]
    
    # Set output file name
    output = outputs[i]

    # Adjust verbosity for the first iteration
    if i == 0:
        settings.quiet.set_value(False)
    else:
        settings.quiet.set_value(True)

    # Set various parameters
    settings.project_name.set_value(output)
    settings.extension.set_value("SiOz")
    settings.export_directory.set_value(f"tests/results/sio2-1008at-11frames-all")
    settings.path_to_xyz_file.set_value(trajectory)
    
    settings.number_of_atoms.set_value(1008)
    settings.header.set_value(2)
    settings.structure.set_value([
                    {"element": "Si", "alias": 2, "number": 336},
                    {"element": "O" , "alias": 1, "number": 672},
                ])
    
    settings.temperature.set_value(300) 
    settings.pressure.set_value(pressures[i])
    
    settings.overwrite_results.set_value(True)
    
    # Set cluster parameter 'bond' criteria
    settings.cluster_settings.set_cluster_parameter("connectivity", ["Si", "O", "Si"])
    settings.cluster_settings.set_cluster_parameter("criteria", "bond")
    settings.cluster_settings.set_cluster_parameter("polyhedra", [[4, 4], [4, 5], [5, 5], [5, 6], [6, 6]])
    
    # Run the main function
    nexus.main(settings)
    
    # Set to quiet mode after first iteration of settings
    settings.quiet.set_value(True)
    settings.overwrite_results.set_value(False)
    
    # Set cluster parameter 'distance' criteria
    settings.cluster_settings.set_cluster_parameter("connectivity", ["O", "O"])
    settings.cluster_settings.set_cluster_parameter("criteria", "distance")
    settings.cluster_settings.set_cluster_parameter("polyhedra", [[2, 2], [2, 3], [3, 3]])
    
    # Run the main function
    nexus.main(settings)
    
    # Set cluster parameter 'distance' criteria
    settings.cluster_settings.set_cluster_parameter("connectivity", ["Si", "Si"])
    settings.cluster_settings.set_cluster_parameter("criteria", "distance")
    settings.cluster_settings.set_cluster_parameter("polyhedra", [[4, 4], [4, 5], [5, 5], [5, 6], [6, 6]])

    # Run the main function
    nexus.main(settings)

print("\n\n\t\tAll trajectories have been processed successfully.")
print(f"\n\t\tResults are saved here \u279c {settings.export_directory.get_value()}\n\n")
    
