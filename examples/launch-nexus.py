# Import the package
import nexus

# Load trajectory data
trajectory = "tests/inputs/SiO2/1008/sio2-1008at-11frames/pos10.xyz"

# Initialize settings
settings = nexus.settings.Settings(extension="SiOz")

# Set project name
settings.project_name.set_value('pos10')

# Set extension
settings.extension.set_value("SiOz")

# Set export directory
settings.export_directory.set_value(f"tests/results/sio2-1008at-1frame/pos10")

# Set path to XYZ file
settings.path_to_xyz_file.set_value(trajectory)

# Set number of atoms
settings.number_of_atoms.set_value(1008)

# Set header
settings.header.set_value(2)

# Set structure
settings.structure.set_value([
                {"element": "Si", "alias": 2, "number": 336},
                {"element": "O" , "alias": 1, "number": 672},
            ])

# Set temperature
settings.temperature.set_value(300) 

# Set pressure
settings.pressure.set_value(0.00)

# Set to print cluster positions
settings.print_clusters_positions.set_value(True)

# Set cluster settings
settings.cluster_settings.set_cluster_parameter("connectivity", ["Si", "O", "Si"])
settings.cluster_settings.set_cluster_parameter("criteria", "bond")
# settings.cluster_settings.set_cluster_parameter("polyhedra", [[4, 4], [5, 5], [6, 6]])
settings.cluster_settings.set_cluster_parameter("polyhedra", [])

# Execute the main function with the provided settings
nexus.main(settings)
