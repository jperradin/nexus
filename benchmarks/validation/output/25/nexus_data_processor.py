import numpy as np
import os
from natsort import natsorted
import sys


class DataPoint:
    """Represents a single data point from the new CSV format"""

    def __init__(
        self, dens_dir, file_name, connectivity_type, concentration, value, error
    ):
        self.dens_dir = dens_dir
        self.file_name = file_name
        self.connectivity_type = connectivity_type
        self.concentration = concentration
        self.value = value
        self.error = error
        self.pressure = 0.0
        self.temperature = 0.0
        self.box_size = 0.0

    def set_conditions(self, pressure, temperature, box_size):
        self.pressure = pressure
        self.temperature = temperature
        self.box_size = box_size

    def get_density_value(self):
        """Extract numerical density from directory name like 'dens2.240'"""
        return float(self.dens_dir.split("p")[1])


class DataCollection:
    """Manages collection of data points and organizes them for export"""

    def __init__(self):
        self.data_points = []

    def add_data_point(self, data_point):
        self.data_points.append(data_point)

    def get_unique_files(self):
        return list(set(dp.file_name for dp in self.data_points))

    def get_connectivity_types_for_file(self, file_name):
        """Get all connectivity types present in a specific file"""
        types = set(
            dp.connectivity_type for dp in self.data_points if dp.file_name == file_name
        )
        return sorted(list(types))

    def get_data_by_variable(self, file_name, connectivity_type, variable):
        """Extract data for a specific file, connectivity type, and variable (pressure, density, etc.)"""
        relevant_points = [
            dp
            for dp in self.data_points
            if dp.file_name == file_name and dp.connectivity_type == connectivity_type
        ]

        data = []
        for dp in relevant_points:
            if variable == "pressure":
                x_val = dp.pressure
            elif variable == "temperature":
                x_val = dp.temperature
            elif variable == "density":
                x_val = dp.get_density_value()
            elif variable == "box":
                x_val = dp.box_size
            elif variable == "concentration":
                x_val = dp.concentration
            else:
                continue

            data.append((x_val, dp.value, dp.error))

        # Sort by x variable
        data.sort(key=lambda x: x[0])
        return data


class DataExporter:
    """Handles exporting organized data to files with proper headers"""

    def __init__(self, data_collection, output_dir="./export"):
        self.data_collection = data_collection
        self.output_dir = output_dir
        self._ensure_output_directories()

    def _ensure_output_directories(self):
        """Create output directories if they don't exist"""
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)

        concentrations_dir = os.path.join(self.output_dir, "concentrations")
        if not os.path.exists(concentrations_dir):
            os.makedirs(concentrations_dir)

    def export_by_variable(self, variable_name):
        """Export data organized by a specific variable (pressure, density, etc.)"""

        for file_name in self.data_collection.get_unique_files():
            connectivity_types = self.data_collection.get_connectivity_types_for_file(
                file_name
            )

            if not connectivity_types:
                continue

            # Collect data for all connectivity types
            all_data = []
            for conn_type in connectivity_types:
                data = self.data_collection.get_data_by_variable(
                    file_name, conn_type, variable_name
                )
                all_data.append(data)

            if not all_data or not all_data[0]:
                continue

            # Create output filename
            output_filename = f"{variable_name}-{file_name}"
            output_path = os.path.join(self.output_dir, output_filename)

            # Export main data file
            self._write_data_file(
                output_path, variable_name, connectivity_types, all_data
            )

            # Export error file
            error_path = output_path.replace(".dat", "--errors.dat")
            self._write_error_file(
                error_path, variable_name, connectivity_types, all_data
            )

            # For concentration data, also create individual files
            if variable_name == "concentration":
                self._export_individual_concentration_files(
                    file_name, connectivity_types, all_data
                )

    def _write_data_file(
        self, output_path, variable_name, connectivity_types, all_data
    ):
        """Write main data file with indexed headers"""

        # Find the maximum number of data points across all connectivity types
        max_points = (
            max(len(data_series) for data_series in all_data) if all_data else 0
        )

        if max_points == 0:
            print(f"Warning: No data to write for {output_path}")
            return

        n_columns = len(connectivity_types) + 2  # +1 for variable column
        data_matrix = np.zeros((n_columns, max_points))

        # Collect all unique x-values (variable values) and sort them
        all_x_values = set()
        for data_series in all_data:
            i = 0
            for x, y, err in data_series:
                data_matrix[0, i] = x  # Variable column
                data_matrix[1, i] = y # Initialize data columns
                data_matrix[2, i] = err # Initialize error columns
                all_x_values.add(x)
                i += 1

        sorted_x_values = sorted(list(all_x_values))

        # # Create a mapping from x-value to index
        # x_to_index = {x_val: idx for idx, x_val in enumerate(sorted_x_values)}

        # # Fill the variable column
        # for i, x_val in enumerate(sorted_x_values):
        #     data_matrix[0, i] = x_val

        # # Fill data for each connectivity type
        # for conn_idx, data_series in enumerate(all_data):
        #     # Create a dictionary for quick lookup of y-values by x-value
        #     y_dict = {x: y for x, y, err in data_series}

        #     # Fill the column, using 0.0 for missing data points
        #     for i, x_val in enumerate(sorted_x_values):
        #         if x_val in y_dict:
        #             data_matrix[conn_idx + 1, i] = y_dict[x_val]
        #         else:
        #             data_matrix[conn_idx + 1, i] = 0.0

        # Write file
        with open(output_path, "w") as f:
            # Write headers with column indices
            f.write(f"# 1 {variable_name}\n")
            for i, conn_type in enumerate(connectivity_types):
                f.write(f"# {i + 2} {conn_type}\n")

            # Write data (only non-empty rows)
            for point_idx in range(len(sorted_x_values)):
                for col_idx in range(n_columns):
                    f.write(f"{data_matrix[col_idx, point_idx]:^10.5f}\t")
                f.write("\n")

        print(f"Exported: {output_path}")

    def _write_error_file(
        self, output_path, variable_name, connectivity_types, all_data
    ):
        """Write error file with indexed headers"""

        # Find the maximum number of data points and collect all x-values
        max_points = (
            max(len(data_series) for data_series in all_data) if all_data else 0
        )

        if max_points == 0:
            return

        # Collect all unique x-values and sort them
        all_x_values = set()
        for data_series in all_data:
            for x, y, err in data_series:
                all_x_values.add(x)

        sorted_x_values = sorted(list(all_x_values))
        n_columns = len(connectivity_types) + 1

        error_matrix = np.zeros((n_columns, len(sorted_x_values)))

        # Fill the variable column
        for i, x_val in enumerate(sorted_x_values):
            error_matrix[0, i] = x_val

        # Fill error data for each connectivity type
        for conn_idx, data_series in enumerate(all_data):
            # Create a dictionary for quick lookup of error values by x-value
            err_dict = {x: err for x, y, err in data_series}

            # Fill the column, using 0.0 for missing data points
            for i, x_val in enumerate(sorted_x_values):
                if x_val in err_dict:
                    error_matrix[conn_idx + 1, i] = err_dict[x_val]
                else:
                    error_matrix[conn_idx + 1, i] = 0.0

        # Write file
        with open(output_path, "w") as f:
            # Write headers with column indices
            f.write(f"# 1 {variable_name}\n")
            for i, conn_type in enumerate(connectivity_types):
                f.write(f"# {i + 2} {conn_type}\n")

            # Write data
            for point_idx in range(len(sorted_x_values)):
                for col_idx in range(n_columns):
                    f.write(f"{error_matrix[col_idx, point_idx]:^10.5f}\t")
                f.write("\n")

    def export_concentration_files(self):
        """Export individual concentration files for each connectivity type"""

        concentrations_dir = os.path.join(self.output_dir, "concentrations")
        if not os.path.exists(concentrations_dir):
            os.makedirs(concentrations_dir)

        for file_name in self.data_collection.get_unique_files():
            connectivity_types = self.data_collection.get_connectivity_types_for_file(
                file_name
            )

            # Get the property name from the file (remove .dat extension)
            property_name = file_name.replace(".dat", "")

            for conn_type in connectivity_types:
                # Get concentration data for this specific connectivity type
                concentration_data = self.data_collection.get_data_by_variable(
                    file_name, conn_type, "concentration"
                )

                if not concentration_data:
                    continue

                # Sort by concentration
                concentration_data.sort(key=lambda x: x[0])

                # Create filename: concentration-{property}-{connectivity}.dat
                output_filename = f"concentration-{property_name}-{conn_type}.dat"
                output_path = os.path.join(concentrations_dir, output_filename)

                # Write the file
                with open(output_path, "w") as f:
                    f.write("# 1 concentration\n")
                    f.write("# 2 data\n")
                    f.write("# 3 error\n")

                    for concentration, value, error in concentration_data:
                        f.write(
                            f"{concentration:^10.5f}\t{value:^10.5f}\t{error:^10.5f}\n"
                        )

                print(f"Exported concentration file: {output_path}")


def parse_csv_file(filepath):
    """Parse the new CSV-style format files"""
    results = []

    try:
        with open(filepath, "r") as f:
            lines = f.readlines()

        # Find data lines (not comments, contains commas)
        for line in lines:
            line = line.strip()
            if line and not line.startswith("#") and "," in line:
                parts = line.split(",")
                if len(parts) >= 4:
                    connectivity_type = parts[0].strip()
                    concentration = float(parts[1].strip())
                    value = float(parts[2].strip())
                    std = parts[3].strip()
                    error_str = parts[4].strip()
                    error = 0.0 if error_str == "nan" else float(error_str)

                    results.append(
                        {
                            "connectivity_type": connectivity_type,
                            "concentration": concentration,
                            "value": value,
                            "std": std,
                            "error": error,
                        }
                    )

    except Exception as e:
        print(f"Error parsing {filepath}: {e}")
        return []

    return results


def load_condition_files(list_dens, unloading=False):
    """Load pressure, temperature, and box size data - skip missing files"""

    conditions = {"pressures": {}, "temperatures": {}, "boxes": {}}

    # Load data files - only process files that exist
    condition_files = [
        ("boxes", "boxes"),
        ("temperature", "temperatures"),
        ("pressure", "pressures"),
    ]

    for filename, key in condition_files:
        if os.path.exists(filename):
            try:
                with open(filename, "r") as f:
                    for li, line in enumerate(f):
                        if li < len(list_dens):
                            idx = -li - 1 if unloading else li
                            conditions[key][list_dens[idx]] = float(line.strip())
                print(f"Loaded {filename} successfully")
            except Exception as e:
                print(f"Error reading {filename}: {e}")
        else:
            print(f"Skipping {filename} - file not found")

    return conditions


def main():
    # Configuration
    rootdir = "./"
    pattern = "p"
    unloading = False
    connectivity = sys.argv[1]

    # Files to process
    target_files = [
        "average_cluster_size.dat",
        "spanning_cluster_size.dat",
        "correlation_length.dat",
        "order_parameter.dat",
        "percolation_probability.dat",
        "largest_cluster_size.dat",
    ]

    # Find density directories
    dirs = natsorted(os.listdir(rootdir))
    list_dens = [d for d in dirs if pattern in d]

    if not list_dens:
        print(f"No directories matching pattern '{pattern}' found")
        return

    # Load condition data
    conditions = load_condition_files(list_dens, unloading)

    # Variables to export by (only include variables that have data)
    available_variables = ["density"]  # density is always available

    if conditions["pressures"]:  # Check if dictionary has any entries
        available_variables.append("pressure")
    if conditions["temperatures"]:  # Check if dictionary has any entries
        available_variables.append("temperature")
    if conditions["boxes"]:  # Check if dictionary has any entries
        available_variables.append("box")

    export_variables = available_variables

    print(f"Available export variables: {export_variables}")

    # Initialize data collection
    data_collection = DataCollection()

    # Process each density directory
    for dens_dir in list_dens:
        if not os.path.isdir(dens_dir):
            continue

        print(f"Processing directory: {dens_dir}")

        files = os.listdir(dens_dir)
        for filename in files:
            if filename in target_files:
                filepath = os.path.join(dens_dir, filename)

                # Parse the CSV file
                parsed_data = parse_csv_file(filepath)

                if parsed_data:
                    print(f"  Found {len(parsed_data)} data points in {filename}")

                    # Create data points and add to collection
                    for data_dict in parsed_data:
                        data_point = DataPoint(
                            dens_dir=dens_dir,
                            file_name=filename,
                            connectivity_type=data_dict["connectivity_type"],
                            concentration=data_dict["concentration"],
                            value=data_dict["value"],
                            error=data_dict["error"],
                        )

                        # Set experimental conditions (only if data exists)
                        data_point.set_conditions(
                            pressure=conditions["pressures"].get(dens_dir, 0.0),
                            temperature=conditions["temperatures"].get(dens_dir, 0.0),
                            box_size=conditions["boxes"].get(dens_dir, 0.0),
                        )

                        if data_dict["connectivity_type"] == connectivity:
                            data_collection.add_data_point(data_point)

    print(f"\nTotal data points collected: {len(data_collection.data_points)}")
    print(f"Files processed: {data_collection.get_unique_files()}")

    # Export data organized by different variables
    exporter = DataExporter(data_collection)

    for variable in export_variables:
        print(f"\nExporting data organized by {variable}...")
        exporter.export_by_variable(variable)

    # Export concentration files separately (each connectivity gets its own file)
    print(f"\nExporting individual concentration files...")
    exporter.export_concentration_files()

    print("\nExport complete!")


if __name__ == "__main__":
    main()
