from ...config.settings import Settings
from ...io.writer.base_writer import BaseWriter

import os
from typing import List, Dict, Tuple

class MultipleFilesSummaryWriter(BaseWriter):
    """
    Aggregates per-file analysis results into summary files.

    Walks the export directory tree, collects ``.dat`` result files grouped by
    analysis type, and writes combined summary files. Supports ``"all"`` mode
    (single summary per analysis type) and ``"connectivity"`` mode (one summary per
    connectivity label).

    Attributes:
        _mode (str): Output mode (``"all"`` or ``"connectivity"``).
    """

    def __init__(self, settings: Settings, mode: str = "all") -> None:
        """
        Initialize the summary writer.

        Args:
            settings (Settings): Configuration settings.
            mode (str): Output mode controlling summary file layout.
        """
        super().__init__(settings)
        self._settings: Settings = settings
        self._mode: str = mode

    def write(self) -> None:
        """Collect result files and write combined summaries."""
        average_cluster_size_files = []
        correlation_length_files = []
        largest_cluster_size_files = []
        order_parameter_files = []
        percolation_probability_files = []
        spanning_cluster_size_files = []

        for subdir, dirs, files in os.walk(self._settings.export_directory):
            dirs.sort()
            for file in files:
                if file.endswith(".dat") and "unwrapped_clusters" not in file and "summary" not in file:
                    if "average_cluster_size" in file:
                        average_cluster_size_files.append(os.path.join(subdir, file))
                    elif "correlation_length" in file:
                        correlation_length_files.append(os.path.join(subdir, file))
                    elif "largest_cluster_size" in file:
                        largest_cluster_size_files.append(os.path.join(subdir, file))
                    elif "order_parameter" in file:
                        order_parameter_files.append(os.path.join(subdir, file))
                    elif "percolation_probability" in file:
                        percolation_probability_files.append(os.path.join(subdir, file))
                    elif "spanning_cluster_size" in file:
                        spanning_cluster_size_files.append(os.path.join(subdir, file))

        average_cluster_size_results = {}
        correlation_length_results = {}
        largest_cluster_size_results = {}
        order_parameter_results = {}
        percolation_probability_results = {}
        spanning_cluster_size_results = {}
        
        if average_cluster_size_files:
            average_cluster_size_results = self._get_results(average_cluster_size_files)
            self._write_summary(average_cluster_size_results, len(average_cluster_size_files), "average_cluster_size_summary.dat", self._mode)
        if correlation_length_files:
            correlation_length_results = self._get_results(correlation_length_files)
            self._write_summary(correlation_length_results, len(correlation_length_files), "correlation_length_summary.dat", self._mode)
        if largest_cluster_size_files:
            largest_cluster_size_results = self._get_results(largest_cluster_size_files)
            self._write_summary(largest_cluster_size_results, len(largest_cluster_size_files), "largest_cluster_size_summary.dat", self._mode)
        if order_parameter_files:
            order_parameter_results = self._get_results(order_parameter_files)
            self._write_summary(order_parameter_results, len(order_parameter_files), "order_parameter_summary.dat", self._mode)
        if percolation_probability_files:
            percolation_probability_results = self._get_results(percolation_probability_files)
            self._write_summary(percolation_probability_results, len(percolation_probability_files), "percolation_probability_summary.dat", self._mode)
        if spanning_cluster_size_files:
            spanning_cluster_size_results = self._get_results(spanning_cluster_size_files)
            self._write_summary(spanning_cluster_size_results, len(spanning_cluster_size_files), "spanning_cluster_size_summary.dat", self._mode)


    def _get_results(self, files: List[str]) -> Dict[str, Tuple[float, float, float]]:
        """
        Parse result files and aggregate values by connectivity type.

        Args:
            files (List[str]): Paths to ``.dat`` result files to parse.

        Returns:
            Dict[str, Tuple[float, float, float]]: Mapping of connectivity labels to
                lists of (concentration, value, standard deviation) tuples.
        """
        results = {}
        for file in files:
            with open(file, 'r') as f:
                lines = f.readlines()
                for line in lines:
                    if line.startswith("#"):
                        continue
                    try:
                        parts = line.split(',')
                        connectivity_type = parts[0]
                        concentration = float(parts[1])
                        average_cluster_size = float(parts[2])
                        standard_deviation = float(parts[3])
                    except:
                        print(f"Error parsing line: {line}")
                        continue
                    if connectivity_type not in results:
                        results[connectivity_type] = []
                    results[connectivity_type].append((concentration, average_cluster_size, standard_deviation))
        return results

    def _write_summary(self, results: Dict[str, Tuple[float, float, float]], n_data: int, file_name: str, mode: str = "all") -> None:
        """
        Write a combined summary file from aggregated results.

        Args:
            results (Dict[str, Tuple[float, float, float]]): Aggregated results
                keyed by connectivity type.
            n_data (int): Number of data files that were aggregated.
            file_name (str): Output file name.
            mode (str): Output mode (``"all"`` or ``"connectivity"``).
        """
        if mode == "all":
            path = os.path.join(self._settings.export_directory, file_name)
            with open(path, 'w') as f:
                for i, connectivity_type in enumerate(results.keys()):
                    f.write(f"# {i+1} : {connectivity_type}\n")
                for j, connectivity_type in enumerate(results.keys()):
                    f.write(f"# {j+i} : std {connectivity_type}\n")
                    
                for i in range(n_data):
                    for connectivity_type in results.keys():
                        f.write(f"{results[connectivity_type][i][1]} ")
                    for j, connectivity_type in enumerate(results.keys()):
                        if j < len(results.keys()) - 1:
                            f.write(f"{results[connectivity_type][i][2]} ")
                        else:
                            f.write(f"{results[connectivity_type][i][2]}\n")
        elif mode == "connectivity":
            for connectivity_type in results.keys():
                path = os.path.join(self._settings.export_directory, file_name.replace('summary', f'{connectivity_type}_summary'))
                with open(path, 'w') as f:
                    f.write(f"# Connectivity_type : {connectivity_type}\n")
                    f.write(f"# Concentration Average_cluster_size Standard_deviation\n")
                    for i in range(n_data):
                        f.write(f"{results[connectivity_type][i][0]} ")
                        f.write(f"{results[connectivity_type][i][1]} ")
                        f.write(f"{results[connectivity_type][i][2]}\n")
