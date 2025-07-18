from typing import List, Dict
from ...core.frame import Frame
from .base_analyzer import BaseAnalyzer
from ...config.settings import Settings
from ...utils.aesthetics import remove_duplicate_lines

import numpy as np
import os
from datetime import datetime


class LargestClusterSizeAnalyzer(BaseAnalyzer):
    """
    Analyzes the size of the largest cluster for each connectivity type,
    providing insights into the dominant structures in the system.
    """

    def __init__(self, settings: Settings) -> None:
        super().__init__(settings)
        # Use private attributes to store raw, per-frame data
        self._raw_sizes: Dict[str, List[float]] = {}
        self._raw_concentrations: Dict[str, List[float]] = {}

        # Public attributes to hold the final, aggregated results
        self.largest_cluster_sizes: Dict[str, float] = {}
        self.std: Dict[str, float] = {}
        self.concentrations: Dict[str, float] = {}
        self.fluctuations: Dict[str, float] = {}

        # A flag to ensure final calculations are only performed once
        self._finalized: bool = False

    def analyze(self, frame: Frame, connectivities: List[str]) -> None:
        """
        Analyzes a single frame to find the largest cluster size for each
        connectivity type and stores the raw data.
        """
        clusters = frame.get_clusters()
        concentrations = frame.get_concentration()

        for connectivity in connectivities:
            # Initialize lists if this is the first time seeing this connectivity
            self._raw_sizes.setdefault(connectivity, [])
            self._raw_concentrations.setdefault(connectivity, [])

            sizes = [
                c.get_size() for c in clusters if c.get_connectivity() == connectivity
            ]
            if sizes:
                self._raw_sizes[connectivity].append(np.max(sizes))
                self._raw_concentrations[connectivity].append(
                    concentrations.get(connectivity, 0.0)
                )
            else:
                self._raw_sizes[connectivity].append(0.0)
                self._raw_concentrations[connectivity].append(
                    concentrations.get(connectivity, 0.0)
                )

        self.update_frame_processed()

    def finalize(self) -> Dict[str, Dict[str, float]]:
        """
        Calculates the final mean, standard deviation, and fluctuation for all
        processed frames. This method is now idempotent and will only compute
        the results once.
        """
        if self._finalized:
            return self.get_result()

        for connectivity, sizes in self._raw_sizes.items():
            if sizes:
                self.largest_cluster_sizes[connectivity] = np.mean(sizes)
                if len(sizes) > 1:
                    self.std[connectivity] = np.std(sizes, ddof=1)
                    mean_size = self.largest_cluster_sizes[connectivity]
                    # Avoid division by zero for fluctuations
                    self.fluctuations[connectivity] = (
                        np.var(sizes, ddof=1) / mean_size if mean_size > 0 else 0.0
                    )
                else:
                    self.std[connectivity] = 0.0
                    self.fluctuations[connectivity] = 0.0
            else:
                self.largest_cluster_sizes[connectivity] = 0.0
                self.std[connectivity] = 0.0
                self.fluctuations[connectivity] = 0.0

            # Replace potential NaN values with 0.0 for robustness
            self.std[connectivity] = np.nan_to_num(self.std[connectivity])
            self.fluctuations[connectivity] = np.nan_to_num(
                self.fluctuations[connectivity]
            )

        for connectivity, concs in self._raw_concentrations.items():
            self.concentrations[connectivity] = np.mean(concs) if concs else 0.0

        self._finalized = True
        return self.get_result()

    def get_result(self) -> Dict[str, Dict[str, float]]:
        """Returns the finalized analysis results."""
        return {
            "concentrations": self.concentrations,
            "largest_cluster_size": self.largest_cluster_sizes,
            "std": self.std,
            "fluctuations": self.fluctuations,
        }

    def print_to_file(self) -> None:
        """Writes the finalized results to a data file."""
        # Ensure results are calculated before printing
        output = self.finalize()

        self._write_header()
        path = os.path.join(self._settings.export_directory, "largest_cluster_size.dat")
        with open(path, "a") as f:
            for connectivity in self.largest_cluster_sizes:
                concentration = output["concentrations"].get(connectivity, 0.0)
                largest_size = output["largest_cluster_size"].get(connectivity, 0.0)
                std_dev = output["std"].get(connectivity, 0.0)
                fluct = output["fluctuations"].get(connectivity, 0.0)
                f.write(
                    f"{connectivity},{concentration},{largest_size},{std_dev},{fluct}\n"
                )
        remove_duplicate_lines(path)

    def _write_header(self) -> None:
        """Initializes the output file with a header."""
        path = os.path.join(self._settings.export_directory, "largest_cluster_size.dat")
        number_of_frames = self.frame_processed_count
        file_exists = os.path.exists(path)

        # Write header only if the file is new or overwrite is enabled
        if self._settings.analysis.overwrite or not file_exists:
            mode = "w"
        else:  # If appending, we might still need a header for a new run
            # A simple check: if file is not empty, assume header exists.
            # A more robust solution could involve checking the last line.
            if os.path.getsize(path) > 0:
                return  # Assume header exists, do nothing.
            mode = "a"

        with open(path, mode, encoding="utf-8") as output:
            output.write(f"# Largest Cluster Size Results\n")
            output.write(f"# Date: {datetime.now()}\n")
            output.write(f"# Frames averaged: {number_of_frames}\n")
            output.write(
                "# Connectivity_type,Concentration,Largest_cluster_size,Standard_deviation_ddof=1,Fluctuations_ddof=1\n"
            )

    def __str__(self) -> str:
        return f"{self.__class__.__name__}"

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}()"
