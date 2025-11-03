from typing import List, Dict
from ...core.frame import Frame
from .base_analyzer import BaseAnalyzer
from ...config.settings import Settings
from ...utils.aesthetics import remove_duplicate_lines

import numpy as np
import os
from datetime import datetime


class SpanningClusterSizeAnalyzer(BaseAnalyzer):
    """
    Computes the size of the largest finite (non-percolating) cluster.

    This is a crucial metric in percolation theory, often denoted as S_span
    for the sub-critical regime. It helps characterize the size of the
    largest clusters that do not span the entire system.
    """

    def __init__(self, settings: Settings) -> None:
        """Initializes the analyzer."""
        super().__init__(settings)
        # Private attributes to store raw, per-frame data
        self._raw_spanning_sizes: Dict[str, List[float]] = {}
        self._raw_concentrations: Dict[str, List[float]] = {}

        # Public attributes to hold the final, aggregated results
        self.spanning_cluster_sizes: Dict[str, float] = {}
        self.std: Dict[str, float] = {}
        self.concentrations: Dict[str, float] = {}
        self.error: Dict[str, float] = {}

        # A flag to ensure final calculations are only performed once
        self._finalized: bool = False

    def analyze(self, frame: Frame, connectivities: List[str]) -> None:
        """
        Analyzes a single frame to find the largest non-percolating cluster
        size for each connectivity type and stores the raw data.
        """
        clusters = frame.get_clusters()
        concentrations = frame.get_concentration()

        for connectivity in connectivities:
            # Initialize lists if this is the first time seeing this connectivity
            self._raw_spanning_sizes.setdefault(connectivity, [])
            self._raw_concentrations.setdefault(connectivity, [])

            # Filter for non-percolating clusters and get their sizes
            sizes = [
                c.get_size()
                for c in clusters
                if c.get_connectivity() == connectivity and not c.is_percolating
            ]

            if sizes:
                self._raw_spanning_sizes[connectivity].append(np.max(sizes))
            else:
                # If no non-percolating clusters exist, the size is 0
                self._raw_spanning_sizes[connectivity].append(0.0)

            self._raw_concentrations[connectivity].append(
                concentrations.get(connectivity, 0.0)
            )

        self.update_frame_processed()

    def finalize(self) -> Dict[str, Dict[str, float]]:
        """
        Calculates the final mean, standard deviation, and fluctuation for the
        spanning cluster size across all processed frames. This method is now idempotent.
        """
        if self._finalized:
            return self.get_result()

        for connectivity, sizes in self._raw_spanning_sizes.items():
            if sizes:
                self.spanning_cluster_sizes[connectivity] = np.mean(sizes)
                if len(sizes) > 1:
                    self.std[connectivity] = np.std(sizes, ddof=1)
                    mean_size = self.spanning_cluster_sizes[connectivity]
                    self.error[connectivity] = self.std[connectivity] / np.sqrt(
                        len(sizes)
                    )
                else:
                    self.std[connectivity] = 0.0
                    self.error[connectivity] = 0.0
            else:
                self.spanning_cluster_sizes[connectivity] = 0.0
                self.std[connectivity] = 0.0
                self.error[connectivity] = 0.0

            self.std[connectivity] = np.nan_to_num(self.std[connectivity])
            self.error[connectivity] = np.nan_to_num(self.error[connectivity])

        for connectivity, concs in self._raw_concentrations.items():
            self.concentrations[connectivity] = np.mean(concs) if concs else 0.0

        self._finalized = True
        return self.get_result()

    def get_result(self) -> Dict[str, Dict[str, float]]:
        """Returns the finalized analysis results."""
        return {
            "concentrations": self.concentrations,
            "spanning_cluster_size": self.spanning_cluster_sizes,
            "std": self.std,
            "error": self.error,
        }

    def print_to_file(self) -> None:
        """Writes the finalized results to a data file."""
        output = self.finalize()
        self._write_header()
        path = os.path.join(
            self._settings.export_directory, "spanning_cluster_size.dat"
        )
        with open(path, "a") as f:
            for connectivity in self.spanning_cluster_sizes:
                concentration = output["concentrations"].get(connectivity, 0.0)
                spanning_size = output["spanning_cluster_size"].get(connectivity, 0.0)
                std = output["std"].get(connectivity, 0.0)
                error = output["error"].get(connectivity, 0.0)
                f.write(
                    f"{connectivity},{concentration},{spanning_size},{std},{error}\n"
                )
        remove_duplicate_lines(path)

    def _write_header(self) -> None:
        """Initializes the output file with a header."""
        path = os.path.join(
            self._settings.export_directory, "spanning_cluster_size.dat"
        )
        number_of_frames = self.frame_processed_count

        if self._settings.analysis.overwrite or not os.path.exists(path):
            mode = "w"
        else:
            if os.path.getsize(path) > 0:
                return
            mode = "a"

        with open(path, mode, encoding="utf-8") as output:
            output.write(f"# Spanning Cluster Size Results\n")
            output.write(f"# Date: {datetime.now()}\n")
            output.write(f"# Frames averaged: {number_of_frames}\n")
            output.write(
                "# Connectivity_type,Concentration,Spanning_cluster_size,Standard_deviation_ddof=1,Standard_error_ddof=1\n"
            )

    def __str__(self) -> str:
        return f"{self.__class__.__name__}"

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}()"
