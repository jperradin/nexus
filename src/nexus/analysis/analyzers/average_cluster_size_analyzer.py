from typing import List, Dict
from ...core.frame import Frame
from .base_analyzer import BaseAnalyzer
from ...config.settings import Settings
from ...utils.aesthetics import remove_duplicate_lines

import numpy as np
import os
from datetime import datetime

class AverageClusterSizeAnalyzer(BaseAnalyzer):
    """
    Computes the weight-average cluster size ``<S>``.

    Uses the formula ``<S> = sum(s^2 * n(s)) / sum(s * n(s))`` where *s* is
    the cluster size and *n(s)* the number of clusters of size *s*. This
    quantity diverges at the percolation threshold. Percolating clusters are
    excluded to focus on the finite-cluster distribution.

    Attributes:
        _raw_average_sizes (Dict[str, List[float]]): Per-frame ``<S>`` values
            keyed by connectivity label.
        _raw_concentrations (Dict[str, List[float]]): Per-frame concentrations.
        average_sizes (Dict[str, float]): Ensemble-averaged ``<S>`` per connectivity.
        std (Dict[str, float]): Standard deviation (ddof=1) of ``<S>``.
        error (Dict[str, float]): Standard error of ``<S>``.
        concentrations (Dict[str, float]): Mean concentration per connectivity.
    """

    def __init__(self, settings: Settings) -> None:
        """
        Initialize the analyzer.

        Args:
            settings (Settings): Configuration settings.
        """
        super().__init__(settings)
        # Private attributes to store raw, per-frame data
        self._raw_average_sizes: Dict[str, List[float | np.float64]] = {}
        self._raw_concentrations: Dict[str, List[float | np.float64]] = {}

        # Public attributes to hold the final, aggregated results
        self.average_sizes: Dict[str, float | np.float64] = {}
        self.std: Dict[str, float | np.float64] = {}
        self.error: Dict[str, float | np.float64] = {}
        self.concentrations: Dict[str, float | np.float64] = {}

        # A flag to ensure final calculations are only performed once
        self._finalized: bool = False

    def analyze(self, frame: Frame, connectivities: List[str]) -> None:
        """
        Compute the weight-average cluster size for each connectivity in *frame*.

        Args:
            frame (Frame): The frame to analyze.
            connectivities (List[str]): Connectivity labels to analyze.
        """
        clusters = frame.get_clusters()
        concentrations = frame.get_concentration()

        for connectivity in connectivities:
            # Initialize lists if this is the first time seeing this connectivity
            self._raw_average_sizes.setdefault(connectivity, [])
            self._raw_concentrations.setdefault(connectivity, [])

            sizes = [
                c.get_size()
                for c in clusters
                if c.get_connectivity() == connectivity and not c.is_percolating
            ]

            if sizes:
                unique_sizes, ns = np.unique(sizes, return_counts=True)
                numerator = np.sum(unique_sizes**2 * ns)
                denominator = np.sum(unique_sizes * ns)

                average_size = numerator / denominator if denominator > 0 else 0.0
                self._raw_average_sizes[connectivity].append(average_size)
            else:
                self._raw_average_sizes[connectivity].append(0.0)

            self._raw_concentrations[connectivity].append(
                concentrations.get(connectivity, 0.0)
            )

        self.update_frame_processed()

    def finalize(self) -> Dict[str, Dict[str, float | np.float64]]:
        """
        Compute ensemble averages of ``<S>`` over all processed frames.

        Returns:
            Dict[str, Dict[str, float]]: The finalized results dictionary.
        """
        if self._finalized:
            return self.get_result()

        for connectivity, sizes in self._raw_average_sizes.items():
            if sizes:
                self.average_sizes[connectivity] = np.mean(sizes)
                if len(sizes) > 1:
                    self.std[connectivity] = np.std(sizes, ddof=1)
                    self.error[connectivity] = self.std[connectivity] / np.sqrt(
                        len(sizes)
                    )
                else:
                    self.std[connectivity] = 0.0
                    self.error[connectivity] = 0.0
            else:
                self.average_sizes[connectivity] = 0.0
                self.std[connectivity] = 0.0
                self.error[connectivity] = 0.0

            self.std[connectivity] = np.nan_to_num(self.std[connectivity])
            self.error[connectivity] = np.nan_to_num(self.error[connectivity])

        for connectivity, concs in self._raw_concentrations.items():
            self.concentrations[connectivity] = np.mean(concs) if concs else 0.0

        self._finalized = True
        return self.get_result()

    def get_result(self) -> Dict[str, Dict[str, float | np.float64]]:
        """
        Return the current results dictionary.

        Returns:
            Dict[str, Dict[str, float]]: Keys are ``"concentrations"``,
                ``"average_cluster_size"``, ``"std"``, and ``"error"``.
        """
        return {
            "concentrations": self.concentrations,
            "average_cluster_size": self.average_sizes,
            "std": self.std,
            "error": self.error,
        }

    def print_to_file(self) -> None:
        """Write ensemble-averaged results to ``average_cluster_size.dat``."""
        output = self.finalize()
        self._write_header()
        path = os.path.join(self._settings.export_directory, "average_cluster_size.dat")
        with open(path, "a") as f:
            for connectivity in self.average_sizes:
                concentration = output["concentrations"].get(connectivity, 0.0)
                average_size = output["average_cluster_size"].get(connectivity, 0.0)
                std = output["std"].get(connectivity, 0.0)
                error = output["error"].get(connectivity, 0.0)
                f.write(
                    f"{connectivity},{concentration},{average_size},{std},{error}\n"
                )
        remove_duplicate_lines(path)

    def _write_header(self) -> None:
        """Write the CSV header to the output file if needed."""
        path = os.path.join(self._settings.export_directory, "average_cluster_size.dat")
        number_of_frames = self.frame_processed_count

        if self._settings.analysis.overwrite or not os.path.exists(path):
            mode = "w"
        else:
            if os.path.getsize(path) > 0:
                return
            mode = "a"

        with open(path, mode, encoding="utf-8") as output:
            output.write("# Average Cluster Size Results\n")
            output.write(f"# Date: {datetime.now()}\n")
            output.write(f"# Frames averaged: {number_of_frames}\n")
            output.write(
                "# Connectivity_type,Concentration,Average_cluster_size,Standard_deviation_ddof=1,Standard_error_ddof=1\n"
            )

    def __str__(self) -> str:
        """Return the class name."""
        return f"{self.__class__.__name__}"

    def __repr__(self) -> str:
        """Return a reproducible string representation."""
        return f"{self.__class__.__name__}()"
