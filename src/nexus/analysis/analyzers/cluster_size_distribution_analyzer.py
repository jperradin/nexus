from typing import List, Dict
from ...core.frame import Frame
from .base_analyzer import BaseAnalyzer
from ...config.settings import Settings
from ...utils.aesthetics import remove_duplicate_lines
import numpy as np
import os
from datetime import datetime


class ClusterSizeDistributionAnalyzer(BaseAnalyzer):
    """
    Computes the cluster size distribution ``n(s)`` for each connectivity type.

    Counts how many clusters of each size *s* exist per connectivity across all
    frames. Percolating clusters are excluded from the distribution. Results
    are written to one file per connectivity.

    Attributes:
        _raw_size_distribution (Dict[str, Dict[int, List[int]]]): Per-frame
            counts keyed by connectivity then cluster size.
        _raw_concentrations (Dict[str, List[float]]): Per-frame concentrations.
        size_distribution (Dict[str, Dict[int, float]]): Aggregated total count
            per cluster size per connectivity.
        std (Dict[str, Dict[int, float]]): Standard deviation (ddof=1) per size.
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
        self._raw_size_distribution: Dict[str, Dict[int, List[int]]] = {}
        self._raw_concentrations: Dict[str, List[float]] = {}

        # Public attributes to hold the final, aggregated results
        self.size_distribution: Dict[str, Dict[int, float | np.float64]] = {}
        self.std: Dict[str, Dict[int, float | np.float64]] = {}
        self.concentrations: Dict[str, float | np.float64] = {}

        # A flag to ensure final calculations are only performed once
        self._finalized: bool = False

    def analyze(self, frame: Frame, connectivities: List[str]) -> None:
        """
        Count clusters by size for each connectivity in *frame*.

        Args:
            frame (Frame): The frame to analyze.
            connectivities (List[str]): Connectivity labels to analyze.
        """
        clusters = frame.get_clusters()
        concentrations = frame.get_concentration()

        for connectivity in connectivities:
            # Initialize dictionaries if this is the first time seeing this connectivity
            self._raw_size_distribution.setdefault(connectivity, {})
            self._raw_concentrations.setdefault(connectivity, [])

            sizes = [
                c.get_size()
                for c in clusters
                if c.get_connectivity() == connectivity and not c.is_percolating
            ]

            if sizes:
                unique_sizes, ns = np.unique(sizes, return_counts=True)
                for s, n in zip(unique_sizes, ns):
                    self._raw_size_distribution[connectivity].setdefault(s, []).append(
                        n
                    )

            # Record concentration for this frame
            self._raw_concentrations[connectivity].append(
                concentrations.get(connectivity, 0.0)
            )

        self.update_frame_processed()

    def finalize(self) -> Dict:
        """
        Aggregate size counts and compute standard deviations over all frames.

        Returns:
            Dict: The finalized results dictionary.
        """
        if self._finalized:
            return self.get_result()

        for connectivity, size_data in self._raw_size_distribution.items():
            self.size_distribution.setdefault(connectivity, {})
            self.std.setdefault(connectivity, {})
            for size, counts in size_data.items():
                # The final value is the total count
                total_count = np.sum(counts)
                num_frames = self.frame_processed_count
                self.size_distribution[connectivity][size] = (
                    total_count  # / num_frames if num_frames > 0 else 0.0
                )

                # To calculate std dev, we need to account for frames where a size didn't appear
                all_counts_for_size = counts + [0] * (num_frames - len(counts))
                if len(all_counts_for_size) > 1:
                    self.std[connectivity][size] = np.std(
                        all_counts_for_size, ddof=1)
                else:
                    self.std[connectivity][size] = 0.0

        for connectivity, concs in self._raw_concentrations.items():
            self.concentrations[connectivity] = np.mean(
                concs) if concs else 0.0

        self._finalized = True
        return self.get_result()

    def get_result(self) -> Dict[str, Dict]:
        """
        Return the current results dictionary.

        Returns:
            Dict[str, Dict]: Keys are ``"concentrations"``,
                ``"size_distribution"``, and ``"std"``.
        """
        return {
            "concentrations": self.concentrations,
            "size_distribution": self.size_distribution,
            "std": self.std,
        }

    def print_to_file(self) -> None:
        """Write results to one ``cluster_size_distribution-<conn>.dat`` per connectivity."""
        output = self.finalize()

        for connectivity in self.size_distribution:
            self._write_header(connectivity)
            path = os.path.join(
                self._settings.export_directory,
                f"cluster_size_distribution-{connectivity}.dat",
            )

            # Sort by size in descending order for plotting
            sorted_sizes = sorted(
                self.size_distribution[connectivity].keys(), reverse=True
            )

            with open(path, "a") as f:
                for size in sorted_sizes:
                    concentration = output["concentrations"].get(
                        connectivity, 0.0)
                    n_s = output["size_distribution"][connectivity].get(
                        size, 0.0)
                    std_dev = output["std"][connectivity].get(size, 0.0)
                    f.write(
                        f"{connectivity},{concentration},{size},{n_s},{std_dev}\n")
            remove_duplicate_lines(path)

    def _write_header(self, connectivity: str) -> None:
        """
        Write the CSV header to the output file for a given connectivity.

        Args:
            connectivity (str): Connectivity label for the output file.
        """
        path = os.path.join(
            self._settings.export_directory,
            f"cluster_size_distribution-{connectivity}.dat",
        )
        number_of_frames = self.frame_processed_count

        if self._settings.analysis.overwrite or not os.path.exists(path):
            mode = "w"
        else:
            if os.path.getsize(path) > 0:
                return
            mode = "a"

        with open(path, mode, encoding="utf-8") as output:
            output.write(
                f"# Cluster Size Distribution Results for {connectivity}\n")
            output.write(f"# Date: {datetime.now()}\n")
            output.write(f"# Frames averaged: {number_of_frames}\n")
            output.write(
                "# Connectivity_type,Concentration,Cluster_size,N_clusters,Standard_deviation_ddof=1\n"
            )

    def __str__(self) -> str:
        """Return the class name."""
        return f"{self.__class__.__name__}"

    def __repr__(self) -> str:
        """Return a reproducible string representation."""
        return f"{self.__class__.__name__}()"
