from typing import List, Dict
from ...core.frame import Frame
from .base_analyzer import BaseAnalyzer
from ...config.settings import Settings
from ...utils.aesthetics import remove_duplicate_lines
import numpy as np
import os
from datetime import datetime


class CorrelationLengthAnalyzer(BaseAnalyzer):
    """
    Computes the correlation length (ξ) of the cluster size distribution.

    The correlation length is a measure of the characteristic size of clusters
    and is calculated using the second moment of the gyration radius distribution,
    weighted by cluster size: ξ² = Σ(2 * R_s² * s² * n_s) / Σ(s² * n_s),
    where R_s is the gyration radius of clusters of size s.
    """

    def __init__(self, settings: Settings) -> None:
        """Initializes the analyzer."""
        super().__init__(settings)
        # Private attributes to store raw, per-frame data
        self._raw_correlation_lengths: Dict[str, List[float]] = {}
        self._raw_concentrations: Dict[str, List[float]] = {}

        # Public attributes to hold the final, aggregated results
        self.correlation_length: Dict[str, float] = {}
        self.std: Dict[str, float] = {}
        self.concentrations: Dict[str, float] = {}
        self.error: Dict[str, float] = {}

        # A flag to ensure final calculations are only performed once
        self._finalized: bool = False

    def analyze(self, frame: Frame, connectivities: List[str]) -> None:
        """
        Analyzes a single frame to compute the correlation length for each
        connectivity type and stores the raw data.
        """
        clusters = frame.get_clusters()
        concentrations = frame.get_concentration()

        for connectivity in connectivities:
            # Initialize lists if this is the first time seeing this connectivity
            self._raw_correlation_lengths.setdefault(connectivity, [])
            self._raw_concentrations.setdefault(connectivity, [])

            non_percolating_clusters = [
                c
                for c in clusters
                if c.get_connectivity() == connectivity and not c.is_percolating
            ]

            if non_percolating_clusters:
                sizes = np.array([c.get_size() for c in non_percolating_clusters])
                gyration_radii_sq = np.array(
                    [c.gyration_radius**2 for c in non_percolating_clusters]
                )

                numerator = np.sum(2 * gyration_radii_sq * sizes**2)
                denominator = np.sum(sizes**2)

                correlation_length_sq = (
                    numerator / denominator if denominator > 0 else 0.0
                )
                self._raw_correlation_lengths[connectivity].append(
                    np.sqrt(correlation_length_sq)
                )
            else:
                self._raw_correlation_lengths[connectivity].append(0.0)

            self._raw_concentrations[connectivity].append(
                concentrations.get(connectivity, 0.0)
            )

        self.update_frame_processed()

    def finalize(self) -> Dict[str, Dict[str, float]]:
        """
        Calculates the final mean, standard deviation, and fluctuation for the
        correlation length across all processed frames. This method is now idempotent.
        """
        if self._finalized:
            return self.get_result()

        for connectivity, lengths in self._raw_correlation_lengths.items():
            if lengths:
                self.correlation_length[connectivity] = np.mean(lengths)
                if len(lengths) > 1:
                    self.std[connectivity] = np.std(lengths, ddof=1)
                    mean_length = self.correlation_length[connectivity]
                    self.error[connectivity] = self.std[connectivity] / np.sqrt(
                        len(lengths)
                    )
                else:
                    self.std[connectivity] = 0.0
                    self.error[connectivity] = 0.0
            else:
                self.correlation_length[connectivity] = 0.0
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
            "correlation_length": self.correlation_length,
            "std": self.std,
            "error": self.error,
        }

    def print_to_file(self) -> None:
        """Writes the finalized results to a data file."""
        output = self.finalize()
        self._write_header()
        path = os.path.join(self._settings.export_directory, "correlation_length.dat")
        with open(path, "a") as f:
            for connectivity in self.correlation_length:
                concentration = output["concentrations"].get(connectivity, 0.0)
                correlation_length = output["correlation_length"].get(connectivity, 0.0)
                std = output["std"].get(connectivity, 0.0)
                error = output["error"].get(connectivity, 0.0)
                f.write(
                    f"{connectivity},{concentration},{correlation_length},{std},{error}\n"
                )
        remove_duplicate_lines(path)

    def _write_header(self) -> None:
        """Initializes the output file with a header."""
        path = os.path.join(self._settings.export_directory, "correlation_length.dat")
        number_of_frames = self.frame_processed_count

        if self._settings.analysis.overwrite or not os.path.exists(path):
            mode = "w"
        else:
            if os.path.getsize(path) > 0:
                return
            mode = "a"

        with open(path, mode, encoding="utf-8") as output:
            output.write(f"# Correlation Length Results\n")
            output.write(f"# Date: {datetime.now()}\n")
            output.write(f"# Frames averaged: {number_of_frames}\n")
            output.write(
                "# Connectivity_type,Concentration,Correlation_length,Standard_deviation_ddof=1,Standard_error_ddof=1\n"
            )

    def __str__(self) -> str:
        return f"{self.__class__.__name__}"

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}()"
