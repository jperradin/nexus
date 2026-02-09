from typing import List, Dict
from ...core.frame import Frame
from .base_analyzer import BaseAnalyzer
from ...config.settings import Settings
from ...utils.aesthetics import remove_duplicate_lines

import numpy as np
import os
from datetime import datetime


class PercolationProbabilityAnalyzer(BaseAnalyzer):
    """
    Computes the percolation probability (Pi) for each connectivity type.

    The percolation probability is the fraction of frames in which at least
    one cluster spans the simulation box *in all three dimensions*. A value of
    1.0 indicates percolation in every frame; 0.0 indicates it never occurs.

    Attributes:
        _raw_percolation_prob_x (Dict[str, List[float]]): Per-frame binary
            indicators (0 or 1) for percolation.
        _raw_concentrations (Dict[str, List[float]]): Per-frame concentrations.
        percolation_probabilities (Dict[str, float]): Ensemble-averaged Pi.
        std (Dict[str, float]): Standard deviation (ddof=1).
        error (Dict[str, float]): Standard error.
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
        # We store a boolean (0 or 1) for each frame indicating if percolation occurred
        self._raw_percolation_prob_x: Dict[str, List[float]] = {}
        self._raw_concentrations: Dict[str, List[float]] = {}

        # Public attributes to hold the final, aggregated results
        self.percolation_probabilities: Dict[str, float | np.float64] = {}
        self.std: Dict[str, float | np.float64] = {}
        self.error: Dict[str, float | np.float64] = {}
        self.concentrations: Dict[str, float | np.float64] = {}

        # A flag to ensure final calculations are only performed once
        self._finalized: bool = False

    def analyze(self, frame: Frame, connectivities: List[str]) -> None:
        """
        Check whether percolation occurs for each connectivity in *frame*.

        Args:
            frame (Frame): The frame to analyze.
            connectivities (List[str]): Connectivity labels to analyze.
        """
        clusters = frame.get_clusters()
        concentrations = frame.get_concentration()

        for connectivity in connectivities:
            # Initialize lists if this is the first time seeing this connectivity
            self._raw_percolation_prob_x.setdefault(connectivity, [])
            self._raw_concentrations.setdefault(connectivity, [])

            # Check if any cluster for this connectivity percolates in the 3 directions
            found_percolating_cluster = any(
                "xyz" in c.percolation_probability
                for c in clusters
                if c.get_connectivity() == connectivity
            )

            self._raw_percolation_prob_x[connectivity].append(
                1.0 if found_percolating_cluster else 0.0
            )
            self._raw_concentrations[connectivity].append(
                concentrations.get(connectivity, 0.0)
            )

        self.update_frame_processed()

    def finalize(self) -> Dict[str, Dict[str, float | np.float64]]:
        """
        Compute ensemble-averaged percolation probability over all frames.

        Returns:
            Dict[str, Dict[str, float]]: The finalized results dictionary.
        """
        if self._finalized:
            return self.get_result()

        for connectivity, probs in self._raw_percolation_prob_x.items():
            if probs:
                self.percolation_probabilities[connectivity] = np.mean(probs)
                if len(probs) > 1:
                    self.std[connectivity] = np.std(probs, ddof=1)
                    self.error[connectivity] = self.std[connectivity] / np.sqrt(
                        len(probs)
                    )
                else:
                    self.std[connectivity] = 0.0
                    self.error[connectivity] = 0.0
            else:
                self.percolation_probabilities[connectivity] = 0.0
                self.std[connectivity] = 0.0
                self.error[connectivity] = 0.0

            self.std[connectivity] = np.nan_to_num(self.std[connectivity])

        for connectivity, concs in self._raw_concentrations.items():
            self.concentrations[connectivity] = np.mean(concs) if concs else 0.0

        self._finalized = True
        return self.get_result()

    def get_result(self) -> Dict[str, Dict[str, float | np.float64]]:
        """
        Return the current results dictionary.

        Returns:
            Dict[str, Dict[str, float]]: Keys are ``"concentrations"``,
                ``"percolation_probabilities"``, ``"std"``, and ``"error"``.
        """
        return {
            "concentrations": self.concentrations,
            "percolation_probabilities": self.percolation_probabilities,
            "std": self.std,
            "error": self.error,
        }

    def print_to_file(self) -> None:
        """Write ensemble-averaged results to ``percolation_probability.dat``."""
        output = self.finalize()
        self._write_header()
        path = os.path.join(
            self._settings.export_directory, "percolation_probability.dat"
        )
        with open(path, "a") as f:
            for connectivity in self.percolation_probabilities:
                concentration = output["concentrations"].get(connectivity, 0.0)
                percolation_probability = output["percolation_probabilities"].get(
                    connectivity, 0.0
                )
                std = output["std"].get(connectivity, 0.0)
                error = output["error"].get(connectivity, 0.0)
                f.write(
                    f"{connectivity},{concentration},{percolation_probability},{std},{error}\n"
                )
        remove_duplicate_lines(path)

    def _write_header(self) -> None:
        """Write the CSV header to the output file if needed."""
        path = os.path.join(
            self._settings.export_directory, "percolation_probability.dat"
        )
        number_of_frames = self.frame_processed_count

        if self._settings.analysis.overwrite or not os.path.exists(path):
            mode = "w"
        else:
            if os.path.getsize(path) > 0:
                return
            mode = "a"

        with open(path, mode, encoding="utf-8") as output:
            output.write("# Percolation Probability Results\n")
            output.write(f"# Date: {datetime.now()}\n")
            output.write(f"# Frames averaged: {number_of_frames}\n")
            output.write(
                "# Connectivity_type,Concentration,Percolation_probability,Standard_deviation_ddof=1,Standard_error_ddof=1\n"
            )

    def __str__(self) -> str:
        """Return the class name."""
        return f"{self.__class__.__name__}"

    def __repr__(self) -> str:
        """Return a reproducible string representation."""
        return f"{self.__class__.__name__}()"
