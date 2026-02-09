from typing import List, Dict
from ...core.frame import Frame
from .base_analyzer import BaseAnalyzer
from ...config.settings import Settings
from ...utils.aesthetics import remove_duplicate_lines

import numpy as np
import os
from datetime import datetime


class OrderParameterAnalyzer(BaseAnalyzer):
    """
    Computes the percolation order parameter (P_inf) for each connectivity type.

    The order parameter is the fraction of networking nodes that belong to the
    percolating cluster. 

    Attributes:
        _raw_order_parameters (Dict[str, List[float]]): Per-frame P_inf values.
        _raw_concentrations (Dict[str, List[float]]): Per-frame concentrations.
        order_parameters (Dict[str, float]): Ensemble-averaged P_inf.
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
        self._raw_order_parameters: Dict[str, List[float]] = {}
        self._raw_concentrations: Dict[str, List[float]] = {}

        # Public attributes to hold the final, aggregated results
        self.order_parameters: Dict[str, float | np.float64] = {}
        self.std: Dict[str, float | np.float64] = {}
        self.error: Dict[str, float | np.float64] = {}
        self.concentrations: Dict[str, float | np.float64] = {}

        # A flag to ensure final calculations are only performed once
        self._finalized: bool = False

    def analyze(self, frame: Frame, connectivities: List[str]) -> None:
        """
        Extract the order parameter from the percolating cluster per connectivity.

        Args:
            frame (Frame): The frame to analyze.
            connectivities (List[str]): Connectivity labels to analyze.
        """
        clusters = frame.get_clusters()
        concentrations = frame.get_concentration()

        for connectivity in connectivities:
            # Initialize lists if this is the first time seeing this connectivity
            self._raw_order_parameters.setdefault(connectivity, [])
            self._raw_concentrations.setdefault(connectivity, [])

            # The order parameter is typically defined for the largest cluster if it percolates.
            # We find the largest percolating cluster.
            percolating_clusters = [
                c
                for c in clusters
                if c.get_connectivity() == connectivity and c.is_percolating
            ]

            if percolating_clusters:
                # Assuming the order parameter is associated with the largest of the percolating clusters
                largest_perc_cluster = max(percolating_clusters, key=lambda c: c.size)
                # We are interested in the 1D order parameter (P∞_x) as per the README
                self._raw_order_parameters[connectivity].append(
                    largest_perc_cluster.order_parameter[0]
                )
            else:
                # If no cluster percolates, the order parameter is 0
                self._raw_order_parameters[connectivity].append(0.0)

            self._raw_concentrations[connectivity].append(
                concentrations.get(connectivity, 0.0)
            )

        self.update_frame_processed()

    def finalize(self) -> Dict[str, Dict[str, float | np.float64]]:
        """
        Compute ensemble averages of P_inf over all processed frames.

        Returns:
            Dict[str, Dict[str, float]]: The finalized results dictionary.
        """
        if self._finalized:
            return self.get_result()

        for connectivity, params in self._raw_order_parameters.items():
            if params:
                self.order_parameters[connectivity] = np.mean(params)
                if len(params) > 1:
                    self.std[connectivity] = np.std(params, ddof=1)
                    self.error[connectivity] = self.std[connectivity] / np.sqrt(
                        len(params)
                    )
                else:
                    self.std[connectivity] = 0.0
                    self.error[connectivity] = 0.0
            else:
                self.order_parameters[connectivity] = 0.0
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
                ``"order_parameters"``, ``"std"``, and ``"error"``.
        """
        return {
            "concentrations": self.concentrations,
            "order_parameters": self.order_parameters,
            "std": self.std,
            "error": self.error,
        }

    def print_to_file(self) -> None:
        """Write ensemble-averaged results to ``order_parameter.dat``."""
        output = self.finalize()
        self._write_header()
        path = os.path.join(self._settings.export_directory, "order_parameter.dat")
        with open(path, "a") as f:
            for connectivity in self.order_parameters:
                concentration = output["concentrations"].get(connectivity, 0.0)
                order_parameter = output["order_parameters"].get(connectivity, 0.0)
                std = output["std"].get(connectivity, 0.0)
                error = output["error"].get(connectivity, 0.0)
                f.write(
                    f"{connectivity},{concentration},{order_parameter},{std},{error}\n"
                )
        remove_duplicate_lines(path)

    def _write_header(self) -> None:
        """Write the CSV header to the output file if needed."""
        path = os.path.join(self._settings.export_directory, "order_parameter.dat")
        number_of_frames = self.frame_processed_count

        if self._settings.analysis.overwrite or not os.path.exists(path):
            mode = "w"
        else:
            if os.path.getsize(path) > 0:
                return
            mode = "a"

        with open(path, mode, encoding="utf-8") as output:
            output.write("# Order Parameter Results\n")
            output.write(f"# Date: {datetime.now()}\n")
            output.write(f"# Frames averaged: {number_of_frames}\n")
            output.write(
                "# Connectivity_type,Concentration,Order_parameter,Standard_deviation_ddof=1,Standard_error_ddof=1\n"
            )

    def __str__(self) -> str:
        """Return the class name."""
        return f"{self.__class__.__name__}"

    def __repr__(self) -> str:
        """Return a reproducible string representation."""
        return f"{self.__class__.__name__}()"
