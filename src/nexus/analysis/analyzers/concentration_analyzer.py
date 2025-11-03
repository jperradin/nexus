from typing import List, Dict
from ...core.frame import Frame
from .base_analyzer import BaseAnalyzer
from ...config.settings import Settings
from ...utils.aesthetics import remove_duplicate_lines

import numpy as np
import os
from datetime import datetime


class ConcentrationAnalyzer(BaseAnalyzer):
    """
    Computes the concentration of clusters for each connectivity type.

    This analyzer tracks the concentration of clusters for each connectivity type
    across all processed frames. The concentration is defined as the ratio of the
    number of nodes in clusters of a given connectivity type to the total number of nodes.
    """

    def __init__(self, settings: Settings) -> None:
        """Initializes the analyzer."""
        super().__init__(settings)
        # Private attribute to store raw, per-frame data
        self._raw_concentrations: Dict[str, List[float]] = {}

        # Public attributes to hold the final, aggregated results
        self.concentrations: Dict[str, float] = {}
        self.std: Dict[str, float] = {}
        self.error: Dict[str, float] = {}

        # A flag to ensure final calculations are only performed once
        self._finalized: bool = False

    def analyze(self, frame: Frame, connectivities: List[str]) -> None:
        """
        Analyzes a single frame to get the concentration for each connectivity
        type and stores the raw data.
        """
        concentrations = frame.get_concentration()
        for connectivity in connectivities:
            # Initialize list if this is the first time seeing this connectivity
            self._raw_concentrations.setdefault(connectivity, [])

            # Record concentration for this frame, defaulting to 0.0 if not present
            self._raw_concentrations[connectivity].append(
                concentrations.get(connectivity, 0.0)
            )

        self.update_frame_processed()

    def finalize(self) -> Dict[str, Dict[str, float]]:
        """
        Calculates the final mean, standard deviation, and fluctuation for the
        concentrations across all processed frames. This method is now idempotent.
        """
        if self._finalized:
            return self.get_result()

        for connectivity, concs in self._raw_concentrations.items():
            if concs:
                self.concentrations[connectivity] = np.mean(concs)
                if len(concs) > 1:
                    self.std[connectivity] = np.std(concs, ddof=1)
                    mean_conc = self.concentrations[connectivity]
                    self.error[connectivity] = self.std[connectivity] / np.sqrt(
                        len(concs)
                    )
                else:
                    self.std[connectivity] = 0.0
                    self.error[connectivity] = 0.0
            else:
                self.concentrations[connectivity] = 0.0
                self.std[connectivity] = 0.0
                self.error[connectivity] = 0.0

            self.std[connectivity] = np.nan_to_num(self.std[connectivity])
            self.error[connectivity] = np.nan_to_num(self.error[connectivity])

        self._finalized = True
        return self.get_result()

    def get_result(self) -> Dict[str, Dict[str, float]]:
        """Returns the finalized analysis results."""
        return {
            "concentrations": self.concentrations,
            "std": self.std,
            "error": self.error,
        }

    def print_to_file(self) -> None:
        """Writes the finalized results to a data file."""
        output = self.finalize()
        self._write_header()
        path = os.path.join(self._settings.export_directory, "concentrations.dat")
        with open(path, "a") as f:
            for connectivity in self.concentrations:
                concentration = output["concentrations"].get(connectivity, 0.0)
                std = output["std"].get(connectivity, 0.0)
                error = output["error"].get(connectivity, 0.0)
                f.write(f"{connectivity},{concentration},{std},{error}\n")
        remove_duplicate_lines(path)

    def _write_header(self) -> None:
        """Initializes the output file with a header."""
        path = os.path.join(self._settings.export_directory, "concentrations.dat")
        number_of_frames = self.frame_processed_count

        if self._settings.analysis.overwrite or not os.path.exists(path):
            mode = "w"
        else:
            if os.path.getsize(path) > 0:
                return
            mode = "a"

        with open(path, mode, encoding="utf-8") as output:
            output.write(f"# Concentration Results\n")
            output.write(f"# Date: {datetime.now()}\n")
            output.write(f"# Frames averaged: {number_of_frames}\n")
            output.write(
                "# Connectivity_type,Concentration,Standard_deviation_ddof=1,Standard_error_ddof=1\n"
            )

    def __str__(self) -> str:
        return f"{self.__class__.__name__}"

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}()"
