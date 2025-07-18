from typing import List, Dict
from ...core.frame import Frame
from .base_analyzer import BaseAnalyzer
from ...config.settings import Settings
from ...utils.aesthetics import remove_duplicate_lines
import numpy as np
import os
from datetime import datetime


class GyrationRadiusAnalyzer(BaseAnalyzer):
    """
    Computes the average gyration radius for each connectivity type.

    This analyzer calculates the mean gyration radius of all non-percolating
    clusters for a given connectivity, averaged over all processed frames.
    This provides a measure of the typical spatial extent of finite clusters.
    """

    def __init__(self, settings: Settings) -> None:
        """Initializes the analyzer."""
        super().__init__(settings)
        # Private attributes to store raw, per-frame data
        self._raw_gyration_radii: Dict[str, List[float]] = {}
        self._raw_concentrations: Dict[str, List[float]] = {}

        # Public attributes to hold the final, aggregated results
        self.gyration_radii: Dict[str, float] = {}
        self.std: Dict[str, float] = {}
        self.concentrations: Dict[str, float] = {}
        self.fluctuations: Dict[str, float] = {}

        # A flag to ensure final calculations are only performed once
        self._finalized: bool = False

    def analyze(self, frame: Frame, connectivities: List[str]) -> None:
        """
        Analyzes a single frame to collect the gyration radius of each
        non-percolating cluster for each connectivity type.
        """
        clusters = frame.get_clusters()
        concentrations = frame.get_concentration()

        for connectivity in connectivities:
            # Initialize lists if this is the first time seeing this connectivity
            self._raw_gyration_radii.setdefault(connectivity, [])
            self._raw_concentrations.setdefault(connectivity, [])

            # Collect gyration radii from all non-percolating clusters of this type
            radii_in_frame = [
                c.gyration_radius
                for c in clusters
                if c.get_connectivity() == connectivity and not c.is_percolating
            ]
            if radii_in_frame:
                self._raw_gyration_radii[connectivity].extend(radii_in_frame)

            self._raw_concentrations[connectivity].append(
                concentrations.get(connectivity, 0.0)
            )

        self.update_frame_processed()

    def finalize(self) -> Dict[str, Dict[str, float]]:
        """
        Calculates the final mean, standard deviation, and fluctuation for the
        gyration radius across all processed frames. This method is now idempotent.
        """
        if self._finalized:
            return self.get_result()

        for connectivity, radii in self._raw_gyration_radii.items():
            if radii:
                self.gyration_radii[connectivity] = np.mean(radii)
                if len(radii) > 1:
                    self.std[connectivity] = np.std(radii, ddof=1)
                    mean_radius = self.gyration_radii[connectivity]
                    self.fluctuations[connectivity] = (
                        np.var(radii, ddof=1) / mean_radius if mean_radius > 0 else 0.0
                    )
                else:
                    self.std[connectivity] = 0.0
                    self.fluctuations[connectivity] = 0.0
            else:
                self.gyration_radii[connectivity] = 0.0
                self.std[connectivity] = 0.0
                self.fluctuations[connectivity] = 0.0

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
            "gyration_radii": self.gyration_radii,
            "std": self.std,
            "fluctuations": self.fluctuations,
        }

    def print_to_file(self) -> None:
        """Writes the finalized results to a single data file."""
        output = self.finalize()
        self._write_header()
        path = os.path.join(self._settings.export_directory, "gyration_radius.dat")
        with open(path, "a") as f:
            for connectivity in self.gyration_radii:
                concentration = output["concentrations"].get(connectivity, 0.0)
                gyration_radius = output["gyration_radii"].get(connectivity, 0.0)
                std = output["std"].get(connectivity, 0.0)
                fluctuations = output["fluctuations"].get(connectivity, 0.0)
                f.write(
                    f"{connectivity},{concentration},{gyration_radius},{std},{fluctuations}\n"
                )
        remove_duplicate_lines(path)

    def _write_header(self) -> None:
        """Initializes the output file with a header."""
        path = os.path.join(self._settings.export_directory, "gyration_radius.dat")
        number_of_frames = self.frame_processed_count

        if self._settings.analysis.overwrite or not os.path.exists(path):
            mode = "w"
        else:
            if os.path.getsize(path) > 0:
                return
            mode = "a"

        with open(path, mode, encoding="utf-8") as output:
            output.write(f"# Average Gyration Radius Results\n")
            output.write(f"# Date: {datetime.now()}\n")
            output.write(f"# Frames averaged: {number_of_frames}\n")
            output.write(
                "# Connectivity_type,Concentration,Average_Gyration_radius,Standard_deviation_ddof=1,Fluctuations_ddof=1\n"
            )

    def __str__(self) -> str:
        return f"{self.__class__.__name__}"

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}()"
