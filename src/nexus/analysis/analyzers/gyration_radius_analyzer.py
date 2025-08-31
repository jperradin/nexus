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
    Computes the distribution of gyration radii by cluster size for each connectivity type.

    This analyzer collects the gyration radius of all non-percolating clusters,
    binned by cluster size, and averages them over all processed frames.
    """

    def __init__(self, settings: Settings) -> None:
        """Initializes the analyzer."""
        super().__init__(settings)
        # Raw per-frame storage
        self._raw_gyration_radii: Dict[str, Dict[int, List[float]]] = {}
        self._raw_concentrations: Dict[str, List[float]] = {}

        # Final aggregated results
        self.gyration_radii: Dict[str, Dict[int, float]] = {}
        self.std: Dict[str, Dict[int, float]] = {}
        self.concentrations: Dict[str, float] = {}

        self._finalized: bool = False

    def analyze(self, frame: Frame, connectivities: List[str]) -> None:
        """
        Collects gyration radii of non-percolating clusters for each connectivity type,
        grouped by cluster size.
        """
        clusters = frame.get_clusters()
        concentrations = frame.get_concentration()

        for connectivity in connectivities:
            self._raw_gyration_radii.setdefault(connectivity, {})
            self._raw_concentrations.setdefault(connectivity, [])

            # Get gyration radii grouped by cluster size
            for c in clusters:
                if c.get_connectivity() == connectivity and not c.is_percolating:
                    size = c.get_size()
                    gyr = c.gyration_radius
                    self._raw_gyration_radii[connectivity].setdefault(size, []).append(gyr)

            self._raw_concentrations[connectivity].append(
                concentrations.get(connectivity, 0.0)
            )

        self.update_frame_processed()

    def finalize(self) -> Dict[str, Dict]:
        """
        Computes average gyration radius (and std) per cluster size for each connectivity.
        """
        if self._finalized:
            return self.get_result()

        for connectivity, size_dict in self._raw_gyration_radii.items():
            self.gyration_radii[connectivity] = {}
            self.std[connectivity] = {}
            for size, radii in size_dict.items():
                arr = np.array(radii)
                self.gyration_radii[connectivity][size] = float(np.mean(arr))
                if len(arr) > 1:
                    self.std[connectivity][size] = float(np.std(arr, ddof=1))
                else:
                    self.std[connectivity][size] = 0.0

        for connectivity, concs in self._raw_concentrations.items():
            self.concentrations[connectivity] = np.mean(concs) if concs else 0.0

        self._finalized = True
        return self.get_result()

    def get_result(self) -> Dict[str, Dict]:
        """Returns the finalized results."""
        return {
            "concentrations": self.concentrations,
            "gyration_radii": self.gyration_radii,
            "std": self.std,
        }

    def print_to_file(self) -> None:
        """Writes the finalized results to per-connectivity data files."""
        output = self.finalize()
        self._write_header()

        for connectivity in self.gyration_radii:
            path = os.path.join(
                self._settings.export_directory,
                f"gyration_radius_distribution-{connectivity}.dat",
            )
            with open(path, "a") as f:
                for size, gyr in sorted(
                    output["gyration_radii"][connectivity].items(),
                    key=lambda x: x[0],
                    reverse=True,
                ):
                    std = output["std"][connectivity][size]
                    conc = output["concentrations"].get(connectivity, 0.0)
                    f.write(f"{connectivity},{conc},{size},{gyr},{std}\n")
            remove_duplicate_lines(path)

    def _write_header(self) -> None:
        """Initializes each output file with a header."""
        number_of_frames = self.frame_processed_count
        for connectivity in self._raw_gyration_radii:
            path = os.path.join(
                self._settings.export_directory,
                f"gyration_radius_distribution-{connectivity}.dat",
            )

            if self._settings.analysis.overwrite or not os.path.exists(path):
                mode = "w"
            else:
                if os.path.getsize(path) > 0:
                    continue
                mode = "a"

            with open(path, mode, encoding="utf-8") as output:
                output.write(f"# Gyration Radius Distribution Results\n")
                output.write(f"# Date: {datetime.now()}\n")
                output.write(f"# Frames averaged: {number_of_frames}\n")
                output.write(
                    "# Connectivity_type,Concentration,Cluster_size,Gyration_radius,Standard_deviation_ddof=1\n"
                )

    def __str__(self) -> str:
        return f"{self.__class__.__name__}"

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}()"
