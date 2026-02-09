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
    Computes the mean gyration radius binned by cluster size for each connectivity.

    Collects the gyration radius of all non-percolating clusters, groups them
    by cluster size, and averages over all processed frames. Results are
    written to one file per connectivity.

    Attributes:
        _raw_gyration_radii (Dict[str, Dict[int, List[float]]]): Per-cluster
            gyration radii keyed by connectivity then cluster size.
        _raw_concentrations (Dict[str, List[float]]): Per-frame concentrations.
        gyration_radii (Dict[str, Dict[int, float]]): Mean gyration radius per
            cluster size per connectivity.
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
        # Raw per-frame storage
        self._raw_gyration_radii: Dict[str, Dict[int, List[float]]] = {}
        self._raw_concentrations: Dict[str, List[float]] = {}

        # Final aggregated results
        self.gyration_radii: Dict[str, Dict[int, float | np.float64]] = {}
        self.std: Dict[str, Dict[int, float | np.float64]] = {}
        self.concentrations: Dict[str, float | np.float64] = {}

        self._finalized: bool = False

    def analyze(self, frame: Frame, connectivities: List[str]) -> None:
        """
        Collect gyration radii of non-percolating clusters, grouped by size.

        Args:
            frame (Frame): The frame to analyze.
            connectivities (List[str]): Connectivity labels to analyze.
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
                    self._raw_gyration_radii[connectivity].setdefault(size, []).append(
                        gyr
                    )

            self._raw_concentrations[connectivity].append(
                concentrations.get(connectivity, 0.0)
            )

        self.update_frame_processed()

    def finalize(self) -> Dict[str, Dict]:
        """
        Compute mean gyration radius and standard deviation per cluster size.

        Returns:
            Dict[str, Dict]: The finalized results dictionary.
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
        """
        Return the current results dictionary.

        Returns:
            Dict[str, Dict]: Keys are ``"concentrations"``,
                ``"gyration_radii"``, and ``"std"``.
        """
        return {
            "concentrations": self.concentrations,
            "gyration_radii": self.gyration_radii,
            "std": self.std,
        }

    def print_to_file(self) -> None:
        """Write results to one ``gyration_radius_distribution-<conn>.dat`` per connectivity."""
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
        """Write the CSV header to each per-connectivity output file if needed."""
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
                output.write("# Gyration Radius Distribution Results\n")
                output.write(f"# Date: {datetime.now()}\n")
                output.write(f"# Frames averaged: {number_of_frames}\n")
                output.write(
                    "# Connectivity_type,Concentration,Cluster_size,Gyration_radius,Standard_deviation_ddof=1\n"
                )

    def __str__(self) -> str:
        """Return the class name."""
        return f"{self.__class__.__name__}"

    def __repr__(self) -> str:
        """Return a reproducible string representation."""
        return f"{self.__class__.__name__}()"
