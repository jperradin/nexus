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
    Analyzer that computes the average cluster size for each connectivity type.

    The average cluster size <S> is calculated as the weighted mean of cluster sizes,
    where each cluster size is weighted by the number of nodes in that cluster size:

    <S> = Σ(s² * n(s)) / Σ(s * n(s))

    where s is the cluster size, n(s) is the number of clusters of size s.
    This is also known as the weight-average cluster size and is commonly
    used in percolation theory to characterize cluster size distributions.

    Attributes:
        average_sizes (Dict[str, List[float]]): Dictionary mapping connectivity types to
            lists of average cluster sizes across frames.
        std (Dict[str, float]): Dictionary mapping connectivity types to
            standard deviations of average cluster sizes.
        concentrations (Dict[str, List[float]]): Dictionary mapping connectivity types to
            lists of concentration values across frames.
    """

    def __init__(self, settings: Settings) -> None:
        """
        Initialize the AverageClusterSizeAnalyzer.

        Args:
            settings (Settings): Configuration settings for the analyzer.
        """
        super().__init__(settings)
        self.average_sizes: Dict[str, List[float]] = {}
        self.std: Dict[str, float] = {}
        self.fluctuations: Dict[str, float] = {}
        self.concentrations: Dict[str, List[float]] = {}

    def analyze(self, frame: Frame, connectivities: List[str]) -> None:
        """
        Analyze a frame to compute the average cluster size for each connectivity type.

        This method calculates the average cluster size <S> using the formula:
        <S> = Σ(s² * n(s)) / Σ(s * n(s))
        where s is the cluster size and n(s) is the number of clusters of size s.

        Only non-percolating clusters are considered in this calculation.

        Args:
            frame (Frame): The frame to analyze.
        """
        clusters = frame.get_clusters()
        concentrations = frame.get_concentration()

        # get all sizes per connectivity
        for connectivity in connectivities:
            sizes = [
                c.get_size()
                for c in clusters
                if c.get_connectivity() == connectivity and not c.is_percolating
            ]
            if sizes:
                sizes, ns = np.unique(sizes, return_counts=True)

                # calculate average cluster size
                average_size = 0.0
                for s, n in zip(sizes, ns):
                    average_size += (s**2 * n) / np.sum(sizes * ns)

                if connectivity not in self.average_sizes:
                    self.average_sizes[connectivity] = []
                self.average_sizes[connectivity].append(average_size)

                if connectivity not in self.concentrations:
                    self.concentrations[connectivity] = []
                self.concentrations[connectivity].append(concentrations[connectivity])

            else:
                if connectivity not in self.average_sizes:
                    self.average_sizes[connectivity] = []
                self.average_sizes[connectivity].append(0.0)

                if connectivity not in self.concentrations:
                    self.concentrations[connectivity] = []
                if connectivity not in concentrations:
                    self.concentrations[connectivity].append(0.0)
                else:
                    self.concentrations[connectivity].append(
                        concentrations[connectivity]
                    )

        # update frame processed
        self.update_frame_processed(frame)

    def update_frame_processed(self, frame: Frame) -> None:
        """
        Update the list of processed frames.

        Args:
            frame (Frame): The frame that has been processed.
        """
        self.frame_processed.append(frame)

    def finalize(self) -> Dict[str, float]:
        """
        Finalize the analysis by calculating means and standard deviations.

        This method computes final average values across all frames for each connectivity type.

        Returns:
            Dict[str, float]: Dictionary containing results with keys:
                - 'concentrations': Mean concentrations for each connectivity
                - 'average_cluster_size': Mean average cluster sizes for each connectivity
                - 'std': Standard deviations of average cluster sizes for each connectivity
        """
        for connectivity, sizes in self.average_sizes.items():
            if len(sizes) == 1:
                self.std[connectivity] = 0.0
                self.fluctuations[connectivity] = 0.0
            else:
                self.std[connectivity] = np.std(sizes, ddof=1)
                self.fluctuations[connectivity] = np.var(sizes, ddof=1) / np.mean(sizes)
            self.average_sizes[connectivity] = np.mean(sizes)
            # replace eventual nan with 0.0
            self.std[connectivity] = np.nan_to_num(self.std[connectivity])

        for connectivity, concentrations in self.concentrations.items():
            self.concentrations[connectivity] = np.mean(concentrations)

        return {
            "concentrations": self.concentrations,
            "average_cluster_size": self.average_sizes,
            "std": self.std,
            "fluctuations": self.fluctuations,
        }

    def get_result(self) -> Dict[str, float]:
        """
        Get the current analysis results.

        Returns:
            Dict[str, float]: Dictionary containing results with keys:
                - 'concentrations': Mean concentrations for each connectivity
                - 'average_cluster_size': Average cluster sizes for each connectivity
                - 'std': Standard deviations of average cluster sizes for each connectivity
        """
        return {
            "concentrations": self.concentrations,
            "average_cluster_size": self.average_sizes,
            "std": self.std,
            "fluctuations": self.fluctuations,
        }

    def print_to_file(self) -> None:
        """
        Write the analysis results to a file.

        This method writes the average cluster size results to a file named
        'average_cluster_size.dat' in the export directory specified in settings.
        """
        self._write_header()
        self._write_data()

    def get_std(self) -> Dict[str, float]:
        """
        Get the standard deviations of average cluster sizes.

        Returns:
            Dict[str, float]: Dictionary mapping connectivity types to standard deviations.
        """
        return self.std

    def get_fluctuations(self) -> Dict[str, float]:
        """
        Get the fluctuations of average cluster sizes.

        Returns:
            Dict[str, float]: Dictionary mapping connectivity types to fluctuations.
        """
        return self.fluctuations

    def _write_header(self) -> None:
        """
        Initialize the output file with a header.

        This method creates or appends to the output file and writes header information
        including date, number of frames processed, and column descriptions.

        The file will be created at the path specified in settings.export_directory.
        If settings.analysis.overwrite is False and the file already exists,
        data will be appended; otherwise, the file will be overwritten.
        """
        path = os.path.join(self._settings.export_directory, "average_cluster_size.dat")
        number_of_frames = len(self.frame_processed)
        overwrite = self._settings.analysis.overwrite
        if not overwrite and os.path.exists(path):
            with open(path, "a", encoding="utf-8") as output:
                output.write(f"# Average Cluster Size Results\n")
                output.write(f"# Date: {datetime.now()}\n")
                output.write(f"# Frames averaged: {number_of_frames}\n")
                output.write(
                    "# Connectivity_type,Concentration,Average_cluster_size,Standard_deviation_ddof=1,Fluctuations_ddof=1\n"
                )
            output.close()
        else:
            with open(path, "w", encoding="utf-8") as output:
                output.write(f"# Average Cluster Size Results\n")
                output.write(f"# Date: {datetime.now()}\n")
                output.write(f"# Frames averaged: {number_of_frames}\n")
                output.write(
                    "# Connectivity_type,Concentration,Average_cluster_size,Standard_deviation_ddof=1,Fluctuations_ddof=1\n"
                )
            output.close()

    def _write_data(self) -> None:
        """
        Write the analysis data to the output file.

        This method appends the analysis results (connectivity types, concentrations,
        average cluster sizes, and standard deviations) to the output file, and
        removes any duplicate lines that might have been introduced.
        """
        output = self.finalize()
        path = os.path.join(self._settings.export_directory, "average_cluster_size.dat")
        with open(path, "a") as f:
            for connectivity in self.average_sizes:
                concentration = output["concentrations"][connectivity]
                average_size = output["average_cluster_size"][connectivity]
                std = output["std"][connectivity]
                fluctuations = output["fluctuations"][connectivity]
                f.write(
                    f"{connectivity},{concentration},{average_size},{std},{fluctuations}\n"
                )
        remove_duplicate_lines(path)

    def __str__(self) -> str:
        """
        Return a string representation of the analyzer.

        Returns:
            str: The name of the analyzer class.
        """
        return f"{self.__class__.__name__}"

    def __repr__(self) -> str:
        """
        Return a string representation for debugging.

        Returns:
            str: A string representation that could be used to recreate the analyzer.
        """
        return f"{self.__class__.__name__}()"

