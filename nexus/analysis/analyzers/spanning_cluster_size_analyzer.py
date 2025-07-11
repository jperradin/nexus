from typing import List, Dict
from ...core.frame import Frame
from ...core.cluster import Cluster
from .base_analyzer import BaseAnalyzer
from ...config.settings import Settings
from ...utils.aesthetics import remove_duplicate_lines

import numpy as np
import os
from datetime import datetime


class SpanningClusterSizeAnalyzer(BaseAnalyzer):
    def __init__(self, settings: Settings) -> None:
        super().__init__(settings)
        self.spanning_cluster_sizes: Dict[str, List[float]] = {}
        self.std: Dict[str, float] = {}
        self.concentrations: Dict[str, List[float]] = {}
        self.fluctuations: Dict[str, float] = {}

    def analyze(self, frame: Frame, connectivities: List[str]) -> None:
        clusters = frame.get_clusters()
        concentrations = frame.get_concentration()

        # get all sizes per connectivity
        for connectivity in connectivities:
            sizes = [c.get_size() for c in clusters if c.get_connectivity() == connectivity and not c.is_percolating]
            if sizes:
                max_size = np.max(sizes)

                if connectivity not in self.spanning_cluster_sizes:
                    self.spanning_cluster_sizes[connectivity] = []
                self.spanning_cluster_sizes[connectivity].append(max_size)

                if connectivity not in self.concentrations:
                    self.concentrations[connectivity] = []
                self.concentrations[connectivity].append(concentrations[connectivity])

            else:
                if connectivity not in self.spanning_cluster_sizes:
                    self.spanning_cluster_sizes[connectivity] = []
                self.spanning_cluster_sizes[connectivity].append(0.0)

                if connectivity not in self.concentrations:
                    self.concentrations[connectivity] = []
                if connectivity not in concentrations:
                    self.concentrations[connectivity].append(0.0)
                else:
                    self.concentrations[connectivity].append(concentrations[connectivity])

        self.update_frame_processed(frame)

    def update_frame_processed(self, frame: Frame) -> None:
        self.frame_processed.append(frame)

    def finalize(self) -> Dict[str, float]:
        for connectivity, sizes in self.spanning_cluster_sizes.items():
            self.spanning_cluster_sizes[connectivity] = np.mean(sizes)
            if len(sizes) == 1:
                self.std[connectivity] = 0.0
                self.fluctuations[connectivity] = 0.0
            else:
                self.std[connectivity] = np.std(sizes, ddof=1)
                self.fluctuations[connectivity] = np.var(sizes, ddof=1) / np.mean(sizes)
            # replace eventual nan with 0.0
            self.std[connectivity] = np.nan_to_num(self.std[connectivity])
            self.fluctuations[connectivity] = np.nan_to_num(self.fluctuations[connectivity])

        for connectivity, concentrations in self.concentrations.items():
            self.concentrations[connectivity] = np.mean(concentrations)

        return {"concentrations": self.concentrations, "spanning_cluster_size": self.spanning_cluster_sizes, "std": self.std, "fluctuations": self.fluctuations}

    def get_result(self) -> Dict[str, float]:
        return {"concentrations": self.concentrations, "spanning_cluster_size": self.spanning_cluster_sizes, "std": self.std, "fluctuations": self.fluctuations}

    def print_to_file(self) -> None:
        self._write_header()
        self._write_data()
        

    def get_std(self) -> Dict[str, float]:
        return self.std

    def get_fluctuations(self) -> Dict[str, float]:
        return self.fluctuations

    def _write_header(self) -> None:
        """
        Initializes the output file with a header.

        Parameters:
        -----------
            overwrite (bool): Whether to overwrite the existing file.
            path_to_directory (str): The directory where the output file will be saved.
            number_of_frames (int): The number of frames used in averaging.
        """
        path = os.path.join(self._settings.export_directory, "spanning_cluster_size.dat")
        number_of_frames = len(self.frame_processed)
        overwrite = self._settings.analysis.overwrite
        if not overwrite and os.path.exists(path):
            with open(path, 'a', encoding='utf-8') as output:
                output.write(f"# Spanning Cluster Size Results\n")
                output.write(f"# Date: {datetime.now()}\n")
                output.write(f"# Frames averaged: {number_of_frames}\n")
                output.write("# Connectivity_type,Concentration,Spanning_cluster_size,Standard_deviation_ddof=1,Fluctuations_ddof=1\n")
            output.close()
        else:
            with open(path, 'w', encoding='utf-8') as output:
                output.write(f"# Spanning Cluster Size Results\n")
                output.write(f"# Date: {datetime.now()}\n")
                output.write(f"# Frames averaged: {number_of_frames}\n")
                output.write("# Connectivity_type,Concentration,Spanning_cluster_size,Standard_deviation_ddof=1,Fluctuations_ddof=1\n")
            output.close()

    def _write_data(self) -> None:
        output = self.finalize()
        path = os.path.join(self._settings.export_directory, "spanning_cluster_size.dat")
        with open(path, "a") as f:
            for connectivity in self.spanning_cluster_sizes:
                concentration = output["concentrations"][connectivity]
                spanning_cluster_size = output["spanning_cluster_size"][connectivity]
                std = output["std"][connectivity]
                fluctuations = output["fluctuations"][connectivity]
                f.write(f"{connectivity},{concentration},{spanning_cluster_size},{std},{fluctuations}\n")
        remove_duplicate_lines(path)

    def __str__(self) -> str:
        return f"{self.__class__.__name__}"

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}()"