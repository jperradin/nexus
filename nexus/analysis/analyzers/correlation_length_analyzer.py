from typing import List, Dict
from ...core.frame import Frame
from .base_analyzer import BaseAnalyzer
from ...config.settings import Settings
from ...utils.aesthetics import remove_duplicate_lines
import numpy as np
import os
from datetime import datetime

class CorrelationLengthAnalyzer(BaseAnalyzer):
    def __init__(self, settings: Settings) -> None:
        super().__init__(settings)
        self.correlation_length = {}
        self.std = {}
        self.concentrations = {}
        self.fluctuations = {}

    def analyze(self, frame: Frame, connectivities: List[str]) -> None:
        clusters = frame.get_clusters()
        concentrations = frame.get_concentration()
        
        # get all gyration radii per connectivity
        for connectivity in connectivities:
            gyration_radii = [c.gyration_radius for c in clusters if c.get_connectivity() == connectivity and not c.is_percolating]
            sizes = [c.get_size() for c in clusters if c.get_connectivity() == connectivity and not c.is_percolating]
            if gyration_radii:
                unique_sizes, ns = np.unique(sizes, return_counts=True)
                distribution = {}
                std = {}
                for size, gyration_radius in zip(sizes, gyration_radii):
                    if size not in distribution:
                        distribution[size] = []
                        std[size] = []
                    distribution[size].append(gyration_radius)
                    std[size].append(gyration_radius)
                
                A = 0
                B = 0
                for size, gyration_radii in distribution.items():
                    A += 2 * np.mean(gyration_radii)**2 * size**2
                    B += size**2 * len(gyration_radii)
                if B == 0: B = 1 # avoid zero division

                correlation_length = np.sqrt(A / B)
                if connectivity not in self.concentrations:
                    self.concentrations[connectivity] = []
                self.concentrations[connectivity].append(concentrations[connectivity])
                if connectivity not in self.correlation_length:
                    self.correlation_length[connectivity] = []
                self.correlation_length[connectivity].append(correlation_length)
                self.std[connectivity] = std
            else:
                if connectivity not in self.correlation_length:
                    self.correlation_length[connectivity] = []
                self.correlation_length[connectivity].append(0.0)

                if connectivity not in self.concentrations:
                    self.concentrations[connectivity] = []
                if connectivity not in concentrations:
                    self.concentrations[connectivity].append(0.0)
                else:
                    self.concentrations[connectivity].append(concentrations[connectivity])

        self.update_frame_processed(frame)

    def update_frame_processed(self, frame: Frame) -> None:
        self.frame_processed.append(frame)

    def finalize(self) -> None:
        for connectivity, sizes in self.correlation_length.items():
            if len(sizes) == 1:
                self.std[connectivity] = 0.0
                self.fluctuations[connectivity] = 0.0
            else:
                self.std[connectivity] = np.std(sizes, ddof=1)
                self.fluctuations[connectivity] = np.var(sizes, ddof=1) / np.mean(sizes)
            self.correlation_length[connectivity] = np.mean(sizes)
            # replace eventual nan with 0.0
            self.std[connectivity] = np.nan_to_num(self.std[connectivity])

        for connectivity, concentrations in self.concentrations.items():
            self.concentrations[connectivity] = np.mean(concentrations)

        return {"concentrations": self.concentrations, "correlation_length": self.correlation_length, "std": self.std, "fluctuations": self.fluctuations}
                

    def get_result(self) -> Dict[str, float]:
        return {"concentrations": self.concentrations, "correlation_length": self.correlation_length, "std": self.std, "fluctuations": self.fluctuations}

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
        path = os.path.join(self._settings.export_directory, "correlation_length.dat")
        number_of_frames = len(self.frame_processed)
        overwrite = self._settings.analysis.overwrite
        if not overwrite and os.path.exists(path):
            with open(path, 'a', encoding='utf-8') as output:
                output.write(f"# Correlation Length Results\n")
                output.write(f"# Date: {datetime.now()}\n")
                output.write(f"# Frames averaged: {number_of_frames}\n")
                output.write("# Connectivity_type,Concentration,Correlation_length,Standard_deviation_ddof=1,Fluctuations_ddof=1\n")
            output.close()
        else:
            with open(path, 'w', encoding='utf-8') as output:
                output.write(f"# Correlation Length Results\n")
                output.write(f"# Date: {datetime.now()}\n")
                output.write(f"# Frames averaged: {number_of_frames}\n")
                output.write("# Connectivity_type,Concentration,Correlation_length,Standard_deviation_ddof=1,Fluctuations_ddof=1\n")
            output.close()

    def _write_data(self) -> None:
        output = self.finalize()
        path = os.path.join(self._settings.export_directory, "correlation_length.dat")
        with open(path, "a") as f:
            for connectivity in self.correlation_length:
                concentration = output["concentrations"][connectivity]
                correlation_length = output["correlation_length"][connectivity]
                std = output["std"][connectivity]
                fluctuations = output["fluctuations"][connectivity]
                f.write(f"{connectivity},{concentration},{correlation_length},{std},{fluctuations}\n")
        remove_duplicate_lines(path)

    def __str__(self) -> str:
        return f"{self.__class__.__name__}"

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}()"