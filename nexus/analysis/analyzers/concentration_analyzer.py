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
    Analyzer that computes the concentration of clusters for each connectivity type.
    
    This analyzer tracks the concentration of clusters for each connectivity type
    across all processed frames. The concentration is defined as the ratio of the
    number of clusters of a given connectivity type to the total number of clusters.
    
    Attributes:
        concentrations (Dict): Dictionary mapping connectivity types to their concentrations.
    """
    def __init__(self, settings: Settings) -> None:
        """
        Initialize the ConcentrationAnalyzer.
        
        Args:
            settings (Settings): Configuration settings for the analyzer.
        """
        super().__init__(settings)
        self.concentrations: Dict[str, List[float]] = {}
        self.std: Dict[str, float] = {}
        self.fluctuations: Dict[str, float] = {}

    def analyze(self, frame: Frame, connectivities: List[str]) -> None:
        concentrations = frame.get_concentration()
        for connectivity in connectivities:
            if connectivity not in self.concentrations:
                self.concentrations[connectivity] = []
                self.std[connectivity] = []
                self.fluctuations[connectivity] = []
            if connectivity not in concentrations:
                self.concentrations[connectivity].append(0.0)
                self.std[connectivity].append(0.0)
                self.fluctuations[connectivity].append(0.0)
            else:
                self.concentrations[connectivity].append(concentrations[connectivity])
                self.std[connectivity].append(concentrations[connectivity])
            
        self.update_frame_processed(frame)     

    def update_frame_processed(self, frame: Frame) -> None:
        self.frame_processed.append(frame)   

    def finalize(self) -> Dict[str, float]:
        for connectivity, concentrations in self.concentrations.items():
            if len(concentrations) == 1:
                self.std[connectivity] = 0.0
                self.fluctuations[connectivity] = 0.0
            else:
                self.std[connectivity] = np.std(concentrations, ddof=1)
                self.fluctuations[connectivity] = np.var(concentrations, ddof=1) / np.mean(concentrations)
            self.concentrations[connectivity] = np.mean(concentrations)
            # replace eventual nan with 0.0
            self.std[connectivity] = np.nan_to_num(self.std[connectivity])
        return {"concentrations": self.concentrations, "std": self.std, "fluctuations": self.fluctuations}

    def get_result(self) -> Dict[str, float]:
        return {"concentrations": self.concentrations, "std": self.std, "fluctuations": self.fluctuations}

    def print_to_file(self) -> None:
        self._write_header()
        self._write_data()

    def get_std(self) -> Dict[str, float]:
        return self.std

    def get_fluctuations(self) -> Dict[str, float]:
        return self.fluctuations

    def _write_header(self) -> None:
        path = os.path.join(self._settings.export_directory, "concentrations.dat")
        number_of_frames = len(self.frame_processed)
        overwrite = self._settings.analysis.overwrite
        if not overwrite and os.path.exists(path):
            with open(path, 'a', encoding='utf-8') as output:
                output.write(f"# Concentration Results\n")
                output.write(f"# Date: {datetime.now()}\n")
                output.write(f"# Frames averaged: {number_of_frames}\n")
                output.write("# Connectivity_type,Concentration,Standard_deviation_ddof=1,Fluctuations_ddof=1\n")
            output.close()
        else:
            with open(path, 'w', encoding='utf-8') as output:
                output.write(f"# Concentration Results\n")
                output.write(f"# Date: {datetime.now()}\n")
                output.write(f"# Frames averaged: {number_of_frames}\n")
                output.write("# Connectivity_type,Concentration,Standard_deviation_ddof=1,Fluctuations_ddof=1\n")
            output.close()
        
    def _write_data(self) -> None:
        output = self.finalize()
        path = os.path.join(self._settings.export_directory, "concentrations.dat")
        with open(path, "a") as f:
            for connectivity in self.concentrations:
                concentration = output["concentrations"][connectivity]
                std = output["std"][connectivity]
                fluctuations = output["fluctuations"][connectivity]
                f.write(f"{connectivity},{concentration},{std},{fluctuations}\n")
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