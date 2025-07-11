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
    Computes the percolation order parameter (P∞) for each connectivity type.

    The order parameter is the fraction of networking nodes that belong to a
    percolating cluster. It is a key metric for identifying the percolation
    threshold in a system.
    """
    def __init__(self, settings: Settings) -> None:
        """Initializes the analyzer."""
        super().__init__(settings)
        # Private attributes to store raw, per-frame data
        self._raw_order_parameters: Dict[str, List[float]] = {}
        self._raw_concentrations: Dict[str, List[float]] = {}

        # Public attributes to hold the final, aggregated results
        self.order_parameters: Dict[str, float] = {}
        self.std: Dict[str, float] = {}
        self.concentrations: Dict[str, float] = {}
        
        # A flag to ensure final calculations are only performed once
        self._finalized: bool = False
    
    def analyze(self, frame: Frame, connectivities: List[str]) -> None:
        """
        Analyzes a single frame to get the order parameter for each connectivity
        type and stores the raw data.
        """
        clusters = frame.get_clusters()
        concentrations = frame.get_concentration()

        for connectivity in connectivities:
            # Initialize lists if this is the first time seeing this connectivity
            self._raw_order_parameters.setdefault(connectivity, [])
            self._raw_concentrations.setdefault(connectivity, [])

            # The order parameter is typically defined for the largest cluster if it percolates.
            # We find the largest percolating cluster.
            percolating_clusters = [c for c in clusters if c.get_connectivity() == connectivity and c.is_percolating]
            
            if percolating_clusters:
                # Assuming the order parameter is associated with the largest of the percolating clusters
                largest_perc_cluster = max(percolating_clusters, key=lambda c: c.size)
                # We are interested in the 1D order parameter (P∞_x) as per the README
                self._raw_order_parameters[connectivity].append(largest_perc_cluster.order_parameter[0])
            else:
                # If no cluster percolates, the order parameter is 0
                self._raw_order_parameters[connectivity].append(0.0)

            self._raw_concentrations[connectivity].append(concentrations.get(connectivity, 0.0))

        self.update_frame_processed(frame)

    def finalize(self) -> Dict[str, Dict[str, float]]:
        """
        Calculates the final mean and standard deviation for the order parameter
        across all processed frames. This method is now idempotent.
        """
        if self._finalized:
            return self.get_result()

        for connectivity, params in self._raw_order_parameters.items():
            if params:
                self.order_parameters[connectivity] = np.mean(params)
                if len(params) > 1:
                    self.std[connectivity] = np.std(params, ddof=1)
                else:
                    self.std[connectivity] = 0.0
            else:
                self.order_parameters[connectivity] = 0.0
                self.std[connectivity] = 0.0
            
            self.std[connectivity] = np.nan_to_num(self.std[connectivity])

        for connectivity, concs in self._raw_concentrations.items():
            self.concentrations[connectivity] = np.mean(concs) if concs else 0.0

        self._finalized = True
        return self.get_result()

    def get_result(self) -> Dict[str, Dict[str, float]]:
        """Returns the finalized analysis results."""
        return {
            "concentrations": self.concentrations, 
            "order_parameters": self.order_parameters, 
            "std": self.std
        }

    def update_frame_processed(self, frame: Frame) -> None:
        self.frame_processed.append(frame)

    def print_to_file(self) -> None:
        """Writes the finalized results to a data file."""
        output = self.finalize()
        self._write_header()
        path = os.path.join(self._settings.export_directory, "order_parameter.dat")
        with open(path, "a") as f:
            for connectivity in self.order_parameters:
                concentration = output["concentrations"].get(connectivity, 0.0)
                order_parameter = output["order_parameters"].get(connectivity, 0.0)
                std = output["std"].get(connectivity, 0.0)
                f.write(f"{connectivity},{concentration},{order_parameter},{std}\n")
        remove_duplicate_lines(path)

    def _write_header(self) -> None:
        """Initializes the output file with a header."""
        path = os.path.join(self._settings.export_directory, "order_parameter.dat")
        number_of_frames = len(self.frame_processed)
        
        if self._settings.analysis.overwrite or not os.path.exists(path):
            mode = 'w'
        else:
            if os.path.getsize(path) > 0: return
            mode = 'a'

        with open(path, mode, encoding='utf-8') as output:
            output.write(f"# Order Parameter Results\n")
            output.write(f"# Date: {datetime.now()}\n")
            output.write(f"# Frames averaged: {number_of_frames}\n")
            output.write("# Connectivity_type,Concentration,Order_parameter,Standard_deviation_ddof=1\n")

    def __str__(self) -> str:
        return f"{self.__class__.__name__}"

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}()"
