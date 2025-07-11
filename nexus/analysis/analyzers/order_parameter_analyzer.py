from typing import List, Dict
from ...core.frame import Frame
from .base_analyzer import BaseAnalyzer
from ...config.settings import Settings
from ...utils.aesthetics import remove_duplicate_lines

import numpy as np
import os 
from datetime import datetime


class OrderParameterAnalyzer(BaseAnalyzer):
    def __init__(self, settings: Settings) -> None:
        super().__init__(settings)
        self.order_parameters: Dict[str, List[float]] = {}
        self.std: Dict[str, float] = {}
        self.concentrations: Dict[str, List[float]] = {}
    
    def analyze(self, frame: Frame, connectivities: List[str]) -> None:
        clusters = frame.get_clusters()
        concentrations = frame.get_concentration()

        # get all order parameters per connectivity
        for connectivity in connectivities:
            percolating_clusters = [c for c in clusters if c.get_connectivity() == connectivity and c.is_percolating]
            
            if percolating_clusters:
                order_parameters = [c.order_parameter for c in percolating_clusters]
                # order_parameters is a list of 3 floats
                if connectivity not in self.order_parameters:
                    self.order_parameters[connectivity] = []
                self.order_parameters[connectivity].extend(order_parameters)
                self.std[connectivity] = [0.0, 0.0, 0.0]
                if connectivity not in self.concentrations:
                    self.concentrations[connectivity] = []
                self.concentrations[connectivity].append(concentrations[connectivity])
            else:
                if connectivity not in self.order_parameters:
                    self.order_parameters[connectivity] = []
                self.order_parameters[connectivity].extend([[0.0, 0.0, 0.0]])
                self.std[connectivity] = [0.0, 0.0, 0.0]
                if connectivity not in self.concentrations:
                    self.concentrations[connectivity] = []
                if connectivity not in concentrations:
                    self.concentrations[connectivity].append(0.0)
                else:
                    self.concentrations[connectivity].append(concentrations[connectivity])

        self.update_frame_processed(frame)

    def finalize(self) -> None:
        for connectivity, order_parameters in self.order_parameters.items():
            self.order_parameters[connectivity] = np.mean(order_parameters, axis=0)
            if len(order_parameters) == 1:
                self.std[connectivity] = [0.0, 0.0, 0.0]
            else:
                self.std[connectivity] = np.std(order_parameters, ddof=1, axis=0)
            # replace eventual nan with 0.0
            self.std[connectivity] = np.nan_to_num(self.std[connectivity])

        for connectivity, concentrations in self.concentrations.items():
            self.concentrations[connectivity] = np.mean(concentrations)

        return {"concentrations": self.concentrations, "order_parameters": self.order_parameters, "std": self.std}

    def get_result(self) -> Dict[str, float]:
        return {"concentrations": self.concentrations, "order_parameters": self.order_parameters, "std": self.std}

    def update_frame_processed(self, frame: Frame) -> None:
        self.frame_processed.append(frame)

    def print_to_file(self) -> None:
        self._write_header()
        self._write_data()
        

    def get_std(self) -> Dict[str, float]:
        return self.std

    def _write_header(self) -> None:
        """
        Initializes the output file with a header.

        Parameters:
        -----------
            overwrite (bool): Whether to overwrite the existing file.
            path_to_directory (str): The directory where the output file will be saved.
            number_of_frames (int): The number of frames used in averaging.
        """
        path = os.path.join(self._settings.export_directory, "order_parameter.dat")
        number_of_frames = len(self.frame_processed)
        overwrite = self._settings.analysis.overwrite
        if not overwrite and os.path.exists(path):
            with open(path, 'a', encoding='utf-8') as output:
                output.write(f"# Order Parameter Results\n")
                output.write(f"# Date: {datetime.now()}\n")
                output.write(f"# Frames averaged: {number_of_frames}\n")
                output.write("# Connectivity_type,Concentration,Order_parameter,Standard_deviation_ddof=1\n")
            output.close()
        else:
            with open(path, 'w', encoding='utf-8') as output:
                output.write(f"# Order Parameter Results\n")
                output.write(f"# Date: {datetime.now()}\n")
                output.write(f"# Frames averaged: {number_of_frames}\n")
                output.write("# Connectivity_type,Concentration,Order_parameter,Standard_deviation_ddof=1\n")
            output.close()

    def _write_data(self) -> None:
        output = self.finalize()
        path = os.path.join(self._settings.export_directory, "order_parameter.dat")
        with open(path, "a") as f:
            for connectivity in self.order_parameters:
                concentration = output["concentrations"][connectivity]
                order_parameter = output["order_parameters"][connectivity]
                std = output["std"][connectivity]
                f.write(f"{connectivity},{concentration},{order_parameter[0]},{std[0]}\n")
        remove_duplicate_lines(path)

    def __str__(self) -> str:
        return f"{self.__class__.__name__}"

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}()"