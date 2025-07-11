from typing import List, Dict
from ...core.frame import Frame
from .base_analyzer import BaseAnalyzer
from ...config.settings import Settings
from ...utils.aesthetics import remove_duplicate_lines
import numpy as np
import os
from datetime import datetime

class ClusterSizeDistributionAnalyzer(BaseAnalyzer):
    """
    Analyzer that computes the distribution of cluster sizes for each connectivity type.
    
    This analyzer tracks how many clusters of each size exist for each connectivity
    type across all processed frames. It excludes percolating clusters from the analysis.
    The distribution data can be used to understand the frequency of different cluster
    sizes in the system and how they vary with connectivity.
    
    Attributes:
        size_distribution (Dict): Nested dictionary mapping connectivity types to dictionaries
            that map cluster sizes to counts of clusters of that size.
        std (Dict): Nested dictionary mapping connectivity types to dictionaries
            that map cluster sizes to standard deviations of cluster counts across frames.
        concentrations (Dict): Dictionary mapping connectivity types to their concentrations.
    """
    def __init__(self, settings: Settings) -> None:
        """
        Initialize the ClusterSizeDistributionAnalyzer.
        
        Args:
            settings (Settings): Configuration settings for the analyzer.
        """
        super().__init__(settings)
        self.size_distribution = {}
        self.std = {}
        self.concentrations = {}

    def analyze(self, frame: Frame, connectivities: List[str]) -> None:
        """
        Analyze a frame to compute the cluster size distribution for each connectivity type.
        
        This method extracts clusters from the frame, groups them by connectivity type,
        and counts the number of clusters of each size. Only non-percolating clusters
        are considered in this analysis.
        
        Args:
            frame (Frame): The frame to analyze.
            connectivities (List[str]): List of connectivities to analyze.
        """
        clusters = frame.get_clusters()
        concentrations = frame.get_concentration()

        # get all sizes per connectivity
        for connectivity in connectivities:
            sizes = [c.get_size() for c in clusters if c.get_connectivity() == connectivity and not c.is_percolating]
            sizes, ns = np.unique(sizes, return_counts=True)
            if sizes.any():
                for s in sizes:
                    if connectivity not in self.size_distribution:
                        self.size_distribution[connectivity] = {}
                        self.std[connectivity] = {}
                    if connectivity not in self.concentrations:
                        self.concentrations[connectivity] = []
                    if s not in self.size_distribution[connectivity]:
                        self.size_distribution[connectivity][s] = []
                        self.std[connectivity][s] = []
                    self.size_distribution[connectivity][s].append(ns[sizes==s])
                    self.std[connectivity][s].append(ns[sizes==s])
                    self.concentrations[connectivity].append(concentrations[connectivity])
            else:
                key = np.int64(0)
                value = np.array([0])
                if connectivity not in self.size_distribution:
                    self.size_distribution[connectivity] = {key: [value]}
                    self.std[connectivity] = {key: [value]}
                    if connectivity not in self.concentrations:
                        self.concentrations[connectivity] = []
                    if connectivity not in concentrations:
                        self.concentrations[connectivity].append(0.0)
                    else:
                        self.concentrations[connectivity].append(concentrations[connectivity])
                else:
                    if key not in self.size_distribution[connectivity]:
                        self.size_distribution[connectivity][key] = [value]
                        self.std[connectivity][key] = [value]
                    else:
                        self.size_distribution[connectivity][key].append(value)
                        self.std[connectivity][key].append(value)
                    if connectivity not in self.concentrations:
                        self.concentrations[connectivity] = []
                    if connectivity not in concentrations:
                        self.concentrations[connectivity].append(0.0)
                    else:
                        self.concentrations[connectivity].append(concentrations[connectivity])

        self.update_frame_processed(frame)

    def update_frame_processed(self, frame: Frame) -> None:
        """
        Update the list of processed frames.
        
        Args:
            frame (Frame): The frame that has been processed.
        """
        self.frame_processed.append(frame)

    def finalize(self) -> Dict:
        """
        Finalize the analysis by calculating sums and standard deviations.
        
        This method computes final values across all frames for each connectivity
        type and cluster size, including:
        - The total count of clusters of each size for each connectivity type
        - The standard deviation of cluster counts across frames
        
        Returns:
            Dict: Dictionary containing results with keys:
                - 'concentrations': Concentrations for each connectivity
                - 'size_distribution': Total counts of clusters by size for each connectivity
                - 'std': Standard deviations of cluster counts for each size and connectivity
        """
        for connectivity, sizes in self.size_distribution.items():
            for size, ns in sizes.items():
                self.size_distribution[connectivity][size] = np.sum(ns)
                if len(ns) == 1:
                    self.std[connectivity][size] = 0.0
                else:
                    self.std[connectivity][size] = np.std(ns, ddof=1)
        
        if self.concentrations:
            for connectivity in self.concentrations:
                self.concentrations[connectivity] = np.mean(self.concentrations[connectivity])

        return {"concentrations": self.concentrations, "size_distribution": self.size_distribution, "std": self.std}

    def get_result(self) -> Dict[str, float]:
        """
        Get the current analysis results.
        
        Returns:
            Dict: Dictionary containing results with keys:
                - 'concentrations': Concentrations for each connectivity
                - 'size_distribution': Counts of clusters by size for each connectivity
                - 'std': Standard deviations of cluster counts for each size and connectivity
        """
        return {"concentrations": self.concentrations, "size_distribution": self.size_distribution, "std": self.std}

    def print_to_file(self) -> None:
        """
        Write the analysis results to a file.
        
        This method writes the cluster size distribution results to a file for
        each connectivity type in the export directory specified in settings.
        """
        self._write_header()
        self._write_data()

    def _write_header(self) -> None:
        """
        Initialize the output files with headers.
        
        This method creates or appends to output files for each connectivity type and writes
        header information including date, number of frames processed, and column descriptions.
        
        The files will be created in the directory specified in settings.export_directory.
        If settings.analysis.overwrite is False and a file already exists, data will be appended;
        otherwise, the file will be overwritten.
        """
        for connectivity in self.size_distribution:
            path = os.path.join(self._settings.export_directory, f"cluster_size_distribution-{connectivity}.dat")
            number_of_frames = len(self.frame_processed)
            overwrite = self._settings.analysis.overwrite
            if not overwrite and os.path.exists(path):
                with open(path, 'a', encoding='utf-8') as output:
                    output.write(f"# Cluster Size Distribution \u279c {number_of_frames} frames averaged.\n")
                    output.write(f"# Date: {datetime.now()}\n")
                    output.write(f"# Frames averaged: {number_of_frames}\n")
                    output.write("# Connectivity_type,Concentration,Cluster_size,N_clusters,Standard_deviation_ddof=1\n")
                output.close()
            else:
                with open(path, 'w', encoding='utf-8') as output:
                    output.write(f"# Cluster Size Distribution \u279c {number_of_frames} frames averaged.\n")
                    output.write(f"# Date: {datetime.now()}\n")
                    output.write(f"# Frames averaged: {number_of_frames}\n")
                    output.write("# Connectivity_type,Concentration,Cluster_size,N_clusters,Standard_deviation_ddof=1\n")
                output.close()

    def _write_data(self) -> None:
        """
        Write the analysis data to the output files.
        
        This method appends the analysis results (connectivity types, concentrations,
        cluster sizes, cluster counts, and standard deviations) to the output files,
        sorted by descending cluster size, and removes any duplicate lines that might
        have been introduced.
        """
        output = self.finalize()
        # sort cluster sizes by descending size
        for connectivity in output["size_distribution"]:
            path = os.path.join(self._settings.export_directory, f"cluster_size_distribution-{connectivity}.dat")
            with open(path, "a") as f:
                for size, ns in sorted(output["size_distribution"][connectivity].items(), key=lambda x: x[0], reverse=True):
                    std = output["std"][connectivity][size]
                    f.write(f"{connectivity},{output['concentrations'][connectivity]},{size},{ns},{std}\n")
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