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
    Computes the distribution of cluster sizes, n(s), for each connectivity type.

    This analyzer tracks how many clusters of each size exist for each connectivity
    type across all processed frames, which is fundamental for percolation theory.
    It excludes percolating clusters from the analysis to focus on the finite clusters.
    """
    def __init__(self, settings: Settings) -> None:
        """Initializes the analyzer."""
        super().__init__(settings)
        # Private attributes to store raw, per-frame data
        self._raw_size_distribution: Dict[str, Dict[int, List[int]]] = {}
        self._raw_concentrations: Dict[str, List[float]] = {}

        # Public attributes to hold the final, aggregated results
        self.size_distribution: Dict[str, Dict[int, float]] = {}
        self.std: Dict[str, Dict[int, float]] = {}
        self.concentrations: Dict[str, float] = {}
        
        # A flag to ensure final calculations are only performed once
        self._finalized: bool = False

    def analyze(self, frame: Frame, connectivities: List[str]) -> None:
        """
        Analyzes a single frame to compute the cluster size distribution for each
        connectivity type and stores the raw data.
        """
        clusters = frame.get_clusters()
        concentrations = frame.get_concentration()

        for connectivity in connectivities:
            # Initialize dictionaries if this is the first time seeing this connectivity
            self._raw_size_distribution.setdefault(connectivity, {})
            self._raw_concentrations.setdefault(connectivity, [])

            sizes = [c.get_size() for c in clusters if c.get_connectivity() == connectivity and not c.is_percolating]
            
            if sizes:
                unique_sizes, ns = np.unique(sizes, return_counts=True)
                for s, n in zip(unique_sizes, ns):
                    self._raw_size_distribution[connectivity].setdefault(s, []).append(n)
            
            # Record concentration for this frame
            self._raw_concentrations[connectivity].append(concentrations.get(connectivity, 0.0))

        self.update_frame_processed(frame)

    def update_frame_processed(self, frame: Frame) -> None:
        self.frame_processed.append(frame)

    def finalize(self) -> Dict:
        """
        Calculates the final mean and standard deviation for the cluster size
        distribution across all processed frames. This method is now idempotent.
        """
        if self._finalized:
            return self.get_result()

        for connectivity, size_data in self._raw_size_distribution.items():
            self.size_distribution.setdefault(connectivity, {})
            self.std.setdefault(connectivity, {})
            for size, counts in size_data.items():
                # The final value is the total count divided by the number of frames
                total_count = np.sum(counts)
                num_frames = len(self.frame_processed)
                self.size_distribution[connectivity][size] = total_count / num_frames if num_frames > 0 else 0.0
                
                # To calculate std dev, we need to account for frames where a size didn't appear
                all_counts_for_size = counts + [0] * (num_frames - len(counts))
                if len(all_counts_for_size) > 1:
                    self.std[connectivity][size] = np.std(all_counts_for_size, ddof=1)
                else:
                    self.std[connectivity][size] = 0.0

        for connectivity, concs in self._raw_concentrations.items():
            self.concentrations[connectivity] = np.mean(concs) if concs else 0.0

        self._finalized = True
        return self.get_result()

    def get_result(self) -> Dict[str, Dict]:
        """Returns the finalized analysis results."""
        return {
            "concentrations": self.concentrations,
            "size_distribution": self.size_distribution,
            "std": self.std
        }

    def print_to_file(self) -> None:
        """Writes the finalized results to data files, one per connectivity."""
        output = self.finalize()
        
        for connectivity in self.size_distribution:
            self._write_header(connectivity)
            path = os.path.join(self._settings.export_directory, f"cluster_size_distribution-{connectivity}.dat")
            
            # Sort by size in descending order for plotting
            sorted_sizes = sorted(self.size_distribution[connectivity].keys(), reverse=True)
            
            with open(path, "a") as f:
                for size in sorted_sizes:
                    concentration = output["concentrations"].get(connectivity, 0.0)
                    n_s = output["size_distribution"][connectivity].get(size, 0.0)
                    std_dev = output["std"][connectivity].get(size, 0.0)
                    f.write(f"{connectivity},{concentration},{size},{n_s},{std_dev}\n")
            remove_duplicate_lines(path)

    def _write_header(self, connectivity: str) -> None:
        """Initializes the output file with a header for a given connectivity."""
        path = os.path.join(self._settings.export_directory, f"cluster_size_distribution-{connectivity}.dat")
        number_of_frames = len(self.frame_processed)
        
        if self._settings.analysis.overwrite or not os.path.exists(path):
            mode = 'w'
        else:
            if os.path.getsize(path) > 0: return
            mode = 'a'

        with open(path, mode, encoding='utf-8') as output:
            output.write(f"# Cluster Size Distribution Results for {connectivity}\n")
            output.write(f"# Date: {datetime.now()}\n")
            output.write(f"# Frames averaged: {number_of_frames}\n")
            output.write("# Connectivity_type,Concentration,Cluster_size,N_clusters_per_frame,Standard_deviation_ddof=1\n")

    def __str__(self) -> str:
        return f"{self.__class__.__name__}"

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}()"
