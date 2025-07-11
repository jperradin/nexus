from typing import List, Dict, TextIO
import os
from ...core.cluster import Cluster
from .base_writer import BaseWriter
from ...config.settings import Settings

class ClustersWriter(BaseWriter):
    def __init__(self, settings: Settings) -> None:
        super().__init__(settings)
        self._settings: Settings = settings

    def set_clusters(self, clusters: List[Cluster]) -> None:
        # Sort clusters by size
        self._clusters: List[Cluster] = sorted(clusters, key=lambda cluster: cluster.size, reverse=True)
    
    def write(self) -> None:
        # check options
        if self._settings.clustering.print_mode == "none":
            return
        elif self._settings.clustering.print_mode == "all":
            self._write_all()
        elif self._settings.clustering.print_mode == "connectivity":
            self._write_connectivity()
        elif self._settings.clustering.print_mode == "individual":
            self._write_individual()

    def _write_all(self) -> None:
        path = os.path.join(self._settings.export_directory, "unwrapped_clusters")
        if not os.path.exists(path):
            os.makedirs(path)
        frame_id = self._clusters[0].frame_id
        path = os.path.join(path, f'all_unwrapped_clusters-frame_{frame_id}.xyz')
        with open(path, 'w') as f:
            for cluster in self._clusters:
                self._write_header(f, cluster)
                self._write_cluster(f, cluster)
        f.close()

    def _write_connectivity(self) -> None:
        clusters_per_connectivity: Dict[str, List[Cluster]] = {}
        for cluster in self._clusters:
            connectivity = cluster.connectivity
            if connectivity not in clusters_per_connectivity:
                clusters_per_connectivity[connectivity] = []
            clusters_per_connectivity[connectivity].append(cluster)

        for connectivity in clusters_per_connectivity.keys():
            path = os.path.join(self._settings.export_directory, "unwrapped_clusters")
            if not os.path.exists(path):
                os.makedirs(path)
            frame_id = clusters_per_connectivity[connectivity][0].frame_id
            path = os.path.join(path, f'{connectivity}_unwrapped_clusters-frame_{frame_id}.xyz')
            with open(path, 'w') as f:
                for cluster in clusters_per_connectivity[connectivity]:
                    self._write_header(f, cluster)
                    self._write_cluster(f, cluster)
            f.close()
            
    def _write_individual(self) -> None:
        path = os.path.join(self._settings.export_directory, "unwrapped_clusters")
        if not os.path.exists(path):
            os.makedirs(path)
        frame_id = self._clusters[0].frame_id
        for cluster in self._clusters:
            path = os.path.join(path, f'individual_unwrapped_clusters-frame_{frame_id}-cluster_{cluster.root_id}.xyz')
            with open(path, 'w') as f:
                self._write_header(f, cluster)
                self._write_cluster(f, cluster)
            f.close()

    def _write_header(self, f: TextIO, cluster: Cluster) -> None:
        f.write(f'{cluster.size}\n')
        lxx, lxy, lxz = cluster.lattice[0]
        lyx, lyy, lyz = cluster.lattice[1]
        lzx, lzy, lzz = cluster.lattice[2]
        lattice_line = f"Lattice=\"{lxx} {lxy} {lxz} {lyx} {lyy} {lyz} {lzx} {lzy} {lzz}\" Properties=species:S:1:index:I:1:pos:R:3:cluster_id:I:1:order_parameter:R:1:center_of_mass:R:3:gyration_radius:R:1"
        f.write(lattice_line + '\n')

    def _write_cluster(self, f: TextIO, cluster: Cluster) -> None:
        for symbol, index, position in zip(cluster.symbols, cluster.indices, cluster.unwrapped_positions):
            f.write(f'{symbol} {index} {position[0]:.5f} {position[1]:.5f} {position[2]:.5f} {cluster.root_id} {cluster.order_parameter[0]:5.5f} {cluster.center_of_mass[0]:5.5f} {cluster.center_of_mass[1]:5.5f} {cluster.center_of_mass[2]:5.5f} {cluster.gyration_radius}\n')