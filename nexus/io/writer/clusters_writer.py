from typing import List, Dict, TextIO
import os
from ...core.cluster import Cluster
from .base_writer import BaseWriter
from ...config.settings import Settings

class ClustersWriter(BaseWriter):
    """
    Writes cluster data to files, including atomic positions and connectivity.

    For each cluster, this writer generates an XYZ file with the unwrapped atomic
    coordinates of both the primary networking atoms and any associated 'decorating'
    atoms (e.g., bridging oxygens). If linkages are defined, it also generates a
    corresponding '.bonds' file.
    """
    def __init__(self, settings: Settings) -> None:
        super().__init__(settings)
        self._settings: Settings = settings

    def set_clusters(self, clusters: List[Cluster]) -> None:
        self._clusters: List[Cluster] = sorted(clusters, key=lambda cluster: cluster.size, reverse=True)
    
    def write(self) -> None:
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
        if not os.path.exists(path): os.makedirs(path)
        if not self._clusters: return
        
        frame_id = self._clusters[0].frame_id
        xyz_path = os.path.join(path, f'all_unwrapped_clusters-frame_{frame_id}.xyz')
        bonds_path = os.path.join(path, f'all_unwrapped_clusters-frame_{frame_id}.bonds')
        
        # Calculate total atoms for the header
        total_atoms = sum(cluster.size + len(cluster.decoration_atoms) for cluster in self._clusters)
        
        with open(xyz_path, 'w') as xyz_file, open(bonds_path, 'w') as bonds_file:
            xyz_file.write(f"{total_atoms}\n")
            self._write_header_comment(xyz_file, self._clusters[0])

            global_node_id_to_local_index = {}
            current_local_index = 1
            
            # Build the complete ID-to-index map for all atoms (networking and decorating)
            for cluster in self._clusters:
                for node_id in cluster.indices:
                    global_node_id_to_local_index[node_id] = current_local_index
                    current_local_index += 1
                # Decorating atoms are also part of the atom count
                current_local_index += len(cluster.decoration_atoms)

            # Reset and write atom positions and bond information
            current_local_index = 1
            for cluster in self._clusters:
                # Write atoms and update the map on the fly
                cluster_map = self._write_cluster_atoms(xyz_file, cluster, current_local_index)
                current_local_index += len(cluster_map)
                
                # Write bonds using the map for this cluster
                self._write_cluster_bonds(bonds_file, cluster, cluster_map)

    def _write_connectivity(self) -> None:
        clusters_per_connectivity: Dict[str, List[Cluster]] = {}
        for cluster in self._clusters:
            clusters_per_connectivity.setdefault(cluster.connectivity, []).append(cluster)

        for connectivity, cluster_list in clusters_per_connectivity.items():
            path = os.path.join(self._settings.export_directory, "unwrapped_clusters", connectivity)
            if not os.path.exists(path): os.makedirs(path)
            if not cluster_list: continue
                
            frame_id = cluster_list[0].frame_id
            xyz_path = os.path.join(path, f'{connectivity}_unwrapped_clusters-frame_{frame_id}.xyz')
            bonds_path = os.path.join(path, f'{connectivity}_unwrapped_clusters-frame_{frame_id}.bonds')
            
            total_atoms = sum(c.size + len(c.decoration_atoms) for c in cluster_list)

            with open(xyz_path, 'w') as xyz_file, open(bonds_path, 'w') as bonds_file:
                xyz_file.write(f"{total_atoms}\n")
                self._write_header_comment(xyz_file, cluster_list[0])

                current_local_index = 1
                for cluster in cluster_list:
                    cluster_map = self._write_cluster_atoms(xyz_file, cluster, current_local_index)
                    current_local_index += len(cluster_map)
                    self._write_cluster_bonds(bonds_file, cluster, cluster_map)
            
    def _write_individual(self) -> None:
        path = os.path.join(self._settings.export_directory, "unwrapped_clusters")
        if not os.path.exists(path): os.makedirs(path)
        if not self._clusters: return
        
        frame_id = self._clusters[0].frame_id
        for cluster in self._clusters:
            xyz_path = os.path.join(path, f'cluster-frame_{frame_id}-id_{cluster.root_id}.xyz')
            bonds_path = os.path.join(path, f'cluster-frame_{frame_id}-id_{cluster.root_id}.bonds')
            
            total_atoms = cluster.size + len(cluster.decoration_atoms)
            with open(xyz_path, 'w') as xyz_file:
                xyz_file.write(f"{total_atoms}\n")
                self._write_header_comment(xyz_file, cluster)
                # For individual files, the local index starts at 1
                node_id_to_local_index = self._write_cluster_atoms(xyz_file, cluster, 1)

            with open(bonds_path, 'w') as bonds_file:
                self._write_cluster_bonds(bonds_file, cluster, node_id_to_local_index)

    def _write_header_comment(self, f: TextIO, cluster: Cluster) -> None:
        """Writes only the comment/properties line of the XYZ header."""
        lxx, lxy, lxz = cluster.lattice[0]
        lyx, lyy, lyz = cluster.lattice[1]
        lzx, lzy, lzz = cluster.lattice[2]
        lattice_line = f'Lattice="{lxx} {lxy} {lxz} {lyx} {lyy} {lyz} {lzx} {lzy} {lzz}" Properties=species:S:1:index:I:1:pos:R:3:cluster_id:I:1\n'
        f.write(lattice_line)

    def _write_cluster_atoms(self, f: TextIO, cluster: Cluster, start_index: int) -> Dict[int, int]:
        """
        Writes all atoms (networking and decorating) to an already open file.
        Returns a map of global node ID to the local index within the file.
        """
        node_id_to_local_index = {}
        local_index = start_index
        
        # Write primary networking atoms
        for symbol, global_id, position in zip(cluster.symbols, cluster.indices, cluster.unwrapped_positions):
            f.write(f'{symbol} {global_id} {position[0]:.5f} {position[1]:.5f} {position[2]:.5f} {cluster.root_id}\n')
            node_id_to_local_index[global_id] = local_index
            local_index += 1
            
        # Write decorating atoms
        for global_id, data in cluster.decoration_atoms.items():
            symbol = data['symbol']
            position = data['position']
            f.write(f'{symbol} {global_id} {position[0]:.5f} {position[1]:.5f} {position[2]:.5f} {cluster.root_id}\n')
            # Add decorating atoms to the map as well for bonding purposes
            node_id_to_local_index[global_id] = local_index
            local_index += 1
            
        return node_id_to_local_index
            
    def _write_cluster_bonds(self, f: TextIO, cluster: Cluster, id_map: Dict[int, int]) -> None:
        """Writes bond data, translating global node IDs to local file indices."""
        if self._settings.clustering.criteria == "distance":
            return
            
        # Write bonds between networking atoms
        for id1, id2 in cluster.linkages:
            if id1 in id_map and id2 in id_map:
                f.write(f"{id_map[id1]} {id_map[id2]}\n")

        # Write bonds between networking atoms and their decorating neighbors
        if self._settings.clustering.criteria == 'bond':
            bridge_symbol = self._settings.clustering.connectivity[1]
            cluster_nodes_map = {node.node_id: node for node in cluster.nodes}
            
            for node_id in cluster.indices:
                original_node = cluster_nodes_map.get(node_id)
                if not original_node: continue
                
                for neighbor in original_node.neighbors:
                    if neighbor.symbol == bridge_symbol:
                        # Check if both atoms are in the current file's map
                        if node_id in id_map and neighbor.node_id in id_map:
                             f.write(f"{id_map[node_id]} {id_map[neighbor.node_id]}\n")
