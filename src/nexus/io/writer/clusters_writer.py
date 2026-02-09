from typing import List, Dict, TextIO
from ...core.cluster import Cluster
from .base_writer import BaseWriter
from ...config.settings import Settings
import os
import numpy as np

class ClustersWriter(BaseWriter):
    """
    Writes unwrapped cluster coordinates and bond connectivity to files.

    For each cluster, generates an XYZ file with unwrapped positions of both the
    primary networking nodes and any decoration nodes (e.g., bridging nodes). When
    linkages are defined, a corresponding ``.bonds`` file is also written. Output
    layout depends on the configured ``print_mode``.

    Attributes:
        _clusters (List[Cluster]): Clusters to write, sorted by size descending.
    """

    def __init__(self, settings: Settings) -> None:
        """
        Initialize the clusters writer.

        Args:
            settings (Settings): Configuration settings.
        """
        super().__init__(settings)
        self._settings: Settings = settings

    def set_clusters(self, clusters: List[Cluster]) -> None:
        """
        Set the clusters to write, sorting them by size in descending order.

        Args:
            clusters (List[Cluster]): Clusters to write.
        """
        self._clusters: List[Cluster] = sorted(clusters, key=lambda cluster: cluster.size, reverse=True)

    def write(self) -> None:
        """Write cluster data according to the configured ``print_mode``."""
        if self._settings.clustering.print_mode == "none":
            return
        elif self._settings.clustering.print_mode == "all":
            self._write_all()
        elif self._settings.clustering.print_mode == "connectivity":
            self._write_connectivity()
        elif self._settings.clustering.print_mode == "individual":
            self._write_individual()

    def _write_all(self) -> None:
        """Write all clusters into a single XYZ and bonds file per frame."""
        path = os.path.join(self._settings.export_directory, "unwrapped_clusters")
        if not os.path.exists(path): os.makedirs(path)
        if not self._clusters: return
        
        frame_id = self._clusters[0].frame_id
        xyz_path = os.path.join(path, f'all_unwrapped_clusters-frame_{frame_id}.xyz')
        bonds_path = os.path.join(path, f'all_unwrapped_clusters-frame_{frame_id}.bonds')
        
        # Calculate total nodes for the header
        atoms_id_list = []
        for c in self._clusters:
            atoms_id_list.extend(c.indices)
            atoms_id_list.extend(c.decoration_atoms.keys())
        # Use np.unique to avoid double-counting any shared decoration atoms
        total_atoms = len(np.unique(atoms_id_list))
        uniques_ids = set()
        
        with open(xyz_path, 'w') as xyz_file, open(bonds_path, 'w') as bonds_file:
            xyz_file.write(f"{total_atoms}\n")
            self._write_header_comment(xyz_file, self._clusters[0])

            global_node_id_to_local_index = {}
            current_local_index = 1
            
            # Build the complete ID-to-index map for all nodes (networking and decorating)
            for cluster in self._clusters:
                for node_id in cluster.indices:
                    global_node_id_to_local_index[node_id] = current_local_index
                    current_local_index += 1
                # Decorating nodes are also part of the node count
                current_local_index += len(cluster.decoration_atoms)

            # Reset and write node positions and bond information
            current_local_index = 1
            for cluster in self._clusters:
                # Write nodes and update the map on the fly
                cluster_map, uniques_ids = self._write_cluster_atoms(xyz_file, cluster, current_local_index, uniques_ids)
                current_local_index += len(cluster_map)
                
                # Write bonds using the map for this cluster
                self._write_cluster_bonds(bonds_file, cluster, cluster_map)

    def _write_connectivity(self) -> None:
        """Write clusters grouped by connectivity, one file pair per connectivity."""
        clusters_per_connectivity: Dict[str, List[Cluster]] = {}
        for cluster in self._clusters:
            clusters_per_connectivity.setdefault(cluster.connectivity, []).append(cluster)

        for connectivity, cluster_list in clusters_per_connectivity.items():
            unique_ids = set()
            path = os.path.join(self._settings.export_directory, "unwrapped_clusters", connectivity)
            if not os.path.exists(path): os.makedirs(path)
            if not cluster_list: continue
                
            frame_id = cluster_list[0].frame_id
            xyz_path = os.path.join(path, f'{connectivity}_unwrapped_clusters-frame_{frame_id}.xyz')
            bonds_path = os.path.join(path, f'{connectivity}_unwrapped_clusters-frame_{frame_id}.bonds')
            
            atoms_id_list = []
            for c in cluster_list:
                atoms_id_list.extend(c.indices)
                atoms_id_list.extend(c.decoration_atoms.keys())
            # Use np.unique to avoid double-counting any shared decoration atoms
            total_atoms = len(np.unique(atoms_id_list))
            # total_atoms = sum(c.size + len(c.decoration_atoms) for c in cluster_list)

            with open(xyz_path, 'w') as xyz_file, open(bonds_path, 'w') as bonds_file:
                xyz_file.write(f"{total_atoms}\n")
                self._write_header_comment(xyz_file, cluster_list[0])

                current_local_index = 1
                for i, cluster in enumerate(cluster_list):
                    cluster.is_spanning = True if i == 0 else False
                    cluster_map, unique_ids = self._write_cluster_atoms(xyz_file, cluster, current_local_index, unique_ids)
                    current_local_index += len(cluster_map)
                    self._write_cluster_bonds(bonds_file, cluster, cluster_map)
        # debug
        hold = 1
            
    def _write_individual(self) -> None:
        """Write each cluster to its own XYZ and bonds file."""
        path = os.path.join(self._settings.export_directory, "unwrapped_clusters")
        if not os.path.exists(path): os.makedirs(path)
        if not self._clusters: return
        
        frame_id = self._clusters[0].frame_id
        for cluster in self._clusters:
            subpath = os.path.join(path, cluster.connectivity)
            if not os.path.exists(subpath): os.makedirs(subpath)
            
            xyz_path = os.path.join(subpath, f'cluster-frame_{frame_id}-id_{cluster.root_id}.xyz')
            bonds_path = os.path.join(subpath, f'cluster-frame_{frame_id}-id_{cluster.root_id}.bonds')
            
            total_atoms = cluster.size + len(cluster.decoration_atoms)
            with open(xyz_path, 'w') as xyz_file:
                xyz_file.write(f"{total_atoms}\n")
                self._write_header_comment(xyz_file, cluster)
                # For individual files, the local index starts at 1
                node_id_to_local_index, _ = self._write_cluster_atoms(xyz_file, cluster, 1)

            with open(bonds_path, 'w') as bonds_file:
                self._write_cluster_bonds(bonds_file, cluster, node_id_to_local_index)

    def _write_header_comment(self, f: TextIO, cluster: Cluster) -> None:
        """
        Write the XYZ comment/properties line with lattice and column descriptors.

        Args:
            f (TextIO): Open file handle to write to.
            cluster (Cluster): Cluster providing the lattice matrix.
        """
        lxx, lxy, lxz = cluster.lattice[0]
        lyx, lyy, lyz = cluster.lattice[1]
        lzx, lzy, lzz = cluster.lattice[2]
        lattice_line = f'Lattice="{lxx} {lxy} {lxz} {lyx} {lyy} {lyz} {lzx} {lzy} {lzz}" Properties=species:S:1:index:I:1:pos:R:3:cluster_id:I:1:coordination:I:1:percolating:I:1:spanning:I:1\n'
        f.write(lattice_line)

    def _write_cluster_atoms(self, f: TextIO, cluster: Cluster, start_index: int, unique_ids: List|None = None) -> Dict[int, int]:
        """
        Write all nodes (networking and decoration) for a cluster to an open file.

        Args:
            f (TextIO): Open file handle to write to.
            cluster (Cluster): The cluster whose nodes are written.
            start_index (int): Starting local index for node numbering.
            unique_ids (List | None): Set of already-written node IDs to avoid
                duplicates, or None to track internally.

        Returns:
            Dict[int, int]: Mapping of global node IDs to local file indices.
        """
        node_id_to_local_index = {}
        local_index = start_index
        if unique_ids is None:
            unique_ids = set()
        
        span = 1 if cluster.is_spanning else 0
        # Write primary networking nodes
        for symbol, global_id, position, node in zip(cluster.symbols, cluster.indices, cluster.unwrapped_positions, cluster.nodes):
            if node.coordination is None:
                coord = 0
            else:
                coord = node.coordination
            f.write(f'{symbol} {global_id} {position[0]:.5f} {position[1]:.5f} {position[2]:.5f} {cluster.root_id} {coord} {len(cluster.percolation_probability)} {span}\n')
            node_id_to_local_index[global_id] = local_index
            local_index += 1
            unique_ids.add(global_id)
            
        # Write decorating nodes
        for global_id, data in cluster.decoration_atoms.items():
            symbol = data['symbol']
            position = data['position']
            coordination = data['coordination']
            # Skip if this decorating atom was already written (shared between clusters ???)
            if global_id in unique_ids:
                continue
            unique_ids.add(global_id)
            if coordination is None:
                coord = 0
            else:
                coord = coordination
            f.write(f'{symbol} {global_id} {position[0]:.5f} {position[1]:.5f} {position[2]:.5f} {cluster.root_id} {coord} {len(cluster.percolation_probability)} {span}\n')
            # Add decorating nodes to the map as well for bonding purposes
            node_id_to_local_index[global_id] = local_index
            local_index += 1
        x = len(unique_ids)    
        return node_id_to_local_index, unique_ids
            
    def _write_cluster_bonds(self, f: TextIO, cluster: Cluster, id_map: Dict[int, int]) -> None:
        """
        Write bond connectivity data, translating global node IDs to local indices.

        Args:
            f (TextIO): Open file handle to write to.
            cluster (Cluster): The cluster whose bonds are written.
            id_map (Dict[int, int]): Mapping of global node IDs to local file indices.
        """
        if self._settings.clustering.criterion == "distance":
            return
        
        _nodes = {node.node_id: node for node in cluster.nodes}

        # Write bonds between networking nodes
        for id1, id2 in cluster.linkages:
            if id1 in id_map and id2 in id_map:
                s1 = _nodes[id1].symbol
                s2 = _nodes[id2].symbol
                f.write(f"{s1}({id_map[id1]})-{s2}({id_map[id2]})\n")
