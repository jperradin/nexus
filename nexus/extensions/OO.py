"""
This file contains all the methods / functions that are specific to O-O clusters.
"""

# external imports
import numpy as np
from tqdm import tqdm

# internal imports
from ..core.atom import Atom
from ..core.box import Box
from ..core.cluster import Cluster
from ..utils.generate_color_gradient import generate_color_gradient


# List of supported elements for the extension OO
LIST_OF_SUPPORTED_ELEMENTS = ["O"]
LIST_OF_EXTRA_CONNECTIVITIES = ["LD", "HD", "VHD", "HV"]
EXTRA_CLUSTERING_METHODS = True


class Oxygen(Atom):
    def __init__(self, element, id, position, frame, cutoffs, extension) -> None:
        super().__init__(element, id, position, frame, cutoffs, extension)

    def calculate_coordination(self) -> int:
        """
        Calculate the coordination number of the atom (ie the number of first neighbours) for the extension OO
        """
        self.coordination = len(
            [
                neighbour
                for neighbour in self.neighbours
                if neighbour.get_element() == "O"
            ]
        )


def transform_into_subclass(atom: Atom) -> object:
    """
    Return a Oxygen object from the subclass Oxygen whether the atom.element is or 'O'.
    """
    if atom.get_element() == "O":
        return Oxygen(
            atom.element,
            atom.id,
            atom.position,
            atom.frame,
            atom.cutoffs,
            atom.extension,
        )
    else:
        raise ValueError(
            f"\tERROR: Atom {atom.element} - {atom.id} can be transformed into Oxygen object."
        )


def get_connectivity(cluster_settings) -> list:
    polyhedra = [v for k, v in cluster_settings.items() if k == "polyhedra"][0]
    list = []
    for poly in polyhedra:
        list.append(f"O{poly[0]}-O{poly[1]}")
    return list


def get_extra_connectivity(cluster_settings) -> list:
    return cluster_settings["extra_clusters"]


def get_default_settings(criteria="distance") -> dict:
    """
    Method that load the default parameters for extension OO.
    """
    # internal imports
    from ..settings.parameter import Parameter, ClusterParameter

    # Structure of the system
    list_of_elements = [
        {"element": "O", "alias": 1, "number": 0},
    ]

    # Cluster settings to be set
    if criteria == "distance":
        dict_cluster_settings = {
            "connectivity": ["O", "O"],
            "criteria": "distance",  # WARNING: if this criteria is set,
            # the pair cutoff O-O will be used as the distance cutoff between the sites.
            "polyhedra": [
                [3, 3],
                [3, 4],
                [4, 4],
                [4, 5],
                [5, 5],
                [5, 6],
                [6, 6],
                [6, 7],
                [7, 7],
                [7, 8],
                [8, 8],
                [8, 9],
                [9, 9],
            ],
            "extra_clusters": ["LD", "HD", "VHD", "HV"],
        }
    else:
        raise ValueError(f'{criteria} not supported. Criteria must be "distance".')

    # Pair cutoffs for the clusters
    list_of_cutoffs = [{"element1": "O", "element2": "O", "value": 3.50}]

    # Settings
    dict_settings = {
        "extension": Parameter("extension", "OO"),
        "structure": Parameter("structure", list_of_elements),
        "cluster_settings": ClusterParameter("cluster_settings", dict_cluster_settings),
        "cutoffs": Parameter("cutoffs", list_of_cutoffs),
    }

    return dict_settings


def calculate_concentrations(atoms: list, criteria: str, quiet: bool) -> dict:
    """
    Calculate the following properties.

    Returns:
    --------
        - dict_concentrations  : dictionnary of concentrations for each cluster connectivity.
    """

    # Initialize the lists
    O2 = []
    O3 = []
    O4 = []
    O5 = []
    O6 = []
    O7 = []
    O8 = []
    O9 = []

    LD = []   # z <= 4
    HD = []   # 4 < z < 8
    VHD = []  # z >= 8
    HV = []   # z > 4

    if criteria == "distance":
        dict_concentrations = {
            "O2-O2": [],
            "O2-O3": [],
            "O3-O3": [],
            "O3-O4": [],
            "O4-O4": [],
            "O4-O5": [],
            "O5-O5": [],
            "O5-O6": [],
            "O6-O6": [],
            "O6-O7": [],
            "O7-O7": [],
            "O7-O8": [],
            "O8-O8": [],
            "O8-O9": [],
            "O9-O9": [],
            "LD": [],
            "HD": [],
            "VHD": [],
            "HV": [],
        }
    # Calculate the proportion of each OO units
    coordination_OO = []
    for atom in atoms:
        counter = len(
            [
                neighbour
                for neighbour in atom.get_neighbours()
                if neighbour.get_element() == "O"
            ]
        )
        coordination_OO.append(counter)
        if counter == 2:
            O2.append(atom)
        if counter == 3:
            O3.append(atom)
        if counter == 4:
            O4.append(atom)
        if counter == 5:
            O5.append(atom)
        if counter == 6:
            O6.append(atom)
        if counter == 7:
            O7.append(atom)
        if counter == 8:
            O8.append(atom)
        if counter == 9:
            O9.append(atom)
        if atom.coordination <= 4:
            LD.append(atom)
        if 4 < atom.coordination < 8:
            HD.append(atom)
        if atom.coordination >= 8:
            VHD.append(atom)
        if atom.coordination > 4:
            HV.append(atom)

    _debug_histogram_proportion_OO = np.histogram(
        coordination_OO, bins=[2, 3, 4, 5, 6, 7, 8, 9, 10], density=True
    )

    if quiet == False:
        progress_bar = tqdm(
            atoms,
            desc="Calculating the concentrations O-O sites",
            colour="BLUE",
            leave=False,
        )
        color_gradient = generate_color_gradient(len(atoms))
        counter = 0
    else:
        progress_bar = atoms

    for atom in progress_bar:
        if quiet == False:
            progress_bar.set_description(
                f"Calculating the concentrations O-O sites ..."
            )
            progress_bar.colour = "#%02x%02x%02x" % color_gradient[counter]
            counter += 1

        for neighbor in atom.get_neighbours():
            if atom.coordination == 2 and neighbor.coordination == 2:
                dict_concentrations["O2-O2"].append(atom.id)
                dict_concentrations["O2-O2"].append(neighbor.id)
            if atom.coordination == 2 and neighbor.coordination == 3:
                dict_concentrations["O2-O3"].append(atom.id)
                dict_concentrations["O2-O3"].append(neighbor.id)
            if atom.coordination == 3 and neighbor.coordination == 3:
                dict_concentrations["O3-O3"].append(atom.id)
                dict_concentrations["O3-O3"].append(neighbor.id)
            if atom.coordination == 3 and neighbor.coordination == 4:
                dict_concentrations["O3-O4"].append(atom.id)
                dict_concentrations["O3-O4"].append(neighbor.id)
            if atom.coordination == 4 and neighbor.coordination == 4:
                dict_concentrations["O4-O4"].append(atom.id)
                dict_concentrations["O4-O4"].append(neighbor.id)
            if atom.coordination == 4 and neighbor.coordination == 5:
                dict_concentrations["O4-O5"].append(atom.id)
                dict_concentrations["O4-O5"].append(neighbor.id)
            if atom.coordination == 5 and neighbor.coordination == 5:
                dict_concentrations["O5-O5"].append(atom.id)
                dict_concentrations["O5-O5"].append(neighbor.id)
            if atom.coordination == 5 and neighbor.coordination == 6:
                dict_concentrations["O5-O6"].append(atom.id)
                dict_concentrations["O5-O6"].append(neighbor.id)
            if atom.coordination == 6 and neighbor.coordination == 6:
                dict_concentrations["O6-O6"].append(atom.id)
                dict_concentrations["O6-O6"].append(neighbor.id)
            if atom.coordination == 6 and neighbor.coordination == 7:
                dict_concentrations["O6-O7"].append(atom.id)
                dict_concentrations["O6-O7"].append(neighbor.id)
            if atom.coordination == 7 and neighbor.coordination == 7:
                dict_concentrations["O7-O7"].append(atom.id)
                dict_concentrations["O7-O7"].append(neighbor.id)
            if atom.coordination == 7 and neighbor.coordination == 8:
                dict_concentrations["O7-O8"].append(atom.id)
                dict_concentrations["O7-O8"].append(neighbor.id)
            if atom.coordination == 8 and neighbor.coordination == 8:
                dict_concentrations["O8-O8"].append(atom.id)
                dict_concentrations["O8-O8"].append(neighbor.id)
            if atom.coordination == 8 and neighbor.coordination == 9:
                dict_concentrations["O8-O9"].append(atom.id)
                dict_concentrations["O8-O9"].append(neighbor.id)
            if atom.coordination == 9 and neighbor.coordination == 9:
                dict_concentrations["O9-O9"].append(atom.id)
                dict_concentrations["O9-O9"].append(neighbor.id)
            if atom.coordination <= 4 and neighbor.coordination <= 4:
                dict_concentrations["LD"].append(atom.id)
                dict_concentrations["LD"].append(neighbor.id)
            if 4 < atom.coordination < 8 and 4 < neighbor.coordination < 8:
                dict_concentrations["HD"].append(atom.id)
                dict_concentrations["HD"].append(neighbor.id)
            if atom.coordination >= 8 and neighbor.coordination >= 8:
                dict_concentrations["VHD"].append(atom.id)
                dict_concentrations["VHD"].append(neighbor.id)
            if atom.coordination > 4 and neighbor.coordination > 4:
                dict_concentrations["HV"].append(atom.id)
                dict_concentrations["HV"].append(neighbor.id)

    # Calculate the concentrations of each connectivity
    for key, value in dict_concentrations.items():
        dict_concentrations[key] = len(np.unique(value)) / len(atoms)
    return dict_concentrations


def find_extra_clusters(
    atoms: list, box: Box, counter_c: int, settings: object
) -> None:
    r"""
        Find LD, HD, VHD and HV clusters in the system.
        LD is defined as a cluster of O atoms with a coordination number <= 4.
        HD is defined as a cluster of O atoms with a coordination number > 4 and < 8.
        VHD is defined as a cluster of O atoms with a coordination number >= 8.
        HV is defined as a cluster of O atoms with a coordination number > 4.

    Parameters:
    -----------
        - atoms      : list of Atom objects
        - box        : Box object
        - counter_c  : int, counter for the cluster id
        - settings   : object, settings object

    Returns:
    --------
        - local_clusters : list of Cluster objects
        - counter_c      : int, updated counter for the cluster id
    """

    # Define function for union-find algorithm
    def find(atom) -> Atom:
        """
        Find the root of the cluster to which the given atom belongs.

        Parameters:
        -----------
            - atom : Atom object

        Returns:
        --------
            - root : Atom object
        """
        if atom.parent != atom:
            atom.parent = find(atom.parent)
        return atom.parent

    def union(atom_1, atom_2) -> None:
        r"""
        Union the two clusters to which the given atoms belong.

        Parameters:
        -----------
            - atom_1 (Atom): First Atom object.
            - atom_2 (Atom): Second Atom object

        Returns:
        --------
            - None.
        """
        root_1 = find(atom_1)
        root_2 = find(atom_2)

        if root_1 != root_2:
            root_2.parent = root_1

    # Get cluster settings and criteria
    cluster_settings = settings.cluster_settings.get_value()
    criteria = cluster_settings["criteria"]

    # Ensure criteria is 'distance'
    if criteria != "distance":
        raise ValueError(
            f"Criteria {criteria} not supported. Criteria must be 'distance'."
        )
    else:
        node_1 = "O"
        node_2 = "O"
        chain = [node_1, node_2]
        connectivities = ["LD", "HD", "VHD", "HV"]

    # Initialize networking atoms dictionary
    networking_atoms = {}
    networking_atoms["LD"] = [atom for atom in atoms if (atom.get_coordination() <= 4)]
    networking_atoms["HD"] = [
        atom for atom in atoms if (4 < atom.get_coordination() < 8)
    ]
    networking_atoms["VHD"] = [atom for atom in atoms if (atom.get_coordination() >= 8)]
    networking_atoms["HV"] = [atom for atom in atoms if (atom.get_coordination() > 4)]

    local_clusters = []

    # Iterate over each connectivity type
    for key in networking_atoms.keys():
        current_network = networking_atoms[key]
        number_of_nodes = 0

        # Generate color gradient for progress bar
        color_gradient = generate_color_gradient(len(current_network))
        if not settings.quiet.get_value():
            progress_bar = tqdm(
                current_network,
                desc=f"Finding {key} clusters",
                colour="BLUE",
                leave=False,
            )
        else:
            progress_bar = current_network
        colour = 0

        # Union-find algorithm to group atoms into clusters
        for atom in progress_bar:
            # Update progress_bar
            if not settings.quiet.get_value():
                progress_bar.set_description(f"Finding {key} clusters ...")
                progress_bar.colour = "#%02x%02x%02x" % color_gradient[colour]
                colour += 1

            for neighbour in atom.get_neighbours():
                if (
                    atom.element == node_1
                    and neighbour.element == node_2
                    and neighbour in current_network
                ):
                    union(neighbour, atom)

        clusters_found = {}
        current_local_clusters = []

        # Group atoms by their root
        for atom in current_network:
            root = find(atom)
            clusters_found.setdefault(root.id, []).append(atom)

        # Generate color gradient for progress bar
        color_gradient = generate_color_gradient(len(clusters_found))
        if not settings.quiet.get_value():
            progress_bar = tqdm(
                range(len(clusters_found)),
                desc=f"Calculating {key} clusters properties ...",
                colour="GREEN",
                leave=False,
            )
        else:
            progress_bar = range(len(clusters_found))

        colour = 0
        for i in progress_bar:
            cluster = list(clusters_found.values())[i]

            # Update progress bar
            if not settings.quiet.get_value():
                progress_bar.set_description(
                    f"Calculating {key} cluster properties ..."
                )
                progress_bar.colour = "#%02x%02x%02x" % color_gradient[colour]
                colour += 1

            for atom in cluster:
                root = find(atom)
                break

            # Create a new Cluster object
            current_cluster = Cluster(
                box=box,
                connectivity=key,
                root_id=counter_c,
                frame=root.frame,
                size=len(cluster),
            )

            # Add atoms to the cluster
            for atom in cluster:
                atom.set_cluster(counter_c, key)
                current_cluster.add_atom(atom)
                if len(cluster) > 1:
                    number_of_nodes += 1

            # Calculate cluster properties
            current_cluster.calculate_unwrapped_positions(
                criteria, chain, settings.quiet.get_value()
            )
            current_cluster.calculate_center_of_mass()
            current_cluster.calculate_gyration_radius()
            current_cluster.calculate_percolation_probability()

            if settings.print_clusters_positions.get_value():
                current_cluster.write_coordinates(settings._output_directory)

            current_local_clusters.append(current_cluster)
            counter_c += 1

            for atom in cluster:
                atom.reset_cluster()

        if number_of_nodes == 0:
            number_of_nodes = 1

        # Finalize cluster properties
        for cluster in current_local_clusters:
            cluster.number_of_nodes = number_of_nodes
            cluster.calculate_order_parameter()
            local_clusters.append(cluster)

    return local_clusters, counter_c
