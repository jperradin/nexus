"""
This file contains all the methods / functions that are specific to SiOz-SiOz clusters.
"""
# external imports
import  numpy as np
from    tqdm import tqdm

# internal imports
from ..core.atom            import Atom
from ..core.box             import Box
from ..core.cluster         import Cluster
from ..utils.generate_color_gradient import generate_color_gradient


# List of supported elements for the extension SiOz
LIST_OF_SUPPORTED_ELEMENTS = ["Si", "O"]
LIST_OF_EXTRA_CONNECTIVITIES = ["stishovite"]
EXTRA_CLUSTERING_METHODS = True

class Silicon(Atom):
    def __init__(self, element, id, position, frame, cutoffs, extension) -> None:
        super().__init__(element, id, position, frame, cutoffs, extension)
        self.number_of_edges = 0 
    
    def get_number_of_edges(self) -> int:
        """
        Return the number of edges sharing
        """
        return self.number_of_edges
    
    def calculate_coordination(self) -> int:
        """
        Calculate the coordination number of the atom (ie the number of first neighbours) for the extension SiOz
        """
        self.coordination = len([neighbour for neighbour in self.neighbours if neighbour.get_element() == "O"])
            

class Oxygen(Atom):
    def __init__(self, element, id, position, frame, cutoffs, extension) -> None:
        super().__init__(element, id, position, frame, cutoffs, extension)
    
    def calculate_coordination(self) -> int:
        """
        Calculate the coordination number of the atom (ie the number of first neighbours) for the extension SiOz
        """
        self.coordination = len([neighbour for neighbour in self.neighbours if neighbour.get_element() == "Si"])
                
def transform_into_subclass(atom:Atom) -> object:
    """
    Return a Silicon object or Oxygen object from the subclass Silicon or Oxygen whether the atom.element is 'Si' or 'O'.  
    """
    if atom.get_element() == 'O':
        return Oxygen(atom.element, atom.id, atom.position, atom.frame, atom.cutoffs, atom.extension)
    elif atom.get_element() == 'Si':
        return Silicon(atom.element, atom.id, atom.position, atom.frame, atom.cutoffs, atom.extension)
    else:
        raise ValueError(f"\tERROR: Atom {atom.element} - {atom.id} can be transformed into Silicon or Oxygen object.")
        

def get_connectivity(cluster_settings) -> list:
        
    polyhedra = [v for k,v in cluster_settings.items() if k == "polyhedra"][0]
    list = []
    connectivity = cluster_settings['connectivity']
    if connectivity == ['Si', 'Si']:
        for poly in polyhedra:
            list.append(f'Si{poly[0]}-Si{poly[1]}')
    elif connectivity == ['Si', 'O', 'Si']:
        for poly in polyhedra:
            list.append(f'SiO{poly[0]}-SiO{poly[1]}')
    elif connectivity == ['O', 'Si', 'O']:
        for poly in polyhedra:
            list.append(f"OSi{poly[0]}-OSi{poly[1]}")
    elif connectivity == ['O', 'O']:
        for poly in polyhedra:
            list.append(f"O{poly[0]}-O{poly[1]}")
    return list

def get_extra_connectivity(cluster_settings) -> list:
    list = []
    extra_clusters = cluster_settings['extra_clusters']
    if extra_clusters == 'stishovite':
        if cluster_settings['criteria'] == 'bond':
            list.append('SiO6-SiO6-stishovite')
        elif cluster_settings['criteria'] == 'distance':
            list.append('Si6-Si6-stishovite')
    return list


def get_default_settings(criteria="bond") -> dict:
    """
    Method that load the default parameters for extension SiOz.
    """
    # internal imports
    from ..settings.parameter import Parameter, ClusterParameter
    
    # Structure of the system
    list_of_elements = [
                {"element": "Si", "alias": 2, "number": 0},
                {"element": "O" , "alias": 1, "number": 0}
            ]
            
    # Cluster settings to be set
    if criteria == "bond":
        dict_cluster_settings = {
            "connectivity"  : ["Si", "O", "Si"],
            "criteria"      : "bond",
            "polyhedra"     :
                [
                    [4, 4],
                    [4, 5],
                    [5, 5],
                    [5, 6],
                    [6, 6],
                    [6, 7]
                ],
            "extra_clusters" : "stishovite" 
                 
        }
    elif criteria == "distance":
        dict_cluster_settings = {
            "connectivity"  : ["Si", "Si"],
            "criteria"      : "distance",   # WARNING: if this criteria is set, 
                                            # the pair cutoff Si-Si will be used as the distance cutoff between the nodes.
            "polyhedra"     :
                [
                    [4, 4],
                    [4, 5],
                    [5, 5],
                    [5, 6],
                    [6, 6],
                    [6, 7]
                ],
            "extra_clusters" : "stishovite" 
        }
    else:
        raise ValueError(f"{criteria} not supported. Criteria must be \"bond\" or \"distance\".")
    
    
    # Pair cutoffs for the clusters
    list_of_cutoffs = [
        { "element1": "O" , "element2": "O" , "value": 3.05},
        { "element1": "Si", "element2": "O" , "value": 2.30},
        { "element1": "Si", "element2": "Si", "value": 3.50}
    ]
    
    # Settings
    dict_settings = {
        "extension": Parameter("extension", "SiOz"),
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
        - SiO4      : list of SiO4 tetrahedra
        - SiO5      : list of SiO5 pentahedra
        - SiO6      : list of SiO6 octahedra
        - SiO7      : list of SiO7 heptahedra 
        - OSi1      : list of OSi1
        - OSi2      : list of OSi2
        - OSi3      : list of OSi3
        - OSi4      : list of OSi4
        - ES_SiO6   : proportion of edge-sharing in SiO6 units
    """
    
    # Initialize the lists 
    SiO4 = []
    SiO5 = []
    SiO6 = []
    SiO7 = []
    OSi1 = []
    OSi2 = []
    OSi3 = []
    OSi4 = []
    ES_SiO6 = []
    
    silicons = [atom for atom in atoms if atom.get_element() == "Si"]
    oxygens  = [atom for atom in atoms if atom.get_element() == "O" ]
    
    number_of_Si = len(silicons)
    number_of_O  = len(oxygens)
    
    if criteria == 'bond':
        dict_concentrations = {
            "SiO4-SiO4" : [],
            "SiO4-SiO5" : [],
            "SiO5-SiO5" : [],
            "SiO5-SiO6" : [],
            "SiO6-SiO6" : [],
            "SiO6-SiO7" : [],
            "SiO7-SiO7" : [],
            "SiO6-SiO6-stishovite" : [],
            "OSi1-OSi2" : [],
            "OSi2-OSi2" : [],
            "OSi2-OSi3" : [],
            "OSi3-OSi3" : [],
            "OSi3-OSi4" : [],
            "OSi4-OSi4" : [] 
        }
    elif criteria == 'distance':
        dict_concentrations = {
            "Si4-Si4" : [],
            "Si4-Si5" : [],
            "Si5-Si5" : [],
            "Si5-Si6" : [],
            "Si6-Si6" : [],
            "Si6-Si7" : [],
            "Si7-Si7" : [],
            "Si6-Si6-stishovite" : [],
            "O1-O2" : [],
            "O2-O2" : [],
            "O2-O3" : [],
            "O3-O3" : [],
            "O3-O4" : [],
            "O4-O4" : [],
        }
    # Calculate the proportion of each SiOz units
    coordination_SiOz = []
    for atom in silicons:
        counter = len([neighbour for neighbour in atom.get_neighbours() if neighbour.get_element() == "O"])
        coordination_SiOz.append(counter)
        if counter == 4:
            SiO4.append(atom)
        if counter == 5:
            SiO5.append(atom)
        if counter == 6:
            SiO6.append(atom)
        if counter == 7:
            SiO7.append(atom)
    
    _debug_histogram_proportion_SiOz = np.histogram(coordination_SiOz, bins=[4,5,6,7,8], density=True) 
    
    # Calculate the proportion of each OSiz units
    coordination_OSiz = []
    for atom in oxygens:
        counter = len([neighbour for neighbour in atom.get_neighbours() if neighbour.get_element() == "Si"])
        coordination_OSiz.append(counter)
        if counter == 1:
            OSi1.append(atom)
        if counter == 2:
            OSi2.append(atom)
        if counter == 3:
            OSi3.append(atom)
        if counter == 4:
            OSi4.append(atom)
            
    _debug_histogram_proportion_OSik = np.histogram(coordination_OSiz, bins=[1,2,3,4,5], density=True) 
    
    # Calculate the number of edge-sharing (2 oxygens shared by 2 silicons)
    if quiet == False:
        progress_bar = tqdm(silicons, desc="Calculating the concentrations SiOz-SiOz sites", colour='BLUE', leave=False)
        color_gradient = generate_color_gradient(len(silicons))
        counter = 0
    else:
        progress_bar = silicons
        
    for silicon in progress_bar:
        if quiet == False:
            progress_bar.set_description(f"Calculating the concentrations SiOz-SiOz sites ...")
            progress_bar.colour = "#%02x%02x%02x" % color_gradient[counter]
            counter += 1
            
        unique_bond = []
        for oxygen in [atom for atom in silicon.get_neighbours() if atom.get_element() == 'O']:
            for second_silicon in [atom for atom in oxygen.get_neighbours() if atom.get_element() == 'Si']:
                if second_silicon.id != silicon.id:
                    unique_bond.append(second_silicon.id)
                if silicon.coordination == 4 and second_silicon.coordination == 4:
                    if criteria == 'bond':
                        dict_concentrations["SiO4-SiO4"].append(silicon.id)
                        dict_concentrations["SiO4-SiO4"].append(second_silicon.id)
                    elif criteria == 'distance':
                        dict_concentrations["Si4-Si4"].append(silicon.id)
                        dict_concentrations["Si4-Si4"].append(second_silicon.id)
                if silicon.coordination == 4 and second_silicon.coordination == 5:
                    if criteria == 'bond':
                        dict_concentrations["SiO4-SiO5"].append(silicon.id)
                        dict_concentrations["SiO4-SiO5"].append(second_silicon.id)
                    elif criteria == 'distance':
                        dict_concentrations["Si4-Si5"].append(silicon.id)
                        dict_concentrations["Si4-Si5"].append(second_silicon.id)
                if silicon.coordination == 5 and second_silicon.coordination == 5:
                    if criteria == 'bond':
                        dict_concentrations["SiO5-SiO5"].append(silicon.id)
                        dict_concentrations["SiO5-SiO5"].append(second_silicon.id)
                    elif criteria == 'distance':
                        dict_concentrations["Si5-Si5"].append(silicon.id)
                        dict_concentrations["Si5-Si5"].append(second_silicon.id)
                if silicon.coordination == 5 and second_silicon.coordination == 6:
                    if criteria == 'bond':
                        dict_concentrations["SiO5-SiO6"].append(silicon.id)
                        dict_concentrations["SiO5-SiO6"].append(second_silicon.id)
                    elif criteria == 'distance':
                        dict_concentrations["Si5-Si6"].append(silicon.id)
                        dict_concentrations["Si5-Si6"].append(second_silicon.id)
                if silicon.coordination == 6 and second_silicon.coordination == 6:
                    if criteria == 'bond':
                        dict_concentrations["SiO6-SiO6"].append(silicon.id)
                        dict_concentrations["SiO6-SiO6"].append(second_silicon.id)
                    elif criteria == 'distance':
                        dict_concentrations["Si6-Si6"].append(silicon.id)
                        dict_concentrations["Si6-Si6"].append(second_silicon.id)
                if silicon.coordination == 6 and second_silicon.coordination == 7:
                    if criteria == 'bond':
                        dict_concentrations["SiO6-SiO7"].append(silicon.id)
                        dict_concentrations["SiO6-SiO7"].append(second_silicon.id)
                    elif criteria == 'distance':
                        dict_concentrations["Si6-Si7"].append(silicon.id)
                        dict_concentrations["Si6-Si7"].append(second_silicon.id)
                if silicon.coordination == 7 and second_silicon.coordination == 7:
                    if criteria == 'bond':
                        dict_concentrations["SiO7-SiO7"].append(silicon.id)
                        dict_concentrations["SiO7-SiO7"].append(second_silicon.id)
                    elif criteria == 'distance':
                        dict_concentrations["Si7-Si7"].append(silicon.id)
                        dict_concentrations["Si7-Si7"].append(second_silicon.id)
                
        unique_bond = np.array(unique_bond)
        
        uniques, counts = np.unique(unique_bond, return_counts=True)
        
        for connectivity in counts:
            if connectivity == 2: # 2 oxygens are shared by 'silicon' and 'second_silicon'
                silicon.number_of_edges += 1

        if silicon.number_of_edges >= 2 and silicon.coordination == 6:
            for oxygen in [atom for atom in silicon.get_neighbours() if atom.get_element() == 'O']:
                for second_silicon in [atom for atom in oxygen.get_neighbours() if atom.get_element() == 'Si']:
                    if second_silicon.id != silicon.id:
                        if second_silicon.number_of_edges >= 2 and second_silicon.coordination == 6:
                            ES_SiO6.append(silicon)
                            ES_SiO6.append(second_silicon)
                            if criteria == 'bond':
                                dict_concentrations["SiO6-SiO6-stishovite"].append(silicon.id)
                                dict_concentrations["SiO6-SiO6-stishovite"].append(second_silicon.id)
                            elif criteria == 'distance':
                                dict_concentrations["Si6-Si6-stishovite"].append(silicon.id)
                                dict_concentrations["Si6-Si6-stishovite"].append(second_silicon.id)

    if quiet == False:
        progress_bar = tqdm(oxygens, desc="Calculating the concentrations SiOz-SiOz sites", colour='BLUE', leave=False)
        color_gradient = generate_color_gradient(len(oxygens))
        counter = 0
    else:
        progress_bar = oxygens
        
    for oxygen in progress_bar:
        if quiet == False:
            progress_bar.set_description(f"Calculating the concentrations OSiz-OSiz sites ...")
            progress_bar.colour = "#%02x%02x%02x" % color_gradient[counter]
            counter += 1
        for silicon in [atom for atom in oxygen.get_neighbours() if atom.get_element() == 'Si']:
            for second_oxygen in [atom for atom in silicon.get_neighbours() if atom.get_element() == 'O']:
                if second_oxygen.id != oxygen.id:
                    if oxygen.coordination == 1 and second_oxygen.coordination == 2:
                        if criteria == 'bond':
                            dict_concentrations["OSi1-OSi2"].append(oxygen.id)
                            dict_concentrations["OSi1-OSi2"].append(second_oxygen.id)
                        elif criteria == 'distance':
                            dict_concentrations["O1-O2"].append(oxygen.id)
                            dict_concentrations["O1-O2"].append(second_oxygen.id)
                    if oxygen.coordination == 2 and second_oxygen.coordination == 2:
                        if criteria == 'bond':
                            dict_concentrations["OSi2-OSi2"].append(oxygen.id)
                            dict_concentrations["OSi2-OSi2"].append(second_oxygen.id)
                        elif criteria == 'distance':
                            dict_concentrations["O2-O2"].append(oxygen.id)
                            dict_concentrations["O2-O2"].append(second_oxygen.id)
                    if oxygen.coordination == 2 and second_oxygen.coordination == 3:
                        if criteria == 'bond':
                            dict_concentrations["OSi2-OSi3"].append(oxygen.id)
                            dict_concentrations["OSi2-OSi3"].append(second_oxygen.id)
                        elif criteria == 'distance':
                            dict_concentrations["O2-O3"].append(oxygen.id)
                            dict_concentrations["O2-O3"].append(second_oxygen.id)
                    if oxygen.coordination == 3 and second_oxygen.coordination == 3:
                        if criteria == 'bond':
                            dict_concentrations["OSi3-OSi3"].append(oxygen.id)
                            dict_concentrations["OSi3-OSi3"].append(second_oxygen.id)
                        elif criteria == 'distance':
                            dict_concentrations["O3-O3"].append(oxygen.id)
                            dict_concentrations["O3-O3"].append(second_oxygen.id)
                    if oxygen.coordination == 3 and second_oxygen.coordination == 4:
                        if criteria == 'bond':
                            dict_concentrations["OSi3-OSi4"].append(oxygen.id)
                            dict_concentrations["OSi3-OSi4"].append(second_oxygen.id)
                        elif criteria == 'distance':
                            dict_concentrations["O3-O4"].append(oxygen.id)
                            dict_concentrations["O3-O4"].append(second_oxygen.id)
                    if oxygen.coordination == 4 and second_oxygen.coordination == 4:
                        if criteria == 'bond':
                            dict_concentrations["OSi4-OSi4"].append(oxygen.id)
                            dict_concentrations["OSi4-OSi4"].append(second_oxygen.id)
                        elif criteria == 'distance':
                            dict_concentrations["O4-O4"].append(oxygen.id)
                            dict_concentrations["O4-O4"].append(second_oxygen.id)
    
    if number_of_Si == 0: number_of_Si = 1 # Avoid division by zero
    if number_of_O == 0: number_of_O = 1 # Avoid division by zero
    
    # Calculate the concentrations of each connectivity
    for key, value in dict_concentrations.items():
        if key[0] == 'S':
            dict_concentrations[key] = len(np.unique(value)) / number_of_Si 
        elif key[0] == 'O':
            dict_concentrations[key] = len(np.unique(value)) / number_of_O
    return dict_concentrations

def find_extra_clusters(atoms:list, box:Box, counter_c:int, settings:object) -> None:
    r"""
    Find stishovite clusters, ie SiO6-SiO6 connected through two edge-sharings.
    
    Parameters:
    -----------
        - atoms (list) : List of Atom objects (ie SiO6 with two edge-sharings units)
        - criteria : str, "bond" or "distance"
        - counter_c : int, counter for the cluster id
        - settings : Settings object
    
    Returns:
    --------
        - local_clusters (list) : List of stishovite clusters (Cluster objects)
        - counter_c (int) : updated counter for the cluster id
    """
    
    # Define function for union-find algorithm
    def find(atom) -> Atom:
        r"""
        Find the root of the cluster to which the given atom belongs.
        
        Parameters:
        -----------
            - atom (Atom): Atom object for which to find the root.
        
        Returns:
        --------
            - root (Atom): Root Atom object of the cluster
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
    
    
    cluster_settings = settings.cluster_settings.get_value()
    criteria = cluster_settings['criteria']
    
    if criteria != 'bond' and criteria != 'distance':
        raise ValueError(f"Criteria {criteria} not supported. Criteria must be 'bond' or 'distance'.")
    
    if criteria == 'bond':
        node_1 = 'Si'
        bridge = 'O'
        node_2 = 'Si'
        chain = [node_1, bridge, node_2]
        connectivity = 'SiO6-SiO6-stishovite'
            
    if criteria == 'distance':
        node_1 = 'Si'
        node_2 = 'Si'
        chain = [node_1, node_2]
        connectivity = 'Si6-Si6-stishovite'
        
    networking_atoms = [atom for atom in atoms if (
        atom.get_element() == 'Si'
        and atom.get_number_of_edges() >= 2
        and atom.get_coordination() == 6
        )]
    
    number_of_nodes = 0
    
    color_gradient = generate_color_gradient(len(networking_atoms))
    if settings.quiet.get_value() == False:
        progress_bar = tqdm(networking_atoms, desc=f"Finding {connectivity} clusters", colour='BLUE', leave=False)
    else:
        progress_bar = networking_atoms
    colour = 0
    for atom in progress_bar:
        # Update progress bar
        if settings.quiet.get_value() == False:
            progress_bar.set_description(f"Finding {connectivity} clusters ...")
            progress_bar.colour = "#%02x%02x%02x" % color_gradient[colour]
            colour += 1
        
        if criteria == 'distance':
            for neighbour in atom.get_neighbours():
                if (atom.element == node_1 and neighbour.element == node_2
                    and atom.coordination == 6 and neighbour.coordination == 6
                    and atom.number_of_edges >= 2 and neighbour.number_of_edges >= 2):
                    union(neighbour, atom)
        
        if criteria == 'bond':
            for neighbour in atom.neighbours:
                if neighbour.element == bridge:
                    for second_neighbour in neighbour.neighbours:
                        if (atom.element == node_1 and second_neighbour.element == node_2
                            and atom.coordination == 6 and second_neighbour.coordination == 6
                            and atom.number_of_edges >= 2 and second_neighbour.number_of_edges >= 2): 
                            union(second_neighbour, atom)
        
    clusters_found = {}
    local_clusters = []
    
    for atom in networking_atoms:
        root = find(atom)
        clusters_found.setdefault(root.id, []).append(atom)
    
    color_gradient = generate_color_gradient(len(clusters_found))
    if settings.quiet.get_value() == False:
        progress_bar = tqdm(range(len(clusters_found)), desc=f"Calculating {connectivity} clusters properties ...", colour='GREEN', leave=False)
    else:
        progress_bar = range(len(clusters_found))
        
    colour = 0
    for i in progress_bar:
        cluster = list(clusters_found.values())[i]
        
        # Update progress bar
        if settings.quiet.get_value() == False:
            progress_bar.set_description(f"Calculating {connectivity} clusters properties ...")
            progress_bar.colour = "#%02x%02x%02x" % color_gradient[colour]
            colour += 1
        
        for atom in cluster:
            root = find(atom)
            break
        
        current_cluster = Cluster(box=box, connectivity=connectivity, root_id=counter_c, frame=root.frame, size=len(cluster))
        
        for atom in cluster:
            atom.set_cluster(counter_c, connectivity)
            current_cluster.add_atom(atom)
            if len(cluster) > 1:
                number_of_nodes += 1
        
        current_cluster.calculate_unwrapped_positions(criteria, chain, settings.quiet.get_value())
        current_cluster.calculate_center_of_mass()
        current_cluster.calculate_gyration_radius()
        current_cluster.calculate_percolation_probability()
        
        # if settings.print_clusters_positions.get_value():
        #     current_cluster.write_coordinates(settings._output_directory)

        local_clusters.append(current_cluster)
        counter_c += 1
        
        for atom in cluster:
            atom.reset_cluster()
    
    if number_of_nodes == 0: number_of_nodes = 1 # Avoid division by zero
    
    for cluster in local_clusters:
        cluster.number_of_nodes = number_of_nodes
        cluster.calculate_order_parameter()
    
    return local_clusters, counter_c         
    
