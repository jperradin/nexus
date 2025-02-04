"""
This file contains all the methods / functions that are specific to O-O clusters.
"""
# external imports
import  numpy as np
from    tqdm import tqdm

# internal imports
from ..core.atom            import Atom
from ..core.box             import Box
from ..core.cluster         import Cluster
from ..utils.generate_color_gradient import generate_color_gradient


# List of supported elements for the extension OO
LIST_OF_SUPPORTED_ELEMENTS = ["O"]
EXTRA_CLUSTERING_METHODS = False

class Oxygen(Atom):
    def __init__(self, element, id, position, frame, cutoffs, extension) -> None:
        super().__init__(element, id, position, frame, cutoffs, extension)
    
    def calculate_coordination(self) -> int:
        """
        Calculate the coordination number of the atom (ie the number of first neighbours) for the extension OO
        """
        self.coordination = len([neighbour for neighbour in self.neighbours if neighbour.get_element() == "O"])
            
                
def transform_into_subclass(atom:Atom) -> object:
    """
    Return a Oxygen object from the subclass Oxygen whether the atom.element is or 'O'.  
    """
    if atom.get_element() == 'O':
        return Oxygen(atom.element, atom.id, atom.position, atom.frame, atom.cutoffs, atom.extension)
    else:
        raise ValueError(f"\tERROR: Atom {atom.element} - {atom.id} can be transformed into Oxygen object.")
        

def get_connectivity(cluster_settings) -> list:
        
    polyhedra = [v for k,v in cluster_settings.items() if k == "polyhedra"][0]
    list = []
    for poly in polyhedra:
        list.append(f'O{poly[0]}-O{poly[1]}')
    return list


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
            "connectivity"  : ["O", "O"],
            "criteria"      : "distance",   # WARNING: if this criteria is set, 
                                            # the pair cutoff O-O will be used as the distance cutoff between the nodes.
            "polyhedra"     :
                [
                    [1, 1],
                    [1, 2],
                    [2, 2],
                    [2, 3],
                    [3, 3],
                    [3, 4],
                    [4, 4],
                    [4, 5],
                    [5, 5],
                ]
        }
    else:
        raise ValueError(f"{criteria} not supported. Criteria must be \"distance\".")
    
    
    # Pair cutoffs for the clusters
    list_of_cutoffs = [
        { "element1": "O", "element2": "O", "value": 3.00}
    ]
    
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
        - O0      : list of O0
        - O1      : list of O1
        - O2      : list of O2 
        - O3      : list of O3 
        - O4      : list of O4 
        - O5      : list of O5 
    """
    
    # Initialize the lists 
    O0 = []
    O1 = []
    O2 = []
    O3 = []
    O4 = []
    O5 = []
    
    if criteria == 'distance':
        dict_concentrations = {
            "O1-O1" : [],
            "O1-O2" : [],
            "O2-O2" : [],
            "O2-O3" : [],
            "O3-O3" : [],
            "O3-O4" : [],
            "O4-O4" : [],
            "O4-O5" : [],
            "O5-O5" : [],
        }
    # Calculate the proportion of each OO units
    coordination_OO = []
    for atom in atoms:
        counter = len([neighbour for neighbour in atom.get_neighbours() if neighbour.get_element() == "O"])
        coordination_OO.append(counter)
        if counter == 0:
            O0.append(atom)
        if counter == 1:
            O1.append(atom)
        if counter == 2:
            O2.append(atom)
        if counter == 3:
            O3.append(atom)
        if counter == 4:
            O4.append(atom)
        if counter == 5:
            O5.append(atom)
            
    _debug_histogram_proportion_OO = np.histogram(coordination_OO, bins=[0,1,2,3,4,5,6], density=True) 
    
    if quiet == False:
        progress_bar = tqdm(atoms, desc="Calculating the concentrations O-O sites", colour='BLUE', leave=False)
        color_gradient = generate_color_gradient(len(atoms))
        counter = 0
    else:
        progress_bar = atoms
        
    for atom in progress_bar:
        if quiet == False:
            progress_bar.set_description(f"Calculating the concentrations O-O sites ...")
            progress_bar.colour = "#%02x%02x%02x" % color_gradient[counter]
            counter += 1
            
        for neighbor in atom.get_neighbours():
            if atom.coordination == 1 and neighbor.coordination == 1:
                dict_concentrations["O1-O1"].append(atom.id)
                dict_concentrations["O1-O1"].append(neighbor.id)
            if atom.coordination == 1 and neighbor.coordination == 2:
                dict_concentrations["O1-O2"].append(atom.id)
                dict_concentrations["O1-O2"].append(neighbor.id)
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
    
    
    # Calculate the concentrations of each connectivity
    for key, value in dict_concentrations.items():
        dict_concentrations[key] = len(np.unique(value)) / len(atoms) 
    return dict_concentrations