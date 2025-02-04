"""
This file contains all the methods / functions that are specific to Si-Si clusters.
"""
# external imports
import  numpy as np
from    tqdm import tqdm

# internal imports
from ..core.atom            import Atom
from ..core.box             import Box
from ..core.cluster         import Cluster
from ..utils.generate_color_gradient import generate_color_gradient


# List of supported elements for the extension SiSi
LIST_OF_SUPPORTED_ELEMENTS = ["Si"]
EXTRA_CLUSTERING_METHODS = False

class Silicon(Atom):
    def __init__(self, element, id, position, frame, cutoffs, extension) -> None:
        super().__init__(element, id, position, frame, cutoffs, extension)
    
    def calculate_coordination(self) -> int:
        """
        Calculate the coordination number of the atom (ie the number of first neighbours) for the extension SiSi
        """
        self.coordination = len([neighbour for neighbour in self.neighbours if neighbour.get_element() == "Si"])
            
                
def transform_into_subclass(atom:Atom) -> object:
    """
    Return a Silicon object from the subclass Silicon.  
    """
    if atom.get_element() == 'Si':
        return Silicon(atom.element, atom.id, atom.position, atom.frame, atom.cutoffs, atom.extension)
    else:
        raise ValueError(f"\tERROR: Atom {atom.element} - {atom.id} can be transformed into Silicon object.")
        

def get_connectivity(cluster_settings) -> list:
        
    polyhedra = [v for k,v in cluster_settings.items() if k == "polyhedra"][0]
    list = []
    for poly in polyhedra:
        list.append(f'Si{poly[0]}-Si{poly[1]}')
    return list


def get_default_settings(criteria="distance") -> dict:
    """
    Method that load the default parameters for extension SiSi.
    """
    # internal imports
    from ..settings.parameter import Parameter, ClusterParameter
    
    # Structure of the system
    list_of_elements = [
                {"element": "Si", "alias": 1, "number": 0},
            ]
            
    # Cluster settings to be set
    if criteria == "distance":
        dict_cluster_settings = {
            "connectivity"  : ["Si", "Si"],
            "criteria"      : "distance",   # WARNING: if this criteria is set, 
                                            # the pair cutoff Si-Si will be used as the distance cutoff between the sites.
            "polyhedra"     :
                [
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
                ]
        }
    else:
        raise ValueError(f"{criteria} not supported. Criteria must be \"distance\".")
    
    
    # Pair cutoffs for the clusters
    list_of_cutoffs = [
        { "element1": "Si", "element2": "Si", "value": 3.00}
    ]
    
    # Settings
    dict_settings = {
        "extension": Parameter("extension", "SiSi"),
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
        - Si3      : list of Si3 
        - Si4      : list of Si4 
        - Si5      : list of Si5 
        - Si6      : list of Si6
        - Si7      : list of Si7  
        - Si8      : list of Si8  
        - Si9      : list of Si9  
    """
    
    # Initialize the lists 
    Si3 = []
    Si4 = []
    Si5 = []
    Si6 = []
    Si7 = []
    Si8 = []
    Si9 = []
    
    if criteria == 'distance':
        dict_concentrations = {
            "Si3-Si3" : [],
            "Si3-Si4" : [],
            "Si4-Si4" : [],
            "Si4-Si5" : [],
            "Si5-Si5" : [],
            "Si5-Si6" : [],
            "Si6-Si6" : [],
            "Si6-Si7" : [],
            "Si7-Si7" : [],
            "Si7-Si8" : [],
            "Si8-Si8" : [],
            "Si8-Si9" : [],
            "Si9-Si9" : [],
        }
    # Calculate the proportion of each SiSi units
    coordination_SiSi = []
    for atom in atoms:
        counter = len([neighbour for neighbour in atom.get_neighbours() if neighbour.get_element() == "Si"])
        coordination_SiSi.append(counter)
        if counter == 3:
            Si3.append(atom)
        if counter == 4:
            Si4.append(atom)
        if counter == 5:
            Si5.append(atom)
        if counter == 6:
            Si6.append(atom)
        if counter == 7:
            Si7.append(atom)
        if counter == 8:
            Si8.append(atom)    
        if counter == 9:
            Si9.append(atom)
            
    _debug_histogram_proportion_SiSi = np.histogram(coordination_SiSi, bins=[3,4,5,6,7,8,9,10], density=True) 
    
    # Calculate the number of edge-sharing (2 oxygens shared by 2 silicons)
    if quiet == False:
        progress_bar = tqdm(atoms, desc="Calculating the concentrations Si-Si sites", colour='BLUE', leave=False)
        color_gradient = generate_color_gradient(len(atoms))
        counter = 0
    else:
        progress_bar = atoms
        
    for atom in progress_bar:
        if quiet == False:
            progress_bar.set_description(f"Calculating the concentrations Si-Si sites ...")
            progress_bar.colour = "#%02x%02x%02x" % color_gradient[counter]
            counter += 1
            
        for neighbor in atom.get_neighbours():
            if atom.coordination == 3 and neighbor.coordination == 3:
                dict_concentrations["Si3-Si3"].append(atom.id)
                dict_concentrations["Si3-Si3"].append(neighbor.id)
            if atom.coordination == 3 and neighbor.coordination == 4:
                dict_concentrations["Si3-Si4"].append(atom.id)
                dict_concentrations["Si3-Si4"].append(neighbor.id)
            if atom.coordination == 4 and neighbor.coordination == 4:
                dict_concentrations["Si4-Si4"].append(atom.id)
                dict_concentrations["Si4-Si4"].append(neighbor.id)
            if atom.coordination == 4 and neighbor.coordination == 5:
                dict_concentrations["Si4-Si5"].append(atom.id)
                dict_concentrations["Si4-Si5"].append(neighbor.id)
            if atom.coordination == 5 and neighbor.coordination == 5:
                dict_concentrations["Si5-Si5"].append(atom.id)
                dict_concentrations["Si5-Si5"].append(neighbor.id)
            if atom.coordination == 5 and neighbor.coordination == 6:
                dict_concentrations["Si5-Si6"].append(atom.id)
                dict_concentrations["Si5-Si6"].append(neighbor.id)
            if atom.coordination == 6 and neighbor.coordination == 6:
                dict_concentrations["Si6-Si6"].append(atom.id)
                dict_concentrations["Si6-Si6"].append(neighbor.id)
            if atom.coordination == 6 and neighbor.coordination == 7:
                dict_concentrations["Si6-Si7"].append(atom.id)
                dict_concentrations["Si6-Si7"].append(neighbor.id)
            if atom.coordination == 7 and neighbor.coordination == 7:
                dict_concentrations["Si7-Si7"].append(atom.id)
                dict_concentrations["Si7-Si7"].append(neighbor.id)
            if atom.coordination == 7 and neighbor.coordination == 8:
                dict_concentrations["Si7-Si8"].append(atom.id)
                dict_concentrations["Si7-Si8"].append(neighbor.id)
            if atom.coordination == 8 and neighbor.coordination == 8:
                dict_concentrations["Si8-Si8"].append(atom.id)
                dict_concentrations["Si8-Si8"].append(neighbor.id)
            if atom.coordination == 8 and neighbor.coordination == 9:
                dict_concentrations["Si8-Si9"].append(atom.id)
                dict_concentrations["Si8-Si9"].append(neighbor.id)
            if atom.coordination == 9 and neighbor.coordination == 9:
                dict_concentrations["Si9-Si9"].append(atom.id)
                dict_concentrations["Si9-Si9"].append(neighbor.id)
    
    
    # Calculate the concentrations of each connectivity
    for key, value in dict_concentrations.items():
        dict_concentrations[key] = len(np.unique(value)) / len(atoms) 
    return dict_concentrations