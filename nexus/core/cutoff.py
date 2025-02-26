"""
Cutoff module for the Nexus project.

This module defines the Cutoff class, which manages cutoff distances for pairs of elements.

Classes:
    - Cutoff: Manages cutoff distances for pairs of elements.
"""

class Cutoff:
    """
    Manages cutoff distances for pairs of elements.

    Attributes:
        cutoffs (dict): Dictionary containing the cutoffs for each pair of elements.
        pairs (list): List of pairs of elements.
        values (list): List of cutoff values.

    Methods:
        __init__: Initializes a Cutoff object.
        get_cutoff: Returns the cutoff for the pair of elements.
        get_max_cutoff: Returns the maximum cutoff in the system.
    """
    
    def __init__(self, cutoffs) -> None:
        """
        Initializes the Cutoff object.

        Parameters:
            cutoffs (dict): Dictionary containing the cutoffs for each pair of elements.
        """
        self.cutoffs : dict = cutoffs   # dictionnary of the cutoffs 
        self.pairs : list = []          # list of the pairs of atoms
        self.values : list = []         # list of the cutoff values
        
        self.load_cutoffs()
        
    def load_cutoffs(self) -> None:
        r"""
        Loads the cutoff values with their associated pair.

        Returns:
        --------
            - None.
        """
        for cutoff in self.cutoffs:
            self.pairs.append([cutoff['element1'], cutoff['element2']])
            self.values.append(cutoff['value'])
        
    def get_cutoff(self, element1, element2) -> float:
        """
        Returns the cutoff for the pair of elements.
        
        Parameters:
        -----------
        - element1 (str): First element.
        - element2 (str): Second element.
        
        Returns:
        --------
        - float: Cutoff for the pair of elements.
        """
        try:
            index = self.pairs.index([element1, element2])
        except:
            index = self.pairs.index([element2, element1])
            
        return self.values[index]
    
    def get_max_cutoff(self) -> float:
        """
        Returns the maximum cutoff in the system.
        
        Returns:
        --------
        - float: Maximum cutoff in the system.
        """
        return max(self.values)
