"""
This file contains all the methods / functions that are specific to Na-Na clusters.
"""

# external imports
import numpy as np

# internal imports
from ..core.atom import Atom

# List of supported elements for the extension Na
LIST_OF_SUPPORTED_ELEMENTS = ["Na", "O"]

class Sodium(Atom):
    def __init__(self, element, id, position, frame, cutoffs, extension) -> None:
        super().__init__(element, id, position, frame, cutoffs, extension)
    
    def calculate_coordination(self) -> int:
        """
        Calculate the coordination number of the atom (ie the number of first neighbours) for the extension Na
        """
        self.coordination = len([neighbour for neighbour in self.neighbours if neighbour.get_element() == "O"])

class Oxygen(Atom):
    def __init__(self, element, id, position, frame, cutoffs, extension) -> None:
        super().__init__(element, id, position, frame, cutoffs, extension)
        
    def calculate_coordination(self) -> int:
        """
        Calculate the coordination number of the atom (ie the number of first neighbours) for the extension Na
        """
        self.coordination = len([neighbour for neighbour in self.neighbours if neighbour.get_element() == "Na"])

def transform_into_subclass(atom:Atom) -> object:
    """
    Return a Sodium object or Oxygen object from the subclass Sodium or Oxygen whether the atom.element is 'Na' or 'O'.  
    """
    if atom.get_element() == 'O':
        return Oxygen(atom.element, atom.id, atom.position, atom.frame, atom.cutoffs, atom.extension)
    elif atom.get_element() == 'Na':
        return Sodium(atom.element, atom.id, atom.position, atom.frame, atom.cutoffs, atom.extension)
    else:
        raise ValueError(f"\tERROR: Atom {atom.element} - {atom.id} can be transformed into Sodium or Oxygen object.")

def get_connectivity(cluster_settings) -> list:
    """
    Return the connectivity between the atoms.
    """
    return ['Na-Na']
    