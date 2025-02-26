"""
Module: write_list_of_files
Description: This module provides a function to write a list of all-in-one unwrapped clusters files to a text file.
"""

import os
from collections import OrderedDict
from natsort import natsorted

def write_list_of_files(dirpath: str) -> None:
    """
    Write a list of the all-in-one unwrapped clusters files to a text file.
    
    Parameters
    ----------
    dirpath : str
        The directory path to search for all-in-one unwrapped clusters files.
        
    Returns
    -------
    None
    """
    
    # Create the list of files
    ordered_files = OrderedDict()
    for root, dirs, files in os.walk(dirpath):
        files = natsorted(files)
        for file in files:
            if file.endswith("all-in-one.xyz"):
                path = os.path.join(root, file)
                ordered_files[path] = None
                
    # Write the list of files to a text file
    with open(os.path.join(dirpath, "list.txt"), "w") as f:
        for file in ordered_files:
            f.write(file + "\n")