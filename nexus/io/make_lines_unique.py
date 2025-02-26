"""
Module: make_lines_unique
Description: This module provides a function to read a file, remove duplicate lines, and rewrite the file with unique lines.
"""

from collections import OrderedDict

def make_lines_unique(filepath: str) -> None:
    """
    Read a file, remove duplicate lines, and rewrite the file with unique lines.
    
    Parameters
    ----------
    filepath : str
        The path to the file to be read and rewritten.
        
    Returns
    -------
    None
    """
    
    # Read the file and store unique lines in an OrderedDict
    unique_lines = OrderedDict()
    with open(filepath, 'r') as file:
        for line in file:
            unique_lines[line] = None

    # Rewrite the unique lines back to the file
    with open(filepath, 'w') as file:
        for line in unique_lines:
            file.write(line)