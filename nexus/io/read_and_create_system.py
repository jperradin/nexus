"""
Module: read_and_create_system
Description: This module provides functions to read an XYZ file and create a System object from the data.
"""

# external imports
from tqdm import tqdm
import numpy as np
import importlib

# internal imports
from ..core.atom import Atom
from ..core.system import System
from ..data import chemical_symbols

def seek_to_line(file, line_number) -> None:
    """
    Seeks to the specified line number in the file.

    Parameters
    ----------
    file : file object
        The file object.
    line_number : int
        The line number to seek to.
    
    Returns
    -------
    None
    """
    
    if line_number == 0:
        file.readline()  # Skip the first line
        return
    
    file.seek(0)  # Go to the beginning of the file
    
    current_line = 0
    
    # Iterate through the file until the desired line is reached
    while current_line != line_number + 1:
        file.readline()
        current_line += 1
    
    return

def read_and_create_system(file_path, frame, frame_size, settings, cutoffs, start, end) -> System:
    """
    Read the XYZ file and return the frame as a System object.
    NOTE: this function is extension dependent.
    
    Parameters
    ----------
    file_path : str
        Path to the XYZ file.
    frame : int
        Frame number to read.
    frame_size : int
        Number of atoms in the frame + number of header lines.
    settings : Settings
        Settings object.
    cutoffs : dict
        Dictionary with the cutoffs for each pair of elements.
    start : int
        ID of the first frame to read.
    end : int
        ID of the last frame to read.
        
    Returns
    -------
    System
        The system object created with the information provided by the input file.
    """
    
    # import extension
    extension = settings.extension.get_value()
    module = importlib.import_module(f"nexus.extensions.{extension}")
    
    system = System(settings)
    
    header = settings.header.get_value()
    
    # Open the file
    with open(file_path, "r") as f:
        seek_to_line(f, frame * frame_size)  # Go to the beginning of the frame
        
        jump = f.readline()  # Skip the comment line
        
        atom_skipped = {}
        sum_skipped = 0
        
        # Read the atoms coordinates in the frame
        if not settings.quiet.get_value():
            progress_bar = tqdm(range(frame_size - header), desc=f"Reading frame {frame}", unit="atoms", leave=False, colour="BLUE")
        else:
            progress_bar = range(frame_size - header)
            
        for i in progress_bar:
            
            line = f.readline()
            
            parts = line.split()  # line is like: "Si 1.234 5.678 9.101"
            
            element = parts[0]
            
            if element in module.LIST_OF_SUPPORTED_ELEMENTS and element in chemical_symbols:
            
                x = float(parts[1])
                y = float(parts[2])
                z = float(parts[3])
                
                # Create the Atom object with the current information
                position = np.array([x, y, z])

                current_atom = Atom(element, i, position, frame, cutoffs, extension=settings.extension.get_value())
                
                # Add the atom to the system
                system.add_atom(current_atom)
            elif element not in module.LIST_OF_SUPPORTED_ELEMENTS and element in chemical_symbols:
                sum_skipped += 1
                if element not in atom_skipped:
                    atom_skipped[element] = 1
                else:
                    atom_skipped[element] += 1
            
    f.close()

    # Check if all the atoms were read
    if len(system.get_atoms()) + sum_skipped != settings.number_of_atoms.get_value():
        raise ValueError(f"\tFrame {frame} does not have the expected number of atoms. Expected: {frame_size - header}, got: {len(system.get_atoms())} stored + {sum_skipped} skipped.")
    
    if len(atom_skipped) > 0:
        expd = settings._output_directory
        if frame == start:
            # Write the header of the file if it is the first frame.
            with open(f"{expd}/skipped_atoms.log", "w") as f:
                f.write(f"Extension: '{extension}'\n")
                f.write(f"Frame: {frame}\n")
                f.write(f"Total number of atoms: {frame_size - header}\n")
                f.write(f"Number of atoms skipped: {len(atom_skipped)}\n")
                for k, v in atom_skipped.items():
                    f.write(f"\u279c {k} : {v}\n")
        else:
            with open(f"{expd}/skipped_atoms.log", "a") as f:
                f.write(f"Extension: '{extension}'\n")
                f.write(f"Frame: {frame}\n")
                f.write(f"Total number of atoms: {frame_size - header}\n")
                f.write(f"Number of atoms skipped: {len(atom_skipped)}\n")
                for k, v in atom_skipped.items():
                    f.write(f"\u279c {k} : {v}\n")

    # End reading the file and return the System object
    return system