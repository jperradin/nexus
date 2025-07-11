import colorsys
from collections import OrderedDict
from colorama import Fore, Style

__all__ = [
    'print_title',
    'generate_color_gradient',
    'remove_duplicate_lines'
]


"""
Module: print_title
-------------------

This module provides a function to print the title and version of the package.

Functions:
----------
    - print_title: Prints the title and the version of the package.
"""

def print_title(__version__) -> None:
    """
    Prints the title and the version of the package.

    Parameters:
    -----------
        __version__ (str): The version of the package.

    Returns:
    --------
        None
    """
    title = Fore.LIGHTBLUE_EX + r"""
                                                                
    :-:      --:   -=+++=-  -:     ==: ::       :-:    :+*##*=: 
   *@@%#-   +@@# -%@@@@@@@=*@%*: :#@@+=@@+     -@@%:  =@@@@@@@+ 
  -@@@@@@+  #@@%:%@@%+==== -%@@@#%@@#:+@@#     =@@@-  %@@%--=-  
  =@@@@@@@= #@@#-@@@%+=-    :*@@@@@+  +@@%:    +@@@-  +@@@@%#*: 
  =@@@=#@@@#@@@+=@@@@@@*      #@@@#   =@@@=    *@@%:   -+#%@@@%-
  =@@% :%@@@@@@==@@@*--     :*@@@@@#: :%@@@*==*@@@+ -##+  :#@@@+
  =@@%  -%@@@@# :@@@#*###*:-%@@@#%@@@- -%@@@@@@@@*  *@@@#*#@@@@-
  -%%+   :+##+:  =%@@@@@@%:=@@%- :#%%-  :+#%%%%*-   :%@@@@@@@#= 
    :              :-----:  :-     :       :::        -=+++=:   
    """ + Style.RESET_ALL
    print(title)
    print(f"__version__ \u279c  {__version__}\n")
    return

def print_title_to_file(__version__, path) -> None:
    """
    Prints the title and the version of the package.

    Parameters:
    -----------
        __version__ (str): The version of the package.

    Returns:
    --------
        None
    """
    title = r"""
                                                                
    :-:      --:   -=+++=-  -:     ==: ::       :-:    :+*##*=: 
   *@@%#-   +@@# -%@@@@@@@=*@%*: :#@@+=@@+     -@@%:  =@@@@@@@+ 
  -@@@@@@+  #@@%:%@@%+==== -%@@@#%@@#:+@@#     =@@@-  %@@%--=-  
  =@@@@@@@= #@@#-@@@%+=-    :*@@@@@+  +@@%:    +@@@-  +@@@@%#*: 
  =@@@=#@@@#@@@+=@@@@@@*      #@@@#   =@@@=    *@@%:   -+#%@@@%-
  =@@% :%@@@@@@==@@@*--     :*@@@@@#: :%@@@*==*@@@+ -##+  :#@@@+
  =@@%  -%@@@@# :@@@#*###*:-%@@@#%@@@- -%@@@@@@@@*  *@@@#*#@@@@-
  -%%+   :+##+:  =%@@@@@@%:=@@%- :#%%-  :+#%%%%*-   :%@@@@@@@#= 
    :              :-----:  :-     :       :::        -=+++=:   
    """
    with open(path, 'w') as f:
        f.write(title)
        f.write("\n")
        f.write(f"__version__ \u279c  {__version__}\n")
    return

def generate_color_gradient(num_iterations):
    """ Generate a gradient of colors to update at each tqdm iteration """

    # Define the start and end colors in RGB
    start_color = (255, 0, 0)  # Red
    end_color = (0, 0, 255)    # Blue
    
    # Check if num_iterations is 0
    if num_iterations == 0:
        return [start_color]  # Return a list with only the start color

    # Check if num_iterations is 1
    elif num_iterations == 1:
        return [start_color, end_color]  # Return a list containing both start and end colors
    else:
        num_iterations += 1 
            
    # Convert RGB to HSV
    start_hsv = colorsys.rgb_to_hsv(*[x / 255.0 for x in start_color])
    end_hsv = colorsys.rgb_to_hsv(*[x / 255.0 for x in end_color])

    # Interpolate between the start and end colors
    color_gradient = []
    for i in range(num_iterations):
        ratio = i / (num_iterations - 1)
        hsv = (
            start_hsv[0] + ratio * (end_hsv[0] - start_hsv[0]),
            start_hsv[1] + ratio * (end_hsv[1] - start_hsv[1]),
            start_hsv[2] + ratio * (end_hsv[2] - start_hsv[2])
        )
        rgb = tuple(int(x * 255) for x in colorsys.hsv_to_rgb(*hsv))
        color_gradient.append(rgb)

    return color_gradient

def remove_duplicate_lines(filepath: str) -> None:
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