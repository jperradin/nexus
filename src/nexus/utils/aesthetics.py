import colorsys
from collections import OrderedDict
from colorama import Fore, Style

__all__ = [
    'print_title',
    'generate_color_gradient',
    'remove_duplicate_lines'
]


def print_title(__version__) -> None:
    """
    Print the ASCII art title and version to the terminal.

    Args:
        __version__ (str): The version string to display.
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
    Write the ASCII art title and version to a file.

    Args:
        __version__ (str): The version string to write.
        path (str): Path to the output file.
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
    """
    Generate an RGB color gradient for progress bar visualization.

    Interpolates in HSV space from red to blue over the given number of steps.

    Args:
        num_iterations (int): Number of gradient steps to produce.

    Returns:
        list: List of (R, G, B) tuples with integer values in [0, 255].
    """

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
    Remove duplicate lines from a file in place, preserving order.

    Args:
        filepath (str): Path to the file to deduplicate.
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