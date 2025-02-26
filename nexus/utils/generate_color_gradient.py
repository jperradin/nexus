"""
Module: generate_color_gradient
-------------------------------

This module provides a function to generate a color gradient between two colors.

Functions:
----------
    - generate_color_gradient: Generates a color gradient between two colors.
"""

import colorsys

def generate_color_gradient(num_iterations):
    """
    Generates a color gradient between red and blue.

    Parameters:
    -----------
        num_iterations (int): The number of colors to generate in the gradient.

    Returns:
    --------
        list: A list of RGB tuples representing the color gradient.
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