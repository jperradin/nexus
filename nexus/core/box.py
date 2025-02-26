"""
Box module for the Nexus project.

This module defines the Box class, which represents a simulation box in three-dimensional space.

Classes:
    - Box: Represents a simulation box in three-dimensional space.
"""

class Box:
    """
    Represents a simulation box in three-dimensional space at each frame of the trajectory.

    Attributes:
        length_x (list): Length of the box in the x-direction.
        length_y (list): Length of the box in the y-direction.
        length_z (list): Length of the box in the z-direction.
        volume (list): Volume of the box.

    Methods:
        __init__: Initializes a Box object.
        add_box: Adds a box to the list of boxes.
        get_volume: Calculates and returns the volume of the box.
        get_box_dimensions: Returns the dimensions of the box.
    """

    def __init__(self) -> None:
        """
        Initializes a Box object.
        """
        self.length_x : list = []       # list of component x of the simulation box size
        self.length_y : list = []       # list of component y of the simulation box size
        self.length_z : list = []       # list of component z of the simulation box size
        self.volume : list = []         # list of volume of the simulation box size
        
    def add_box(self, length_x, length_y, length_z) -> None:
        """
        Adds a box to the list of boxes.

        Parameters
        ----------
            - length_x (float) :Length of the box in the x-direction.
            - length_y (float) :Length of the box in the y-direction.
            - length_z (float) :Length of the box in the z-direction.
        
        Returns
        -------
            - None
        """
        self.length_x.append(length_x)
        self.length_y.append(length_y)
        self.length_z.append(length_z)
        self.volume.append(self.get_volume(len(self.length_x) - 1))

    def get_volume(self, configuration) -> float:
        """
        Calculates and returns the volume of the box assuming the box is always cubic.

        Parameters:
        -----------
            - configuration (int): Index of the configuration.

        Returns:
        --------
            - float: Volume of the box.
        """
        return self.length_x[configuration] * self.length_y[configuration] * self.length_z[configuration]

    def get_box_dimensions(self, configuration) -> list:
        """
        Returns the dimensions of the box.

        Parameters:
        -----------
            - configuration (int): Index of the configuration.

        Returns:
        --------
            - list: Dimensions of the box [length_x, length_y, length_z].
        """
        return [self.length_x[configuration], self.length_y[configuration], self.length_z[configuration]]
