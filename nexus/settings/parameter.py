"""
Module: parameter
-----------------

This module defines classes to represent and handle parameters for the application.

Classes:
--------
    - Parameter: Represents a generic parameter with a name and value.
    - ClusterParameter: Represents a parameter specific to cluster settings.
"""

class Parameter:
    """
    The Parameter class represents a parameter with a name and a value.

    Attributes:
    -----------
        name (str): Name of the parameter.
        value: Value associated with the parameter.

    Methods:
    --------
        __init__(self, name, value): Initializes a Parameter object with a name and value.
        get_name(self): Returns the name of the parameter.
        get_value(self): Returns the value associated with the parameter.
        set_value(self, new_value): Sets a new value for the parameter.
    """

    def __init__(self, name, value) -> None:
        """
        Initializes a Parameter object with a name and value.

        Parameters:
        -----------
            name (str): Name of the parameter.
            value: Value associated with the parameter.
        """
        self.name: str = name
        self.value = value
        self.disable_warnings = False
    
        @property
        def name(self):
            return self.__name
        name.setter
        def name(self, value):
            if not isinstance(value, str):
                raise ValueError(f"Invalid value for 'name': {value}")
            self.__name = value
        
        @property
        def value(self):
            return self.__value
        value.setter
        def value(self, value):
                self.__value = value

    def get_name(self) -> str:
        """
        Returns the name of the parameter.
        
        Returns:
        --------
            str: Name of the parameter.
        """
        return self.name

    def get_value(self):
        """
        Returns the value associated with the parameter.
        
        Returns:
        --------
            value: Value associated with the Parameter.
        """
        return self.value

    def set_value(self, new_value) -> None:
        """
        Sets a new value for the parameter.

        Parameters:
        -----------
            new_value: The new value to be set for the parameter.
        
        Raises:
        -------
            ValueError: If the value provided is invalid for the parameter.

        Note:
        -----
            If the parameter is 'print_clusters_positions', setting it to True will generate a large amount of data.
            This option is not recommended for large systems. If 'print_clusters_positions' is set to True,
            the user will be prompted to confirm the action.
        """
        if self.name == "print_clusters_positions" and self.disable_warnings == False: 
            if new_value == True:
                print(f"\tWARNING: Printing the positions of the clusters will generate a large amount of data.")
                print(f"\t         This option is not recommended for large systems.")
                print(f"\t         Activate the option ? [y/n]")
                response = input()
                if response.lower() == "y":
                    self.value = True
                else:
                    self.value = False
            elif new_value == False:
                self.value = new_value
            else:
                raise ValueError(f"ERROR: Invalid value for 'print_clusters_positions': {new_value}")
        else:
            self.value = new_value
            
class ClusterParameter(Parameter):
    """
    Represents a parameter specific to cluster settings.

    Inherits from Parameter.

    Attributes:
    -----------
        name (str): The name of the cluster parameter.
        value (dict): The value of the cluster parameter, stored as a dictionary.

    Methods:
    --------
        __init__(name, value): Initializes a ClusterParameter object.
        set_cluster_parameter(key, new_value): Replaces a value of the settings.
    """
    def __init__(self, name, value) -> None:
        """
        Initializes a ClusterParameter object.

        Parameters:
        -----------
            name (str): The name of the cluster parameter.
            value (dict): The value of the cluster parameter.
        """
        super().__init__(name, value)
        
    def set_cluster_parameter(self, key, new_value) -> None:
        """
        Replaces a value of the settings.

        Parameters:
        -----------
            key (str): The key of the parameter to be replaced.
            new_value: The new value to be assigned to the parameter.

        Raises:
        -------
            ValueError: If the provided key is not found or the new value is invalid.
        """
        if key not in ["connectivity", "main_former", "criteria", "polyhedra"]:
            raise ValueError(f"\tERROR: Cluster parameter '{key}' couldn't be found in cluster settings.")
        elif key == "criteria":
            if new_value != "distance" and new_value != "bond":
                raise ValueError(f"\tERROR: Cluster parameter '{key}' can only be set as 'distance' or 'bond', not '{new_value}'.")
        elif key == "connectivity":
            if not isinstance(new_value, list):
                raise ValueError(f"\tERROR: Cluster parameter '{key}' can only be set as a list of element (ie ['Si', 'O', 'Si']), not {new_value}.")
        elif key == "polyhedra":
            if not isinstance(new_value, list):
                raise ValueError(f"\tERROR: Cluster parameter '{key}' can only be set as a list of a list of coordinations (ie [[4, 4], [5, 5]]), not {new_value}.")
        for k, v in self.value.items():
            if k == key:
                self.value[key] = new_value
                return



