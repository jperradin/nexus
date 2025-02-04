# internal imports
from .parameter import Parameter
from .parameter import ClusterParameter

# external imports
import importlib


class Settings:
    """
    Represents the settings for a project.

    Attributes:
    ----------
        - project_name (Parameter): Name of the project used for the output directory.
        - export_directory (Parameter): Parent directory where the output files will be saved.
        - build_fancy_recaps (Parameter): Whether to build fancy recaps of the results into a single file.
        - build_fancy_plots (Parameter): Whether to build fancy plots of the results into a single file.
        - path_to_xyz_file (Parameter): Path to the XYZ file containing the atomic coordinates.
        - number_of_atoms (Parameter): Number of atoms in the system.
        - number_of_frames (Parameter): Number of frames in the XYZ file.
        - header (Parameter): Number of lines in the header of the XYZ file.
        - range_of_frames (Parameter): Range of frames to be analysed.
        - frames_to_analyse (Parameter): Frames to be analysed.
        - timestep (Parameter): Timestep of the simulation.
        - lbox (Parameter): Box length.
        - temperature (Parameter): Temperature of the system.
        - pressure (Parameter): Pressure of the system.
        - version (Parameter): Version of the software.
        - quiet (Parameter): Whether to print settings or not.
        - overwrite_results (Parameter): Whether to overwrite files by default.
        - print_clusters_positions (Parameter): Whether to print the positions of the clusters.
        - max_size (Parameter): Maximum size of the clusters for the cluster size distribution.
        - supported_extensions (Parameter): List of supported extensions.

    Methods:
    --------
        - __init__(extension="SiOz"): Initializes a Settings object with default settings.
        - load_default_settings(extension): Loads default settings based on the extension.
        - print_settings(): Prints the current settings.
        - print_all_settings(): Prints all settings, including those not recommended for printing.
    """

    def __init__(self, extension="SiOz") -> None:
        """
        Initializes a Settings object with default settings.

        Parameters:
        -----------
            - extension (str): The extension of the project.
        """
        self.load_default_settings(extension)

    def load_default_settings(self, extension) -> None:
        """
        Loads default settings based on the extension.

        Parameters:
        -----------
            - extension (str): The extension of the project.
        """
        self.project_name : Parameter = Parameter("project_name", "default")                        # Name of the project that will be used for the output directory
        self.export_directory : Parameter = Parameter("export_directory", "export")                 # Parent directory where the output files will be saved in the project directory
        self._output_directory : str = ''                                                           # Output directory (hidden attribute)
        self.build_fancy_recaps : Parameter = Parameter("build_fancy_recaps", False)                # Build fancy recaps of the results into a single file
        self.build_fancy_plots : Parameter = Parameter("build_fancy_plots", False)                  # Build fancy plots of the results into a single file
        self.path_to_xyz_file : Parameter = Parameter("path_to_xyz_file", "input.xyz")              # Path to the XYZ file containing the atomic coordinates
        self.number_of_atoms : Parameter = Parameter("number_of_atoms", 0)                          # Number of atoms in the system
        self.number_of_frames : Parameter = Parameter("number_of_frames", 0)                        # Number of frames in the XYZ file
        self.header : Parameter = Parameter("header", 0)                                            # Number of lines in the header of the XYZ file
        self.range_of_frames : Parameter = Parameter("range_of_frames", None)                       # Range of frames to be analysed
        self.frames_to_analyse : Parameter = Parameter("frames_to_analyse", None)                   # Frames to be analysed
        self.timestep : Parameter = Parameter("timestep", 0.0016)                                   # Timestep of the simulation
        self.lbox : Parameter = Parameter("lbox", 0.0)                                              # Box length
        self.temperature : Parameter = Parameter("temperature", 0.0)                                # Temperature of the system
        self.pressure : Parameter = Parameter("pressure", 0.0)                                      # Pressure of the system
        self.version : Parameter = Parameter("version", "0.1.11")                                    # Version of the software
        self.quiet : Parameter = Parameter("quiet", False)                                          # Do not print any settings
        self.overwrite_results : Parameter = Parameter("overwrite_results", False)                  # Overwrite files by default
        self.print_clusters_positions : Parameter = Parameter("print_clusters_positions", False)    # Print the positions of the clusters
        self.max_size : Parameter = Parameter("max_size", 100)                                      # Maximum size of the clusters for the cluster size distribution

        self.supported_extensions : Parameter = Parameter("extensions", ["SiOz", "SiSi", "OO"]) # Update the list when adding a new extension

        if extension in self.supported_extensions.get_value():
            module = importlib.import_module(f"nexus.extensions.{extension}")
            default_settings = module.get_default_settings()
            self.extension : Parameter = default_settings["extension"]                          # Name of the extension used for the analysis (default is 'SiOz')
            self.structure : Parameter = default_settings["structure"]                          # Detailed chemical composition of the structure of the system.
            self.cluster_settings : ClusterParameter = default_settings["cluster_settings"]     # All the settings for the cluster analysis
            self.cutoffs : Parameter = default_settings["cutoffs"]                              # All the cutoffs of each pair of atoms.
        else:
            raise ValueError(f"Extension '{extension}' is not supported.")

        self.settings_to_print : str = '' # Settings to print on the terminal.

    def print_settings(self) -> None:
        """
        Prints the current settings.
        """
        max_attr_length = max(len(attr) for attr in self.__dict__)
        separator = "\t\t________________________________________________"
        settings_output = f"\tSETTINGS:\n{separator}\n"
        max_attr_length = max(len('Path to input file'), len('Number of frames'))
        settings_output += f"\t\t{'Path to input file'.ljust(max_attr_length)} \u279c\t {self.path_to_xyz_file.get_value()}\n"
        settings_output += f"\t\t{'Number of frames'.ljust(max_attr_length)} \u279c\t {self.number_of_frames.get_value()}\n"
        if self.range_of_frames.get_value() is not None:
            settings_output += f"\t\tRange of frames    \u279c\t {self.range_of_frames.get_value()}\n"
        settings_output += f"{separator}\n"
        settings_output += f"\t\tStructure:\n"
        max_attr_length = max(len("Number of atoms"), len("Species"))
        settings_output += f"\t\t  {'Number of atoms'.ljust(int(max_attr_length/2))} \u279c\t {self.number_of_atoms.get_value()}\n"
        for atom in self.structure.get_value():
            settings_output += f"\t\t  {'Species'.ljust(max_attr_length)} \u279c\t {atom['element']:2}\t|\tNumber of atoms \u279c\t {atom['number']}\n"
        settings_output += f"{separator}\n"
        settings_output += f"\t\tExport directory   \u279c\t {self._output_directory}\n"
        settings_output += f"{separator}\n"
        settings_output += f"\t\tCluster settings:\n"
        max_attr_length = max(len(str(k)) for k,v in self.cluster_settings.get_value().items())
        for k, v in self.cluster_settings.get_value().items():
            settings_output += f"\t\t  {k.ljust(max_attr_length)} \u279c\t {v}\n"

        settings_output += "\n"
        if self.quiet.get_value() == False:
            print("\r" + settings_output, end="")
            self.settings_to_print = settings_output
        else:
            self.settings_to_print = settings_output


    def print_all_settings(self) -> None:
        """
        Prints all settings, including those not recommended for printing.
        """
        max_attr_length = max(len(attr) for attr in self.__dict__)
        separator = "\t\t________________________________________________"
        print(f"\tSETTINGS:")
        print(separator)
        for p, v in self.__dict__.items():
            print(f"\t\t{v.get_name().ljust(max_attr_length)} ➜ {v.get_value()}")
