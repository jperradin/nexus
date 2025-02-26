"""
Module: result
--------------

This module defines various classes to represent and handle different types of results
related to cluster size analysis. Each class provides methods to calculate specific 
metrics and write the results to files.

Classes:
--------
    - Result: Base class for all result types.
    - AverageClusterSize: Represents a result for average cluster sizes.
    - BiggestClusterSize: Represents a result for the biggest cluster size.
    - SpanningClusterSize: Represents a result for the spanning cluster size.
    - ClusterSizeDistribution: Represents a result for the cluster size distribution.
    - GyrationRadiusDistribution: Represents a result for the gyration radius distribution.
    - CorrelationLength: Represents a result for the correlation length.
    - OrderParameter: Represents a result for the order parameter.
    - PercolationProbability: Represents a result for the percolation probability.
"""

# External imports
import numpy as np
import os

# Internal imports
from .make_lines_unique import make_lines_unique

__all__ = [
    'Result',
    'AverageClusterSize',
    'BiggestClusterSize',
    'SpanningClusterSize',
    'ClusterSizeDistribution',
    'GyrationRadiusDistribution',
    'CorrelationLength',
    'OrderParameter',
    'PercolationProbability'
]

class Result:
    """
    Represents a generic result.

    Attributes:
    -----------
        property (str): The property name.
        info (str): Additional information about the result.
        init_frame (int): The initial frame.
        timeline (list): List of data points over time.
        result (float): The computed result.
        error (float): The error associated with the result.
    """
    def __init__(self, name, info, init_frame) -> None:
        """
        Initializes a Result object.

        Parameters:
        -----------
            name (str): The property name.
            info (str): Additional information about the result.
            init_frame (int): The initial frame.
        """
        self.property : str = name
        self.info : str = info
        self.init_frame : int = init_frame
        self.timeline : list = []
        self.timeline_c : list = []
        self.result : float = 0.
        self.error : float = 0.
        self.concentration : float = 0.

class AverageClusterSize(Result):
    """
    Represents a result for average cluster sizes.

    Attributes:
    -----------
        average_size (float): The computed average cluster size.
        filepath (str): The path to the output file.
    """
    def __init__(self, name: str, info: str, init_frame: int) -> None:
        """
        Initializes an AverageClusterSizeResult object.

        Parameters:
        -----------
            name (str): The property name.
            info (str): Additional information about the result.
        """
        super().__init__(name, info, init_frame)
        self.average_size = 0
        self.filepath = None
    
    def add_to_timeline(self, value, concentration) -> None:
        """
        Appends a data point to the timeline.

        Parameters:
        -----------
            value: The data point to add.
            concentration: The concentration associated with the data point.
        """
        self.timeline.append(value)
        self.timeline_c.append(concentration)
    
    def calculate_average_cluster_size(self) -> None:
        """
        Calculates the average cluster size based on the timeline data.
        """
        self.timeline = np.array(self.timeline)
        self.average_size = np.zeros(self.timeline.shape[0])
        for f in range(self.timeline.shape[0]):
            dict_sizes = self.timeline[f]
            list_sizes = dict_sizes[f+self.init_frame]
            _bins = list(range(min(list_sizes), max(list_sizes)+2))
            _hist,_bins = np.histogram(list_sizes, bins=_bins)
            _bins = _bins[:-1]
            _idx = _hist == 0
            appearences = _hist[~_idx]
            sizes = _bins[~_idx]
            
            normalization = np.sum(sizes * appearences)
            if normalization == 0: normalization = 1 # avoid zero division
            for s, ns in zip(sizes, appearences):
                self.average_size[f] += (s**2 * ns) / normalization
        
        self.result = np.mean(self.average_size)
        self.error = np.std(self.average_size) / np.sqrt(len(self.average_size))
        self.concentration = np.mean(self.timeline_c)

    def write_file_header(self, overwrite, path_to_directory, number_of_frames) -> None:
        """
        Initializes the output file with a header.

        Parameters:
        -----------
            overwrite (bool): Whether to overwrite the existing file.
            path_to_directory (str): The directory where the output file will be saved.
            number_of_frames (int): The number of frames used in averaging.
        """
        filename = "average_cluster_size.dat"
        if not os.path.exists(path_to_directory):
            os.makedirs(path_to_directory)

        self.filepath = os.path.join(path_to_directory, filename)
        
        if not overwrite and os.path.exists(self.filepath):
            with open(self.filepath, 'a', encoding='utf-8') as output:
                output.write(f"# Average cluster size \u279c {number_of_frames} frames averaged.\n")
                output.write("# Concentration \u279c Average cluster size +/- Error\n")
            output.close()
        else:
            with open(self.filepath, 'w', encoding='utf-8') as output:
                output.write(f"# Average cluster size \u279c {number_of_frames} frames averaged.\n")
                output.write("# Concentration \u279c Average cluster size +/- Error\n")
            output.close()

    def append_results_to_file(self) -> None:
        """
        Appends the result to the output file.
        """
        with open(self.filepath, 'a', encoding='utf-8') as output:
            output.write(f"{self.concentration:10.6f} \u279c {self.result:10.6f} +/- {self.error:<10.5f} # {self.info}\n")
        output.close()
        
        make_lines_unique(self.filepath)

class BiggestClusterSize(Result):
    """
    Represents a result for the biggest cluster size.

    Attributes:
    -----------
        biggest_size (float): The computed biggest cluster size.
        filepath (str): The path to the output file.
    """
    def __init__(self, name: str, info: str, init_frame: int) -> None:
        """
        Initializes a BiggestClusterSizeResult object.

        Parameters:
        -----------
            name (str): The property name.
            info (str): Additional information about the result.
        """
        super().__init__(name, info, init_frame)
        self.biggest_size = 0
        self.filepath = None

    def add_to_timeline(self, value, concentration) -> None:
        """
        Appends a data point to the timeline.

        Parameters:
        -----------
            value: The data point to add.
            concentration: The concentration associated with the data point.
        """
        self.timeline.append(value)
        self.timeline_c.append(concentration)

    def calculate_biggest_cluster_size(self) -> None:
        """
        Calculates the biggest cluster size based on the timeline data.
        """
        self.timeline = np.array(self.timeline)
        self.biggest_size = np.zeros(self.timeline.shape[0])
        for f in range(self.timeline.shape[0]):
            dict_sizes = self.timeline[f]
            list_sizes = dict_sizes[f+self.init_frame]
            self.biggest_size[f] = max(list_sizes)
        
        self.result = np.mean(self.biggest_size)
        self.error  = np.std(self.biggest_size) / np.sqrt(len(self.biggest_size))
        self.concentration = np.mean(self.timeline_c)
        
    def write_file_header(self, overwrite, path_to_directory, number_of_frames) -> None:
        """
        Initializes the output file with a header.

        Parameters:
        -----------
            overwrite (bool): Whether to overwrite the existing file.
            path_to_directory (str): The directory where the output file will be saved.
            number_of_frames (int): The number of frames used in averaging.
        """
        filename = "biggest_cluster_size.dat"
        if not os.path.exists(path_to_directory):
            os.makedirs(path_to_directory)

        self.filepath = os.path.join(path_to_directory, filename)
        
        if not overwrite and os.path.exists(self.filepath):
            with open(self.filepath, 'a', encoding='utf-8') as output:
                output.write(f"# Biggest cluster size \u279c {number_of_frames} frames averaged.\n")
                output.write("# Concentration \u279c Biggest cluster size +/- Error\n")
            output.close()
        else:
            with open(self.filepath, 'w', encoding='utf-8') as output:
                output.write(f"# Biggest cluster size \u279c {number_of_frames} frames averaged.\n")
                output.write("# Concentration \u279c Biggest cluster size +/- Error\n")
            output.close()
        
    def append_results_to_file(self):
        """
        Appends the result to the output file.
        """
        with open(self.filepath, 'a', encoding='utf-8') as output:
            output.write(f"{self.concentration:10.6f} \u279c {self.result:10.5f} +/- {self.error:<10.5f} # {self.info}\n")
        output.close()
        
        make_lines_unique(self.filepath)
        
class SpanningClusterSize(Result):
    """
    Represents a result for the spanning cluster size.

    Attributes:
    -----------
        spanning_size (float): The computed spanning cluster size.
        filepath (str): The path to the output file.
    """
    def __init__(self, name: str, info: str, init_frame: int) -> None:
        """
        Initializes a SpanningClusterSizeResult object.

        Parameters:
        -----------
            name (str): The property name.
            info (str): Additional information about the result.
        """
        super().__init__(name, info, init_frame)
        self.spanning_size = 0
        self.filepath = None

    def add_to_timeline(self, value, concentration) -> None:
        """
        Appends a data point to the timeline.

        Parameters:
        -----------
            value: The data point to add.
            concentration: The concentration associated with the data point.
        """
        self.timeline.append(value)
        self.timeline_c.append(concentration)

    def calculate_spanning_cluster_size(self) -> None:
        """
        Calculates the spanning cluster size based on the timeline data.
        """
        self.timeline = np.array(self.timeline)
        self.spanning_size = np.zeros(self.timeline.shape[0])
        for f in range(self.timeline.shape[0]):
            dict_sizes = self.timeline[f]
            list_sizes = dict_sizes[f+self.init_frame]
            self.spanning_size[f] = max(list_sizes)
        
        self.result = np.mean(self.spanning_size)
        self.error  = np.std(self.spanning_size) / np.sqrt(len(self.spanning_size))
        self.concentration = np.mean(self.timeline_c)
        
    def write_file_header(self, overwrite, path_to_directory, number_of_frames) -> None:
        """
        Initializes the output file with a header.

        Parameters:
        -----------
            overwrite (bool): Whether to overwrite the existing file.
            path_to_directory (str): The directory where the output file will be saved.
            number_of_frames (int): The number of frames used in averaging.
        """
        filename = "spanning_cluster_size.dat"
        if not os.path.exists(path_to_directory):
            os.makedirs(path_to_directory)

        self.filepath = os.path.join(path_to_directory, filename)
        
        if not overwrite and os.path.exists(self.filepath):
            with open(self.filepath, 'a', encoding='utf-8') as output:
                output.write(f"# Spanning cluster size \u279c {number_of_frames} frames averaged.\n")
                output.write("# Concentration \u279c Spanning cluster size +/- Error\n")
            output.close()
        else:
            with open(self.filepath, 'w', encoding='utf-8') as output:
                output.write(f"# Spanning cluster size \u279c {number_of_frames} frames averaged.\n")
                output.write("# Concentration \u279c Spanning cluster size +/- Error\n")
            output.close()
    
    def append_results_to_file(self):
        """
        Appends the result to the output file.
        """
        with open(self.filepath, 'a', encoding='utf-8') as output:
            output.write(f"{self.concentration:10.6f} \u279c {self.result:10.5f} +/- {self.error:<10.5f} # {self.info}\n")
        output.close()
        
        make_lines_unique(self.filepath)
    
class ClusterSizeDistribution(Result):
    """
    Represents a result for the cluster size distribution.
    
    Attributes:
    -----------
        distribution (dict): The computed cluster size distribution.
        filepath (str): The path to the output file.
    """
    def __init__(self, name: str, info: str, init_frame: int) -> None:
        """
        Initializes a ClusterSizeDistributionResult object.

        Parameters:
        -----------
            name (str): The property name.
            info (str): Additional information about the result.
            init_frame (int): The initial frame.
        """
        super().__init__(name, info, init_frame)
        self.distribution = {}
        self.timeline = {}
        self.filepath = None
        
    def add_to_timeline(self, frame: int, value: list, concentration: float) -> None:
        """
        Appends a data point to the timeline.

        Parameters:
        -----------
            frame (int): The frame number.
            value (list): The data point to add.
            concentration (float): The concentration associated with the data point.
        """
        self.timeline[frame] = value
        self.timeline_c.append(concentration)
    
    def calculate_cluster_size_distribution(self) -> None:
        """
        Calculates the cluster size distribution based on the timeline data.
        """
        # key 1 is the frame number
        # key 2 is the size of the cluster
        for frame, sizes in self.timeline.items():
            for size, ns in sizes.items():
                if size not in self.distribution:
                    self.distribution[size] = []
                if isinstance(ns, list):
                    for n in ns:
                        self.distribution[size].append(n)
                else:
                    self.distribution[size].append(ns)
        
        # Sort the distribution by size decreasingly
        self.distribution = dict(sorted(self.distribution.items(), key=lambda item: item[0], reverse=True))
        self.concentration = np.mean(self.timeline_c)
                            
    def write_file_header(self, path_to_directory, number_of_frames) -> None:
        """
        Initializes the output file with a header.

        Parameters:
        -----------
            path_to_directory (str): The directory where the output file will be saved.
            number_of_frames (int): The number of frames used in averaging.
        """
        filename = f"cluster_size_distribution-{self.info}.dat"
        if not os.path.exists(path_to_directory):
            os.makedirs(path_to_directory)

        self.filepath = os.path.join(path_to_directory, filename)
        
        with open(self.filepath, 'w', encoding='utf-8') as output:
            output.write(f"# Cluster size distribution \u279c {number_of_frames} frames analysed.\n")
            output.write("# Size \u279c Number of clusters\n")
        output.close()
    
    def append_results_to_file(self):
        """
        Appends the result to the output file.
        """
        
        # Get largest list of sizes
        max_length = 0
        for size, ns in self.distribution.items():
            if len(ns) > max_length:
                max_length = len(ns)
        
        # Fill the lists with zeros
        for size, ns in self.distribution.items():
            if len(ns) < max_length:
                self.distribution[size] += [0] * (max_length - len(ns))
                
        with open(self.filepath, 'a', encoding='utf-8') as output:
            output.write(f"# Concentration \u279c {self.concentration}\n")
            for size, ns in self.distribution.items():
                # output.write(f"# size \u279c {size}\n")
                output.write(f"{size:7d} \u279c {np.sum(ns):<5d} # +/- {np.std(ns)/np.sqrt(len(ns)):5.5f}\n")

class GyrationRadiusDistribution(Result):
    """
    Represents a result for the gyration radius distribution.
    
    Attributes:
    -----------
        distribution (dict): The computed gyration radius distribution.
        filepath (str): The path to the output file.
    """
    def __init__(self, name, info, init_frame) -> None:
        """
        Initializes a GyrationRadiusDistributionResult object.
        
        Parameters:
        -----------
            name (str): The property name.
            info (str): Additional information about the result.
            init_frame (int): The initial frame.
        """
        super().__init__(name, info, init_frame)
        self.distribution_rgyr = {}
        self.timeline = {}
        self.filepath_rgyr = None
        
    def add_to_timeline(self, frame: int, value: list, concentration: float) -> None:
        """
        Appends a data point to the timeline.

        Parameters:
        -----------
            frame (int): The frame number.
            value (list): The data point to add.
            concentration (float): The concentration associated with the data point.
        """
        self.timeline[frame] = value
        self.timeline_c.append(concentration)
    
    def calculate_gyration_radius_distribution(self) -> None:
        """
        Calculates the gyration radius distribution based on the timeline data.
        """
        for frame, sizes in self.timeline.items():
            for size, radii in sizes.items():
                if size not in self.distribution_rgyr:
                    self.distribution_rgyr[size] = []
                if radii != []:
                    for r in radii:
                        self.distribution_rgyr[size].append(r)
        
        # Sort the distribution_rgyr by size decreasingly
        self.distribution_rgyr = dict(sorted(self.distribution_rgyr.items(), key=lambda item: item[0], reverse=True))
        self.concentration = np.mean(self.timeline_c)
    
    def write_file_header(self, path_to_directory, number_of_frames) -> None:
        """
        Initializes the output file with a header.

        Parameters:
        -----------
            path_to_directory (str): The directory where the output file will be saved.
            number_of_frames (int): The number of frames used in averaging.
        """
        filename = f"gyration_radius_distribution-{self.info}.dat"
        if not os.path.exists(path_to_directory):
            os.makedirs(path_to_directory)

        self.filepath_rgyr = os.path.join(path_to_directory, filename)
        
        with open(self.filepath_rgyr, 'w', encoding='utf-8') as output:
            output.write(f"# Gyration radius distribution \u279c {number_of_frames} frames analysed.\n")
            output.write("# Size \u279c Gyration radius +/- Error\n")
        output.close()

    def append_results_to_file(self):
        """
        Appends the result to the output file.
        """
        
        # Get largest list of radii
        max_length = 0
        for size, radii in self.distribution_rgyr.items():
            if len(radii) > max_length:
                max_length = len(radii)
        
        # Fill the lists with zeros
        for size, radii in self.distribution_rgyr.items():
            if len(radii) < max_length:
                self.distribution_rgyr[size] += [0.0] * (max_length - len(radii))
        
        with open(self.filepath_rgyr, 'a', encoding='utf-8') as output:
            output.write(f"# Concentration \u279c {self.concentration}\n")
            for size, radii in self.distribution_rgyr.items():
                if size == 0:
                    continue
                if len(radii) > 1 :
                    output.write(f"{size:7d} \u279c {np.mean([r for r in radii if r != 0]):<5.5f} +/- {np.std(radii)/np.sqrt(len(radii)):5.5f} # [")
                    counter = 0
                    for r in radii:
                        if counter <= 6:
                            if r != 0.:
                                output.write(f"{r:5.5f}, ")
                                counter += 1
                    if counter > 6:
                        output.write("... ]\n")
                    else:
                        output.write("]\n")
                elif len(radii) == 1:
                    output.write(f"{size:7d} \u279c {radii[0]:<5.5f} +/- {0.00000} # {radii}\n")
                else:
                    output.write(f"{size:7d} \u279c {0.00000:<5.5f} +/- {0.00000} # {radii}\n")
        output.close()

class CorrelationLength(Result):
    """
    Represents a result for the correlation length.
    
    Attributes:
    -----------
        corre_length (float): The computed correlation length.
        filepath (str): The path to the output file.
    """
    def __init__(self, name, info, init_frame) -> None:
        """
        Initializes a CorrelationLengthDistributionResult object.
        
        Parameters:
        -----------
            name (str): The property name.
            info (str): Additional information about the result.
            init_frame (int): The initial frame.
        """
        super().__init__(name, info, init_frame)
        self.corre_length = 0.
        self.timeline = {}
        self.filepath = None
    
    def add_to_timeline(self, frame, value, concentration) -> None:
        """
        Appends a data point to the timeline.

        Parameters:
        -----------
            frame (int): The frame number.
            value: The data point to add.
            concentration: The concentration associated with the data point.
        """
        self.timeline[frame] = value
        self.timeline_c.append(concentration)
        
    def calculate_correlation_length(self) -> None:
        r"""
        Calculates the gyration radius distribution based on the timeline data.
        
        Mathematically, the correlation length is defined as:
        $\(\xi = \sqrt{\frac{2 \sum_{i=1}^{N} r_{i}^{2} s_{i}^{2}}{\sum_{i=1}^{N} s_{i}^{2}}}$\) where:
            - \(r_{i}\) is the radius of the cluster \(i\),
            - \(s_{i}\) is the size of the cluster \(i\),
            - \(N\) is the number of clusters.
        
        """

        lists_r_s = {}
        lists_s_ns = {}
        correlation_lengths = {}
        numerators = {}
        
        n = 0 # number of all clusters
               
        for frame, sizes in self.timeline.items():
            for size, radii in sizes.items():
                ns = len(radii)
                s = size
                if frame not in lists_r_s:
                    lists_r_s[frame] = [[], []] # r, s, n
                    lists_s_ns[frame] = [[], []] # s, n
                    numerators[frame] = 0
                lists_s_ns[frame][0].append(s)
                lists_s_ns[frame][1].append(ns)
                for r in radii:
                    lists_r_s[frame][0].append(r)
                    lists_r_s[frame][1].append(s)
                n += ns
                
            den = np.sum([s**2 * n for s, n in zip(lists_s_ns[frame][0], lists_s_ns[frame][1])])
            if den == 0: den = 1 # avoid zero division
            
            numerators[frame] = np.sum([2 * r**2 * s**2 for r, s in zip(lists_r_s[frame][0], lists_r_s[frame][1])])
            
            correlation_lengths[frame] = np.sqrt(numerators[frame] / den)
            
        if n == 0: n = 1 # avoid zero division
               
        self.corre_length = np.mean(list(correlation_lengths.values()))
        self.error = np.std(list(correlation_lengths.values())) / np.sqrt(n)
        self.concentration = np.mean(self.timeline_c)

    def write_file_header(self, overwrite, path_to_directory, number_of_frames) -> None:
        """
        Initializes the output file with a header.

        Parameters:
        -----------
            overwrite (bool): Whether to overwrite the existing file.
            path_to_directory (str): The directory where the output file will be saved.
            number_of_frames (int): The number of frames used in averaging.
        """
        filename = f"correlation_length.dat"
        if not os.path.exists(path_to_directory):
            os.makedirs(path_to_directory)

        self.filepath = os.path.join(path_to_directory, filename)
        
        if not overwrite and os.path.exists(self.filepath):
            with open(self.filepath, 'a', encoding='utf-8') as output:
                output.write(f"# Correlation length \u279c {number_of_frames} frames averaged.\n")
                output.write("# Concentration \u279c Correlation length +/- Error\n")
            output.close()
        else:
            with open(self.filepath, 'w', encoding='utf-8') as output:
                output.write(f"# Correlation length \u279c {number_of_frames} frames averaged.\n")
                output.write("# Concentration \u279c Correlation length +/- Error\n")
            output.close()
    
    def append_results_to_file(self):
        """
        Appends the result to the output file.
        """
        
        with open(self.filepath, 'a', encoding='utf-8') as output:
            output.write(f"{self.concentration:10.6f} \u279c {self.corre_length:10.5f} +/- {self.error:<10.5f} # {self.info}\n")
        output.close()
        
        make_lines_unique(self.filepath)
    
class OrderParameter(Result):
    """
    Represents a result for the order parameter.
    
    Attributes:
    -----------
        order_parameter (float): The computed order parameter.
        filepath (str): The path to the output file.
    """
    def __init__(self, name, info, init_frame) -> None:
        """
        Initializes an OrderParameterResult object.
        
        Parameters:
        -----------
            name (str): The property name.
            info (str): Additional information about the result.
            init_frame (int): The initial frame.
        """
        super().__init__(name, info, init_frame)
        self.order_parameter = 0.
        self.timeline = {}
        self.filepath = None
        
    def add_to_timeline(self, frame: int, value: list, concentration: float) -> None:
        """
        Appends a data point to the timeline.

        Parameters:
        -----------
            frame (int): The frame number.
            value (list): The data point to add.
            concentration (float): The concentration associated with the data point.
        """
        self.timeline[frame] = value
        self.timeline_c.append(concentration)
    
    def calculate_order_parameter(self) -> None:
        """
        Calculates the order parameter based on the timeline data.
        """
        order_parameters = []
        for frame, order_parameter in self.timeline.items():
            order_parameters.append(order_parameter)
        
        self.order_parameter = np.mean(order_parameters, axis=0)
        self.error = np.std(order_parameters, axis=0) / np.sqrt(len(order_parameters))
        self.concentration = np.mean(self.timeline_c)
    
    def write_file_header(self, overwrite, path_to_directory, number_of_frames) -> None:
        """
        Initializes the output file with a header.

        Parameters:
        -----------
            overwrite (bool): Whether to overwrite the existing file.
            path_to_directory (str): The directory where the output file will be saved.
            number_of_frames (int): The number of frames used in averaging.
        """
        filename = f"order_parameter.dat"
        if not os.path.exists(path_to_directory):
            os.makedirs(path_to_directory)

        self.filepath = os.path.join(path_to_directory, filename)
        
        if not overwrite and os.path.exists(self.filepath):
            with open(self.filepath, 'a', encoding='utf-8') as output:
                output.write(f"# Order parameter \u279c {number_of_frames} frames averaged.\n")
                output.write("# Concentration \u279c Order parameter +/- Error\n")
            output.close()
        else:
            with open(self.filepath, 'w', encoding='utf-8') as output:
                output.write(f"# Order parameter \u279c {number_of_frames} frames averaged.\n")
                output.write("# Concentration \u279c Order parameter +/- Error\n")
            output.close()
    
    def append_results_to_file(self):
        """
        Appends the result to the output file.
        """
        
        with open(self.filepath, 'a') as output:
            for i in range(len(self.order_parameter)):
                output.write(f"{self.concentration:10.6f} \u279c {self.order_parameter[i]:10.6f} +/- {self.error[i]:<10.6f} # {self.info} - {i+1}D\n")
        output.close()
        
        make_lines_unique(self.filepath)

class PercolationProbability(Result):
    """
    Represents a result for the percolation probability.
    
    Attributes:
    -----------
        percolation_probability (float): The computed percolation probability.
        filepath (str): The path to the output file.
    """
    def __init__(self, name, info, init_frame) -> None:
        """
        Initializes a PercolationProbabilityResult object.
        
        Parameters:
        -----------
            name (str): The property name.
            info (str): Additional information about the result.
            init_frame (int): The initial frame.
        """
        super().__init__(name, info, init_frame)
        self.percolation_probability = 0.
        self.timeline = {}
        self.filepath = None
        
    def add_to_timeline(self, frame: int, value:list, concentration: float) -> None:
        """
        Appends a data point to the timeline.

        Parameters:
        -----------
            frame (int): The frame number.
            value (list): The data point to add.
            concentration (float): The concentration associated with the data point.
        """
        self.timeline[frame] = value
        self.timeline_c.append(concentration)
        
    def calculate_percolation_probability(self) -> None:
        """
        Calculates the percolation probability based on the timeline data.
        """
        percolation_probabilities = []
        for frame, percolation_probability in self.timeline.items():
            percolation_probabilities.append(percolation_probability)
        
        self.percolation_probability = np.mean(percolation_probabilities, axis=0)
        self.error = np.std(percolation_probabilities, axis=0) / np.sqrt(len(percolation_probabilities))
        self.concentration = np.mean(self.timeline_c)
        
    def write_file_header(self, overwrite, path_to_directory, number_of_frames) -> None:
        """
        Initializes the output file with a header.

        Parameters:
        -----------
            overwrite (bool): Whether to overwrite the existing file.
            path_to_directory (str): The directory where the output file will be saved.
            number_of_frames (int): The number of frames used in averaging.
        """
        filename = f"percolation_probability.dat"
        if not os.path.exists(path_to_directory):
            os.makedirs(path_to_directory)

        self.filepath = os.path.join(path_to_directory, filename)
        
        if not overwrite and os.path.exists(self.filepath):
            with open(self.filepath, 'a', encoding='utf-8') as output:
                output.write(f"# Percolation probability \u279c {number_of_frames} frames averaged.\n")
                output.write("# Concentration \u279c Percolation probability +/- Error\n")
            output.close()
        else:
            with open(self.filepath, 'w', encoding='utf-8') as output:
                output.write(f"# Percolation probability \u279c {number_of_frames} frames averaged.\n")
                output.write("# Concentration \u279c Percolation probability +/- Error\n")
            output.close()
        
    def append_results_to_file(self):
        """
        Appends the result to the output file.
        """
        
        with open(self.filepath, 'a') as output:
            for i in range(len(self.percolation_probability)):
                output.write(f"{self.concentration:10.6f} \u279c {self.percolation_probability[i]:10.6f} +/- {self.error[i]:<10.6f} # {self.info} - {i+1}D\n")
        output.close()
        
        make_lines_unique(self.filepath)