"""
Main module for the Nexus project.

This module serves as the entry point for the Nexus application.
It initializes the necessary components and orchestrates the workflow.

Imports:
    - Standard libraries
    - Third-party libraries
    - Internal modules

Functions:
    - main: The main function that initializes and runs the application.
"""

# Standard library imports
import os
import importlib

# Third-party imports
import numpy as np
from tqdm import tqdm

# Internal imports
from . import io
from . import core
from .utils.generate_color_gradient import generate_color_gradient

def main(settings):
    """
    Analyzes a trajectory based on the provided settings.

    Parameters:
    -----------
        settings (Settings): An object containing settings for the analysis.

    This function performs an analysis on a trajectory based on the provided settings.
    It creates necessary directories, counts configurations, imports required modules,
    reads lattice properties, sets up cutoffs, and calculates percolation properties
    for each frame within the specified range. The analysis includes finding clusters
    based on connectivity criteria, calculating cluster sizes, and writing results
    to output files.

    Returns:
    --------
        None
    """
    # Build the output directory
    new_directory = os.path.join(settings.export_directory.get_value(), settings.project_name.get_value())
    
    settings._output_directory = new_directory
    
    # Create the output directory if it does not exist
    if not os.path.exists(settings._output_directory):
        os.makedirs(settings._output_directory)
    
    input_file = settings.path_to_xyz_file.get_value()
    
    # Count the number of configurations in the trajectory
    n_config = io.count_configurations(input_file)
    n_atoms  = settings.number_of_atoms.get_value()
    n_header = settings.header.get_value()
    settings.number_of_frames.set_value(n_config)
    
    settings.print_settings()
    
    # Import the extension
    module = importlib.import_module(f"nexus.extensions.{settings.extension.get_value()}")
    
    # Create the Box object and append lattice for each frame
    box = core.Box()
    io.read_lattices_properties(box, input_file)
    
    # Create the Cutoff object
    cutoffs = core.Cutoff(settings.cutoffs.get_value())
    
    # Setting up the results objects for the percolation properties
    connectivities = module.get_connectivity(settings.cluster_settings.get_value())
    
    # Check if there any extra clustering method in the extension
    extra_function_to_call = module.EXTRA_CLUSTERING_METHODS
        
    # Setting the for loop with user settings
    if settings.range_of_frames.get_value() is not None:
        start = settings.range_of_frames.get_value()[0]
        end   = settings.range_of_frames.get_value()[1]
    else:
        start = 0
        end   = n_config
        
    if end-start == 0:
        raise ValueError(f'\t\tERROR: Range of frames selected is wrong \u279c {settings.range_of_frames.get_value()}.')
    else:
        settings.frames_to_analyse.set_value(end-start)

    color_gradient = generate_color_gradient(end-start)
    if settings.quiet.get_value():
        progress_bar = range(start, end)
    else:
        progress_bar = tqdm(range(start, end), desc="Analysing trajectory ... ", unit='frame', leave=False, colour="YELLOW")
    
    # Create the Result objects
    results_average_cluster_size = {}
    results_biggest_cluster_size = {}
    results_spanning_cluster_size = {}
    results_cluster_size_distribution = {}
    results_gyration_radius_distribution = {}
    results_correlation_length = {}
    results_order_parameter = {}
    results_percolation_probability = {}
    
    # Keep results of previous runs if overwrite_results is False
    overwrite_results = settings.overwrite_results.get_value()
    
    for c in connectivities:
        results_average_cluster_size[c] = io.AverageClusterSize('average_cluster_size', c, start)
        results_biggest_cluster_size[c] = io.BiggestClusterSize('biggest_cluster_size', c, start)
        results_spanning_cluster_size[c] = io.SpanningClusterSize('spanning_cluster_size', c, start)
        results_cluster_size_distribution[c] = io.ClusterSizeDistribution('cluster_size_distribution', c, start)
        results_gyration_radius_distribution[c] = io.GyrationRadiusDistribution('gyration_radius_distribution', c, start)
        results_correlation_length[c] = io.CorrelationLength('correlation_length', c, start)
        results_order_parameter[c] = io.OrderParameter('order_parameter', c, start)
        results_percolation_probability[c] = io.PercolationProbability('percolation_probability', c, start)
        results_average_cluster_size[c].write_file_header(overwrite_results, settings._output_directory, end-start)
        results_biggest_cluster_size[c].write_file_header(overwrite_results, settings._output_directory, end-start)
        results_spanning_cluster_size[c].write_file_header(overwrite_results, settings._output_directory, end-start)
        results_cluster_size_distribution[c].write_file_header(settings._output_directory, end-start)
        results_gyration_radius_distribution[c].write_file_header(settings._output_directory, end-start)
        results_correlation_length[c].write_file_header(overwrite_results, settings._output_directory, end-start)
        results_order_parameter[c].write_file_header(overwrite_results, settings._output_directory, end-start)
        results_percolation_probability[c].write_file_header(overwrite_results, settings._output_directory, end-start)
    
    # Adding the results for extra clustering method
    if extra_function_to_call is not None:
        if extra_function_to_call:
            for c in module.get_extra_connectivity(settings.cluster_settings.get_value()):
                results_average_cluster_size[c] = io.AverageClusterSize('average_cluster_size', c, start)
                results_biggest_cluster_size[c] = io.BiggestClusterSize('biggest_cluster_size', c, start)
                results_spanning_cluster_size[c] = io.SpanningClusterSize('spanning_cluster_size', c, start)
                results_cluster_size_distribution[c] = io.ClusterSizeDistribution('cluster_size_distribution', c, start)
                results_gyration_radius_distribution[c] = io.GyrationRadiusDistribution('gyration_radius_distribution', c, start)
                results_correlation_length[c] = io.CorrelationLength('correlation_length', c, start)
                results_order_parameter[c] = io.OrderParameter('order_parameter', c, start)
                results_percolation_probability[c] = io.PercolationProbability('percolation_probability', c, start)
                results_average_cluster_size[c].write_file_header(overwrite_results, settings._output_directory, end-start)
                results_biggest_cluster_size[c].write_file_header(overwrite_results, settings._output_directory, end-start)
                results_spanning_cluster_size[c].write_file_header(overwrite_results, settings._output_directory, end-start)
                results_cluster_size_distribution[c].write_file_header(settings._output_directory, end-start)
                results_gyration_radius_distribution[c].write_file_header(settings._output_directory, end-start)
                results_correlation_length[c].write_file_header(overwrite_results, settings._output_directory, end-start)
                results_order_parameter[c].write_file_header(overwrite_results, settings._output_directory, end-start)
                results_percolation_probability[c].write_file_header(overwrite_results, settings._output_directory, end-start)
    
    # Loop over the frames in trajectory
    for i in progress_bar:
        # Update progress bar
        if not settings.quiet.get_value():
            progress_bar.set_description(f"Analysing frame n°{i}")
            progress_bar.colour = "#%02x%02x%02x" % color_gradient[i-start]
        
        # Create the System object at the current frame
        system = io.read_and_create_system(input_file, i, n_atoms+n_header, settings, cutoffs, start, end)
        system.frame = i
        
        # Set the Box object to the System object
        system.box = box
        settings.lbox.set_value(system.box.get_box_dimensions(i))
        
        # Calculate the nearest neighbours of all atoms in the system.
        system.calculate_neighbours()
        
        # Calculate the concentrations (ie, the number of sites in the lattice)
        dict_units = system.calculate_concentrations(settings.extension.get_value())
        
        for c in connectivities:
            system.find_clusters(c)
            system.set_concentrations(c)
            
            concentration = system.get_concentration(c)
        
            list_sizes = system.get_filtered_cluster_sizes(c)
            list_all_sizes = system.get_all_cluster_sizes(c)
            
            results_average_cluster_size[c].add_to_timeline({i : list_sizes}, concentration)
            results_biggest_cluster_size[c].add_to_timeline({i : list_all_sizes}, concentration)
            results_spanning_cluster_size[c].add_to_timeline({i : list_sizes}, concentration)
            
            dict_sizes = system.get_cluster_sizes_distribution(c)
            dict_rgyr = system.get_gyration_radius_distribution(c, list_sizes)
            
            results_cluster_size_distribution[c].add_to_timeline(i, dict_sizes, concentration)
            results_gyration_radius_distribution[c].add_to_timeline(i, dict_rgyr, concentration)
            results_correlation_length[c].add_to_timeline(i, dict_rgyr, concentration)
            
            order_parameter = system.calculate_order_parameter(c)
            percolation_probability = system.calculate_percolation_probability(c)
            
            results_order_parameter[c].add_to_timeline(i, order_parameter, concentration)
            results_percolation_probability[c].add_to_timeline(i, percolation_probability, concentration)
        
        if extra_function_to_call is not None:
            if extra_function_to_call:
                for c in module.get_extra_connectivity(settings.cluster_settings.get_value()):
                    system.find_extra_clusters()
                    
                    system.set_concentrations(c)
            
                    concentration = system.get_concentration(c)
                    list_sizes = system.get_filtered_cluster_sizes(c)
                    list_all_sizes = system.get_all_cluster_sizes(c)
                    
                    results_average_cluster_size[c].add_to_timeline({i : list_sizes}, concentration)
                    results_biggest_cluster_size[c].add_to_timeline({i : list_all_sizes}, concentration)
                    results_spanning_cluster_size[c].add_to_timeline({i : list_sizes}, concentration)
                    
                    dict_sizes = system.get_cluster_sizes_distribution(c)
                    dict_rgyr = system.get_gyration_radius_distribution(c, list_sizes)
                    
                    results_cluster_size_distribution[c].add_to_timeline(i, dict_sizes, concentration)
                    results_gyration_radius_distribution[c].add_to_timeline(i, dict_rgyr, concentration)
                    results_correlation_length[c].add_to_timeline(i, dict_rgyr, concentration)
                    
                    order_parameter = system.calculate_order_parameter(c)
                    percolation_probability = system.calculate_percolation_probability(c)
                    
                    results_order_parameter[c].add_to_timeline(i, order_parameter, concentration)
                    results_percolation_probability[c].add_to_timeline(i, percolation_probability, concentration)
    
    for c in connectivities:
        results_average_cluster_size[c].calculate_average_cluster_size()
        results_average_cluster_size[c].append_results_to_file()
        results_biggest_cluster_size[c].calculate_biggest_cluster_size()
        results_biggest_cluster_size[c].append_results_to_file()
        results_spanning_cluster_size[c].calculate_spanning_cluster_size()
        results_spanning_cluster_size[c].append_results_to_file()
        results_cluster_size_distribution[c].calculate_cluster_size_distribution()
        results_cluster_size_distribution[c].append_results_to_file()
        results_gyration_radius_distribution[c].calculate_gyration_radius_distribution()
        results_gyration_radius_distribution[c].append_results_to_file()
        results_correlation_length[c].calculate_correlation_length()
        results_correlation_length[c].append_results_to_file()
        results_order_parameter[c].calculate_order_parameter()
        results_order_parameter[c].append_results_to_file()
        results_percolation_probability[c].calculate_percolation_probability()
        results_percolation_probability[c].append_results_to_file()
    
    if extra_function_to_call is not None:
        if extra_function_to_call:
            for c in module.get_extra_connectivity(settings.cluster_settings.get_value()):
                results_average_cluster_size[c].calculate_average_cluster_size()
                results_average_cluster_size[c].append_results_to_file()
                results_biggest_cluster_size[c].calculate_biggest_cluster_size()
                results_biggest_cluster_size[c].append_results_to_file()
                results_spanning_cluster_size[c].calculate_spanning_cluster_size()
                results_spanning_cluster_size[c].append_results_to_file()
                results_cluster_size_distribution[c].calculate_cluster_size_distribution()
                results_cluster_size_distribution[c].append_results_to_file()
                results_gyration_radius_distribution[c].calculate_gyration_radius_distribution()
                results_gyration_radius_distribution[c].append_results_to_file()
                results_correlation_length[c].calculate_correlation_length()
                results_correlation_length[c].append_results_to_file()
                results_order_parameter[c].calculate_order_parameter()
                results_order_parameter[c].append_results_to_file()
                results_percolation_probability[c].calculate_percolation_probability()
                results_percolation_probability[c].append_results_to_file()
    
    if settings.print_clusters_positions.get_value():
        # Generating a list.txt file with the paths of the all-in-one unwrapped clusters files.
        unwrapped_clusters_dir_path = os.path.join(settings._output_directory, 'unwrapped_clusters')
        io.write_list_of_files(unwrapped_clusters_dir_path)
