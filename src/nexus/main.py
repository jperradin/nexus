from tqdm import tqdm
import numpy as np
import os
from datetime import datetime
import psutil
import uuid
import time

from .config.settings import Settings
from .io.reader.reader_factory import ReaderFactory
from .core.system import System
from .analysis.strategy_factory import StrategyFactory
from .analysis.analyzer_factory import AnalyzerFactory
from .io.writer.writer_factory import WriterFactory
from .utils import *
from .version import __version__

def main(settings: Settings):
    """
    Main function to test the package.
    """

    perf = performance.Performance(
        id=str(uuid.uuid4()),
        name=f"{settings.project_name}-run_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
    )
    start_time = time.time()
    process = psutil.Process(os.getpid())
    init_memory = process.memory_info().rss / (1024 * 1024)  # Convert to MB

    # Print title and settings
    if settings.verbose:
        aesthetics.print_title(__version__)
    if settings.verbose:
        print(settings)

    # Create export directory
    settings.export_directory = os.path.join(settings.export_directory, settings.project_name)
    if not os.path.exists(settings.export_directory):
        os.makedirs(settings.export_directory)

    # Save logs
    if settings.save_logs:
        writer = WriterFactory(settings).get_writer('LogsWriter')
        writer.write()

    # Initialize reader and system
    scan_start = time.time()
    reader = ReaderFactory(settings).get_reader()
    reader.set_verbose(settings.verbose)
    system = System(reader, settings)
    scan_end = time.time()
    
    # Track reader initialization performance
    perf.add_metric('scan_trajectory_time_ms', (scan_end - scan_start) * 1000)
    perf.record_history()

    # Get total number of frames
    if settings.range_of_frames[1] == -1:
        total = system.get_num_frames()
    elif settings.range_of_frames[0] == -1:
        settings.range_of_frames[0] = 1
    elif settings.range_of_frames[0] == settings.range_of_frames[1]:
        total = 1
    else:
        if settings.range_of_frames[1] == -1:
            total = system.get_num_frames() - settings.range_of_frames[0]
        else:
            total = settings.range_of_frames[1] - settings.range_of_frames[0]

    # Initialize analyzers
    analyzers = []
    for analyzer in settings.analysis.get_analyzers():
        analyzers.append(AnalyzerFactory(settings).get_analyzer(analyzer))

    # Initialize progress bar
    progress_bar_kwargs = {
            "disable": not settings.verbose,
            "leave": True,
            "ncols": os.get_terminal_size().columns,
            "colour": "red"
        }

    progress_bar = tqdm(enumerate(system.iter_frames()), desc=f"Processing frames {settings.range_of_frames}...", unit="frame", total=total, **progress_bar_kwargs)

    # Track per-frame metrics
    frame_times = []
    neighbor_times = []
    cluster_times = []
    analysis_times = []
    number_nodes = []

    # Read and process frames
    for i, frame in progress_bar:

        frame_start = time.time()

        if settings.lattice.apply_custom_lattice:
            frame.set_lattice(settings.lattice.custom_lattice)
        
        # Initialize nodes
        frame.initialize_nodes()
        # Find neighbors
        neighbor_start = time.time()
        strategy = StrategyFactory(frame, settings).get_strategy(settings)
        strategy.find_neighbors()
        neighbor_end = time.time()
        neighbor_times.append((neighbor_end - neighbor_start) * 1000)
        
        # Find clusters
        cluster_start = time.time()
        connectivities = strategy.get_connectivities()
        clusters = strategy.build_clusters()
        frame.set_clusters(clusters)
        frame.set_connectivities(connectivities)
        cluster_end = time.time()
        cluster_times.append((cluster_end - cluster_start) * 1000)

        # Analyze clusters
        analysis_start = time.time()
        for analyzer in analyzers:
            analyzer.analyze(frame, connectivities)
        analysis_end = time.time()
        analysis_times.append((analysis_end - analysis_start) * 1000)

        # Print clusters
        if settings.clustering.with_printed_unwrapped_clusters:
            writer = WriterFactory(settings).get_writer('ClustersWriter')
            writer.set_clusters(frame.get_clusters())
            writer.write()

        frame_end = time.time()
        frame_times.append((frame_end - frame_start) * 1000)

        number_nodes.append(len(frame))

        del frame

        # Record per-frame performance
        if i % 10 == 0 or i == total - 1:  # Record every 10 frames or the last frame
            current_memory = process.memory_info().rss / (1024 * 1024)
            cpu_percent = process.cpu_percent()
            
            perf.execution_time_ms = frame_times[-1]
            perf.memory_usage_mb = current_memory
            perf.cpu_usage_percent = cpu_percent
            perf.add_metric("frame_number", i)
            perf.add_metric("neighbor_finding_time_ms", neighbor_times[-1])
            perf.add_metric("cluster_finding_time_ms", cluster_times[-1])
            perf.add_metric("analysis_time_ms", analysis_times[-1])
            perf.add_metric("number_nodes", number_nodes[-1])
            perf.record_history()

    # Print results
    for analyzer in analyzers:
        analyzer.print_to_file()       

    # Record overall performance metrics
    end_time = time.time()
    total_time_ms = (end_time - start_time) * 1000
    final_memory = process.memory_info().rss / (1024 * 1024)
    memory_increase = final_memory - init_memory
    
    perf.execution_time_ms = total_time_ms
    perf.memory_usage_mb = final_memory
    perf.cpu_usage_percent = process.cpu_percent()
    perf.add_metric("total_frames_processed", total)
    perf.add_metric("memory_increase_mb", memory_increase)
    perf.add_metric("avg_frame_time_ms", sum(frame_times) / len(frame_times) if frame_times else 0)
    perf.add_metric("avg_neighbor_time_ms", sum(neighbor_times) / len(neighbor_times) if neighbor_times else 0)
    perf.add_metric("avg_cluster_time_ms", sum(cluster_times) / len(cluster_times) if cluster_times else 0)
    perf.add_metric("avg_analysis_time_ms", sum(analysis_times) / len(analysis_times) if analysis_times else 0)
    perf.add_metric("avg_number_nodes", sum(number_nodes) / len(number_nodes) if number_nodes else 0)
    perf.record_history()

    # Save performance metrics
    if settings.save_performance:
        writer = WriterFactory(settings).get_writer('PerformanceWriter')
        writer.write(perf)
           
        
        
    