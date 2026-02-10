# Features

Nexus-CAT provides a wide range of features for cluster analysis in atomistic simulations. Below is an overview of the key features:

### Cluster Analysis

- **Multiple Clustering Strategies**: Nexus-CAT supports several strategies for defining clusters:
    - **Distance Strategy**: Connects nodes based on a simple distance cutoff.
    - **Bonding Strategy**: Identifies clusters based on a three-node bonding pattern (e.g., Si-O-Si).
    - **Coordination Strategy**: Extends the bonding pattern by adding constraints on the coordination number of the nodes.
    - **Shared Strategy**: A more advanced strategy that connects nodes based on a minimum number of shared neighbors, useful for analyzing polyhedral linkages.

### Percolation Properties

- **Average Cluster Size**: Calculate the average size of clusters in the system.
- **Largest Cluster Size**: Identify the largest cluster in the system.
- **Spanning Cluster Size**: Exclude percolating clusters and calculate the largest remaining cluster.
- **Gyration Radius**: Measure the spatial extent of clusters.
- **Correlation Length**: Calculate the correlation length of clusters.
- **Percolation Probability**: Determine the probability of clusters percolating in 1D, 2D, or 3D.
- **Cluster Size Distribution**: Computes the distribution of cluster sizes, n(s), for each connectivity type.
- **Concentration**: Computes the concentration of clusters for each connectivity type.

### Extensibility

- Customize cluster connectivity criterion and polyhedra.
- The use of factories (`StrategyFactory`, `AnalyzerFactory`, `ReaderFactory`, `WriterFactory`) allows for easy extension of the package with new algorithms and file formats.

### Visualization and Output

- **XYZ Cluster Exporter**: The `ClustersWriter` generates XYZ files containing the unwrapped atomic coordinates of the clusters, which can be used for visualization in other software.
- **Data Summaries**: The `MultipleFilesSummaryWriter` can be used to summarize data from multiple analysis files, making it easier to compare results across different simulations.
- **Performance Metrics**: The `PerformanceWriter` saves detailed performance metrics, which can be useful for benchmarking and optimization.

For a complete list of features, refer to the [API Reference](api_reference.md).