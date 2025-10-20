<div align="center">

# NEXUS-CAT
##### Cluster Analysis Toolkit
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![PyPI version](https://badge.fury.io/py/nexus-cat.svg)](https://badge.fury.io/py/nexus-cat)
[![Documentation Status](https://readthedocs.org/projects/nexus-cat/badge/?version=latest)](https://nexus-cat.readthedocs.io/en/latest/)

<img alt="NEXUS-CAT" width=400 src="./assets/Logo_Nexus-CAT_RVB_1.png" />
</div>

logo made by [Lisap](https://lisaperradinportfolio.framer.website/)

## ⇁ TOC
- [NEXUS-CAT](#nexus-cat)
        - [Cluster Analysis Toolkit](#cluster-analysis-toolkit)
  - [⇁ TOC](#⇁-toc)
  - [⇁ Who is this for?](#⇁-who-is-this-for)
  - [⇁ Description and features](#⇁-description-and-features)
  - [⇁ Installation](#⇁-installation)
  - [⇁ Getting started](#⇁-getting-started)
  - [⇁ Documentation](#⇁-documentation)
  - [⇁ Contributing](#⇁-contributing)
  - [⇁ License](#⇁-license)


## ⇁ Who is this for?

`nexus-cat` is designed for researchers, scientists, and students in computational materials science, condensed matter physics, and physical chemistry. It is particularly well-suited for those who:

* Work with atomistic simulation data from molecular dynamics (MD) or Monte Carlo (MC) simulations, especially in formats like XYZ and LAMMPS.
* Are studying phenomena related to **percolation theory**, such as the determination of percolation thresholds and critical exponents.
* Investigate the structure of **disordered materials**, such as glasses and amorphous solids. The package's advanced clustering strategies can analyze complex network structures, like the connectivity of silicate polyhedra in silica glass.
* Need to characterize the formation and properties of clusters, aggregates, or networks in their systems.

The package provides a flexible and powerful tool for moving beyond simple structural metrics and exploring the rich, multi-scale organization inherent in many-body systems.

Of course. Here is a more detailed feature description section that could be used to update the package's documentation.

## ⇁ Description and Features

`nexus-cat` offers a comprehensive suite of tools for network and cluster analysis, with a strong focus on concepts from percolation theory.

---

### **Clustering Strategies** 

At the core of `nexus-cat` are its flexible clustering strategies, which define how connections are made between nodes to form clusters. You can choose the strategy that best suits the physics of your system.

* **Distance Strategy**: This is the most straightforward approach, connecting any two nodes that are within a specified cutoff distance of each other. It's ideal for identifying simple aggregates or molecular groups.
* **Bonding Strategy**: This strategy identifies clusters based on a three-node bonding pattern, such as `Si-O-Si`. It connects two nodes if they share a common "bridging" atom.
* **Coordination Strategy**: Building on the bonding strategy, this approach adds constraints on the coordination number of the connected nodes. This allows for more specific structural identification, such as linking only 4-coordinated silicon atoms.
* **Shared Strategy**: This advanced strategy connects nodes based on a minimum number of shared neighbors. It is particularly powerful for distinguishing between different types of polyhedral linkages (e.g., corner-, edge-, or face-sharing) in complex, dense materials.

---

### **Percolation and Structural Analysis** 

Once clusters are identified, `nexus-cat` provides a range of analyzers to quantify their properties.

* **Average Cluster Size ($\langle S \rangle$)**: Calculates the weight-average cluster size, a key metric in percolation theory that diverges at the percolation threshold. The calculation is defined as $\langle S(p) \rangle = \sum_s \frac{s^2n_s(p)}{\sum_s s n_s(p)}$, where $n_s$ is the number of clusters of size $s$. To focus on the finite cluster distribution, percolating clusters are excluded from this calculation.
* **Largest Cluster Size ($S_{max}$)**: Identifies the size of the largest cluster in the system, regardless of whether it percolates or not.
* **Spanning Cluster Size ($S_{span}$)**: Determines the size of the largest *finite* (non-percolating) cluster.
* **Gyration Radius ($R_g$)**: Measures the spatial extent of a cluster. It is calculated based on the unwrapped coordinates of the atoms within the cluster to correctly handle periodic boundary conditions.
* **Correlation Length ($\xi$)**: Calculates the characteristic size of clusters using the second moment of the gyration radius distribution.
* **Percolation Probability ($\Pi$)**: Determines the probability of finding a cluster that spans the simulation box along any given dimension. A cluster is considered to be percolating if its span is greater than the corresponding lattice vector.
* **Order Parameter ($P_{\infty}$)**: Calculates the fraction of networking atoms that belong to a percolating cluster, a key metric for identifying the percolation threshold.

---

### **Efficient and Extensible Framework** 

* **I/O**: Includes efficient readers for **XYZ** and **LAMMPS** trajectory files, which scan and index frames for fast, on-demand parsing.
* **Configuration**: A **builder pattern** (`SettingsBuilder`) and **dataclasses** (`GeneralSettings`, `ClusteringSettings`, etc.) provide a clear, robust, and validated way to configure your analysis.
* **Extensibility**: The factory pattern for readers, writers, analyzers, and strategies makes it straightforward to extend the package with new file formats or custom analysis algorithms.
* **Performance Tracking**: The package can track and save detailed performance metrics, including execution time, memory usage, and CPU usage, helping you optimize your analysis workflows.

## ⇁ Installation

### Basic installation

To install `nexus-cat` as a package, you can use pip:

```bash
pip install nexus-cat
```

Note: the package does not auto upgrade itself, please run the following command to upgrade to the latest version:

```bash
pip install nexus-cat --upgrade
```

### Installation from the source code

If you want to install the package from the source code to implement your extensions for example, you can clone the repository:

```bash
git clone git@github.com:jperradin/nexus.git
```

Then install the package in development mode:

```bash
cd nexus
pip install -e .
```

## ⇁ Getting started

As a first example you can follow the steps of the [Getting started](https://nexus-cat.readthedocs.io/en/latest/getting_started.html) section of the documentation.

## ⇁ Documentation

The documentation is available [here](https://nexus-cat.readthedocs.io/en/latest/).

## ⇁ Notes

Portions of the documentation and this code were assisted by large language models or AI coding assistants and have been reviewed for accuracy and compliance.

## ⇁ License

This project is licensed under the [MIT License](https://opensource.org/licenses/MIT).
