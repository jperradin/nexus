<div align="center">

  # NEXUS-CAT
  ##### Cluster Analysis Toolkit
  [![License](https://img.shields.io/badge/license-MIT-blue.svg)](https://opensource.org/licenses/MIT)
  [![PyPI version](https://badge.fury.io/py/nexus-cat.svg)](https://badge.fury.io/py/nexus-cat)
  [![Documentation Status](https://readthedocs.org/projects/nexus-cat/badge/?version=latest)](https://nexus-cat.readthedocs.io/en/latest/)

  <img alt="NEXUS-CAT" width=400 src="./assets/Logo_Nexus-CAT_RVB_1.png" />
</div>

## ⇁ TOC
- [Description and features](#-description-and-features)
- [Installation](#-installation)
- [Getting started](#-getting-started)
- [Documentation](#-documentation)
- [Contributing](#-contributing)
- [License](#-license)

## ⇁ Description and features

`nexus-cat` is a package designed to find clusters of connected polyhedra in an atomistic simulation trajectory. It provides functionality to analyze cluster properties according to the percolation theory:
- *Note: Here the notion of size refers to the number of polyhedra in a cluster, not the physical size of the cluster, ie its radius nor its volume.*
- **Average cluster size** $\langle s \rangle$: $\langle s(p) \rangle = \sum_s \frac{s^2n_s(p)}{\sum_s s n_s(p)}$
  - with $n_s$ the number of clusters of size $s$ (ie number of polyhedra in the cluster).
  - 1 sized clusters and percolating clusters are not taken into account in the calculation.
- **Biggest cluster size** $s_{max}$: largest cluster size in the system no matter the percolation threshold.
- **Spanning cluster size** $s_{\infty}$ : largest cluster size in the system excluding the percolating cluster.
- **Gyration radius** $R_g$ : $R_s² = \frac{1}{2s^2}\sum_{i,j}|\overrightarrow{r_i}-\overrightarrow{r_j}|^2$
  - with $r_i$ the **unwrapped** coordinates of the atom $_i$ in the cluster of size $s$. 
  - 1 sized clusters and percolating clusters are not taken into account in the calculation.
- **Correlation length** $\xi$ : $\xi^2 = \frac{\sum_s 2R_s²s²n_s(p)}{\sum_ss²n_s(p)}$
  - with $n_s$ the number, $R_s$ the average gyration radius of clusters of size $s$ (ie number of polyhedra in the cluster).
  - 1 sized clusters and percolating clusters are not taken into account in the calculation.
- **Percolation probability** $\Pi$ :
```math
\Pi = \begin{cases}
0 & \text{if } R_g < L_{box} \\
1 & \text{if } R_g \geq L_{box} 
\end{cases}
```
  - with $L_{box}$ is the length of the simulation box.
  - Note: The percolation probability is calculated for each direction of the simulation box, a cluster can percolate in 1D, 2D or 3D. 

- **Order parameter $P_∞$** : 
```math
P_∞ = \begin{cases}0 & \text{if } \Pi = 0 \\\frac{s_{max}}{N} & \text{if } \Pi = 1 
\end{cases}
```
 
  - with $s_{max}$ the number of polyhedra in the biggest cluster, $N$ the total number of **connected** polyhedra in the system (1 sized clusters excluded).
  - Note : the order parameter is calculated with $\Pi$ in 1D. 

## ⇁ Installation

To install `nexus`, first clone this repository as you please, for example with SSH:

```bash
git clone git@github.com:jperradin/nexus.git
```
Then you can use pip, it will install dependencies and the main package in your Python environment:

```bash
pip install nexus-cat
```

## ⇁ Getting started

As a first example you can follow the steps of the [Getting started](https://nexus-cat.readthedocs.io/en/latest/getting_started.html) section of the documentation.

## ⇁ Documentation

The documentation is available [here](https://nexus-cat.readthedocs.io/en/latest/)

## ⇁ Contributing

Contributions to `Nexus-CAT` are welcome! You can contribute by submitting bug reports, feature requests, new extension requests, or pull requests through GitHub.

## ⇁ License

This project is licensed under the [MIT License](https://opensource.org/licenses/MIT).

