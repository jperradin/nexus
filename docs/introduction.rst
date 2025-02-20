.. _introduction:

Introduction
============

Welcome to the Nexus-Cat Python package documentation. This section provides an overview of the project, its goals, and its features.


Project Goals
-------------
The Nexus-cat package aims to provide an assets of tools for the analysis of clusters within atomistic simulations trajectories.
The package is designed to be user-friendly, efficient, and flexible, allowing easy analysis and visualize their simulation data.
Extensions can be easily added to the package to be able to analyze different systems and cluster connectivities.

Percolation Analysis
--------------------

The package includes calculations of percolation properties of clusters.

Here is a non-exhaustive list of the percolation properties that can be calculated:

Average Cluster Size
====================

**Average cluster size** :math:`\langle s \rangle`:

.. math::

    \langle s(p) \rangle = \sum_s \frac{s^2 n_s(p)}{\sum_s s n_s(p)}

- with :math:`n_s` the number of clusters of size :math:`s` (i.e., number of polyhedra in the cluster).
- 1-sized clusters and percolating clusters are not taken into account in the calculation.

Biggest Cluster Size
====================

**Biggest cluster size** :math:`s_{\max}`:  
Largest cluster size in the system, regardless of the percolation threshold.

Spanning Cluster Size
=====================

**Spanning cluster size** :math:`s_{\infty}`:  
Largest cluster size in the system excluding the percolating cluster.

Gyration Radius
===============

**Gyration radius** :math:`R_g`:

.. math::

    R_s^2 = \frac{1}{2s^2} \sum_{i,j} |\overrightarrow{r_i} - \overrightarrow{r_j}|^2

- with :math:`r_i` the **unwrapped** coordinates of atom :math:`i` in the cluster of size :math:`s`.
- 1-sized clusters and percolating clusters are not taken into account in the calculation.

Correlation Length
==================

**Correlation length** :math:`\xi`:

.. math::

    \xi^2 = \frac{\sum_s 2 R_s^2 s^2 n_s(p)}{\sum_s s^2 n_s(p)}

- with :math:`n_s` the number and :math:`R_s` the average gyration radius of clusters of size :math:`s` (i.e., number of polyhedra in the cluster).
- 1-sized clusters and percolating clusters are not taken into account in the calculation.

Percolation Probability
=======================

**Percolation probability** :math:`\Pi`:

.. math::

    \Pi =
    \begin{cases}
        0 & \text{if } R_g < L_{\text{box}} \\
        1 & \text{if } R_g \geq L_{\text{box}}
    \end{cases}

- with :math:`L_{\text{box}}` the length of the simulation box.
- **Note**: The percolation probability is calculated for each direction of the simulation box; a cluster can percolate in 1D, 2D, or 3D.

Order Parameter
===============

**Order parameter** :math:`P_{\infty}`:

.. math::

    P_{\infty} =
    \begin{cases}
        0 & \text{if } \Pi = 0 \\
        \frac{s_{\max}}{N} & \text{if } \Pi = 1
    \end{cases}

- with :math:`s_{\max}` the number of polyhedra in the biggest cluster, and :math:`N` the total number of **connected** polyhedra in the system (excluding 1-sized clusters).
- **Note**: The order parameter is calculated with :math:`\Pi` in 1D.
