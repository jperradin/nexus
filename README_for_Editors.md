# NEXUS-CAT — Cluster Analysis Toolkit
## Code Verification Guide for Editors

**Package name:** nexus-cat  
**Version:** 0.1.0  
**License:** MIT  
**Python compatibility:** 3.9, 3.10, 3.11, 3.12  
**Documentation:** https://nexus-cat.readthedocs.io/en/latest/  
**Source code:** https://github.com/jperradin/nexus

---

## Table of Contents

1. [Overview](#1-overview)
2. [System Requirements](#2-system-requirements)
3. [Installation](#3-installation)
4. [Repository Structure](#4-repository-structure)
5. [Quick Verification](#5-quick-verification)
6. [Detailed Example: SiO2 Network Analysis](#6-detailed-example-sio2-network-analysis)
7. [Output Files and Expected Results](#7-output-files-and-expected-results)
8. [Clustering Strategies Reference](#8-clustering-strategies-reference)
9. [Analyzers Reference](#9-analyzers-reference)
10. [Configuration Reference](#10-configuration-reference)
11. [Notes on Performance](#11-notes-on-performance)

---

## 1. Overview

`nexus-cat` is a Python package for network and cluster analysis of atomistic simulation trajectories. It is designed for researchers working with molecular dynamics (MD) or Monte Carlo (MC) data, with a particular focus on percolation theory applied to disordered materials (e.g., silicate glasses).

The single entry point is `nexus.main(settings)`. All configuration is handled through a `SettingsBuilder` object, which validates all parameters before execution. Results are exported to CSV files in a user-specified directory.

**Core capabilities:**

- Four clustering strategies: distance-based, bond-based (three-node bridging), coordination-number-constrained, and shared-neighbor-constrained.
- Nine analyzers: average cluster size, largest cluster size, spanning cluster size, gyration radius, correlation length, percolation probability, order parameter, concentration, and cluster size distribution.
- Two trajectory formats: extended XYZ and LAMMPS dump.
- Percolation detection via the period-vector algorithm (Livraghi et al., J. Chem. Theory Comput. 2021, 17, 4025).
- JIT-compiled geometry kernels (Numba) and generator-based frame iteration for large-scale trajectories.

---

## 2. System Requirements

### Hardware

- Any modern CPU with at least 4 GB of RAM for the example datasets.
- No GPU required.

### Software

| Dependency | Minimum version | Purpose |
|---|---|---|
| Python | 3.9 | Runtime |
| NumPy | 1.20 | Array operations |
| SciPy | 1.12 | KD-tree spatial queries |
| Numba | 0.53 | JIT-compiled geometry kernels |
| tqdm | 4.50 | Progress bars |
| colorama | 0.4.4 | Terminal color output |
| psutil | 5.8 | Memory and CPU monitoring |

All dependencies are automatically installed by pip (see Section 3).

**Note on Numba:** On first execution, Numba compiles the JIT-decorated functions and caches them on disk. This causes a one-time overhead of a few seconds, which will not occur on subsequent runs.

---

## 3. Installation

### Option A — Install from PyPI (recommended for editors)

```
pip install nexus-cat
```

To upgrade to the reviewed version:

```
pip install nexus-cat==0.1.0
```

### Option B — Install from source

```
git clone https://github.com/jperradin/nexus.git
cd nexus
pip install -e .
```

### Verify the installation

```python
import nexus
print(nexus.__version__)
```

Expected output: `0.1.0`

### Recommended: use a virtual environment

```
python -m venv nexus-env
source nexus-env/bin/activate      # Linux / macOS
nexus-env\Scripts\activate         # Windows
pip install nexus-cat
```

---

## 4. Repository Structure

```
nexus/
├── README.md
├── README_for_Editors.md          # This file
├── pyproject.toml                 # Build configuration and dependencies
├── setup.py
│
├── examples/
│   └── inputs/
│       ├── quickstart-SiO2.py               # Main example script
│       ├── example-SiO2-1008at.xyz          # Small test trajectory (1008 atoms)
│       ├── example-SiO2-1008at-high_pressure.xyz
│       ├── example-SiO2-27216at.xyz         # Medium trajectory (27216 atoms)
│       ├── example-SiO2-27216at-low_pressure.xyz
│       ├── example-SiO2-96000at.xyz         # Large trajectory
│       └── example-SiO2-1056000at.xyz       # Very large trajectory
│
├── benchmarks/
│   ├── performance/
│   │   ├── benchmarks.py                    # Performance benchmark script
│   │   └── data/
│   │       ├── generate_simple_cubic_lattice.py
│   │       └── generate_intricated_simple_cubic_lattice.py
│   └── validation/
│       ├── validation-perco.py              # Validation against percolation theory
│       ├── data/
│       │   ├── generate.py                  # Generate site-percolation lattices
│       │   └── generate_single_frame.py
│       └── output/                          # Pre-computed validation results
│           ├── 20/, 25/, 30/, 35/, 40/, 45/, 50/
│           └── data_collapse/               # Finite-size scaling data
│
└── src/nexus/                     # Package source (see Section 4.1 below)
    ├── __init__.py
    ├── main.py
    ├── version.py
    ├── config/settings.py
    ├── core/
    ├── analysis/
    ├── io/
    └── utils/
```

### 4.1 Source package layout

```
src/nexus/
├── __init__.py                    # Public API: SettingsBuilder, main, __version__
├── main.py                        # Pipeline entry point
├── version.py                     # Package version string
│
├── config/
│   └── settings.py                # GeneralSettings, LatticeSettings,
│                                  # ClusteringSettings, AnalysisSettings,
│                                  # Settings, SettingsBuilder
│
├── core/
│   ├── node.py                    # Node (atom) dataclass with __slots__
│   ├── frame.py                   # Single trajectory frame
│   ├── cluster.py                 # Cluster with period-vector percolation detection
│   └── system.py                  # Lazy trajectory iterator
│
├── analysis/
│   ├── analyzer_factory.py
│   ├── strategy_factory.py
│   ├── analyzers/                 # Nine analyzer classes
│   └── strategies/                # Four clustering strategy classes + KD-tree searcher
│
├── io/
│   ├── parser/parser.py
│   ├── reader/                    # XYZReader, LAMMPSReader, ReaderFactory
│   └── writer/                    # ClustersWriter, LogsWriter, PerformanceWriter
│
└── utils/
    ├── aesthetics.py              # Terminal display
    ├── geometry.py                # Numba JIT geometry functions
    └── performance.py             # Resource tracking
```

---

## 5. Quick Verification

The following script confirms that the package imports correctly and that the settings validation machinery works. It requires no trajectory file.

```python
from nexus import SettingsBuilder
import nexus.config.settings as c

# Verify that invalid settings raise informative errors
try:
    settings = (
        SettingsBuilder()
        .with_general(c.GeneralSettings(
            project_name="test",
            export_directory="./test_out",
            file_location="nonexistent.xyz",
            apply_pbc=True,
        ))
        .with_lattice(c.LatticeSettings(apply_custom_lattice=False))
        .with_clustering(c.ClusteringSettings(
            criterion="distance",
            node_types=["Si"],
            node_masses=[28.0855],
            connectivity=["Si", "Si"],
            cutoffs=[c.Cutoff(type1="Si", type2="Si", distance=3.5)],
        ))
        .with_analysis(c.AnalysisSettings(with_average_cluster_size=True))
        .build()
    )
    print("Settings object built successfully.")
    print(f"Output directory: {settings.output_directory}")
    print(f"Criterion: {settings.clustering.criterion}")
except Exception as e:
    print(f"Unexpected error: {e}")
```

Expected output:
```
Settings object built successfully.
Output directory: ./test_out/test
Criterion: distance
```

---

## 6. Detailed Example: SiO2 Network Analysis

This is the primary example demonstrating all major features. It is located at `examples/inputs/quickstart-SiO2.py` and uses the trajectory `examples/inputs/example-SiO2-1008at.xyz` (1008 atoms, multiple frames).

Run from the repository root:

```
python examples/inputs/quickstart-SiO2.py
```

### What this example does

1. Reads an SiO2 molecular dynamics trajectory in extended XYZ format.
2. Applies the **Coordination strategy**: nodes are silicon (Si) atoms connected through bridging oxygen (O) atoms; Si nodes are filtered by their coordination number (4 to 6 oxygen neighbors).
3. Computes **pairwise** connectivities (SiO4-SiO4, SiO5-SiO5, SiO6-SiO6).
4. Runs **all nine analyzers** on each frame and exports ensemble averages.
5. Writes unwrapped cluster coordinates (XYZ format, one file per connectivity per frame).

### Script content (abbreviated)

```python
from nexus import SettingsBuilder, main
import nexus.config.settings as c

path = "./examples/inputs/example-SiO2-1008at.xyz"

config_general = c.GeneralSettings(
    project_name="quickstart-SiO2",
    export_directory="./examples/outputs",
    file_location=path,
    range_of_frames=(0, -1),        # process all frames
    apply_pbc=True,
    verbose=True,
    save_logs=True,
    save_performance=False,
)

config_clustering = c.ClusteringSettings(
    criterion="bond",
    node_types=["Si", "O"],
    node_masses=[28.0855, 15.9994],
    connectivity=["Si", "O", "Si"],
    cutoffs=[
        c.Cutoff(type1="Si", type2="Si", distance=3.50),
        c.Cutoff(type1="Si", type2="O",  distance=2.30),
        c.Cutoff(type1="O",  type2="O",  distance=3.05),
    ],
    with_coordination_number=True,
    coordination_mode="O",
    coordination_range=[4, 6],
    with_pairwise=True,
    with_printed_unwrapped_clusters=True,
    print_mode="connectivity",
)

config_analysis = c.AnalysisSettings(with_all=True)

settings = (
    SettingsBuilder()
    .with_general(config_general)
    .with_lattice(c.LatticeSettings(apply_custom_lattice=False))
    .with_clustering(config_clustering)
    .with_analysis(config_analysis)
    .build()
)

main(settings)
```

### Expected output tree

After execution, the following directory is created:

```
examples/outputs/quickstart-SiO2/
├── logs.txt
├── Average_cluster_size_SiO4-SiO4.dat
├── Average_cluster_size_SiO5-SiO5.dat
├── Average_cluster_size_SiO6-SiO6.dat
├── Largest_cluster_size_SiO4-SiO4.dat
├── ...  (one .dat file per analyzer per connectivity)
├── Percolation_probability_SiO4-SiO4.dat
├── Order_parameter_SiO4-SiO4.dat
├── Cluster_size_distribution_SiO4-SiO4.dat
└── clusters/                           # only if with_printed_unwrapped_clusters=True
    ├── frame_0_SiO4-SiO4.xyz
    ├── frame_0_SiO5-SiO5.xyz
    └── ...
```

---

## 7. Output Files and Expected Results

### Output format

All analyzer outputs are written as whitespace-delimited `.dat` files with a header row. The format is consistent across all analyzers:

**Average cluster size example** (`Average_cluster_size_SiO4-SiO4.dat`):
```
frame,average_cluster_size
0,12.45
1,13.02
...
ensemble_average,12.74
```

**Percolation probability example** (`Percolation_probability_1-1.dat`):
```
frame,percolation_probability_x,percolation_probability_y,percolation_probability_z,percolation_probability_xyz
0,1,1,1,1
1,1,1,1,1
...
ensemble_average,0.98,0.97,0.99,0.95
```

**Cluster size distribution example** (`Cluster_size_distribution_SiO4-SiO4.dat`):
```
cluster_size,count
1,245
2,87
3,42
...
```

### Logs file

When `save_logs=True`, a `logs.txt` file is written to the output directory containing the full settings configuration used for the run, the package version, and a timestamp.

### Unwrapped cluster coordinates

When `with_printed_unwrapped_clusters=True`, clusters are written in extended XYZ format with unwrapped (non-periodic) coordinates. These files can be visualized with OVITO, VMD, or any XYZ-compatible viewer.

---

## 8. Clustering Strategies Reference

Strategy selection is automatic based on the `ClusteringSettings` fields, following this priority:

| Priority | Strategy | Trigger condition |
|---|---|---|
| 1 | Coordination | `with_coordination_number=True` and `with_number_of_shared=False` |
| 2 | Shared | `with_number_of_shared=True` |
| 3 | Distance | `criterion="distance"` |
| 4 | Bonding | `criterion="bond"` (default fallback) |

### Distance strategy

Connects any two nodes of matching types within the specified cutoff distance. Connectivity must be a two-element list, e.g., `["Si", "Si"]` or `["1", "1"]`.

```python
config_clustering = c.ClusteringSettings(
    criterion="distance",
    node_types=["1"],
    node_masses=[1.0],
    connectivity=["1", "1"],
    cutoffs=[c.Cutoff(type1="1", type2="1", distance=1.1)],
)
```

### Bonding strategy

Connects two networking nodes (first and third elements of connectivity) if they share a common bridging node (second element) within the respective cutoffs. Connectivity must be a three-element list.

```python
config_clustering = c.ClusteringSettings(
    criterion="bond",
    node_types=["Si", "O"],
    node_masses=[28.0855, 15.9994],
    connectivity=["Si", "O", "Si"],
    cutoffs=[
        c.Cutoff(type1="Si", type2="O", distance=2.30),
    ],
)
```

### Coordination strategy

Extends the bonding strategy by constraining the coordination number of networking nodes. Accepted nodes must have a number of bridging neighbors within `coordination_range`.

Three sub-modes for connectivity labeling:
- `with_pairwise=True`: separate connectivities for each matching pair of coordination numbers (e.g., SiO4-SiO4, SiO5-SiO5).
- `with_mixing=True`: cross-coordination connectivities (e.g., SiO4-SiO5).
- `with_alternating=True`: alternating pairs (e.g., SiO4-SiO5, SiO5-SiO6).

### Shared strategy

Extends the coordination strategy by additionally requiring that two connected networking nodes share a minimum number of common bridging neighbors (`shared_threshold`). Use `shared_threshold_mode="exact"` (exactly N shared neighbors) or `"at_least"` (at least N shared neighbors). This distinguishes corner-sharing (1 shared), edge-sharing (2 shared), and face-sharing (3 shared) polyhedral linkages.

```python
config_clustering = c.ClusteringSettings(
    criterion="bond",
    node_types=["Si", "O"],
    node_masses=[28.0855, 15.9994],
    connectivity=["Si", "O", "Si"],
    cutoffs=[
        c.Cutoff(type1="Si", type2="O", distance=2.30),
    ],
    with_coordination_number=True,
    coordination_mode="O",
    coordination_range=[4, 4],
    with_pairwise=True,
    with_number_of_shared=True,
    shared_mode="O",
    shared_threshold=1,
    shared_threshold_mode="exact",     # corner-sharing tetrahedra only
)
```

---

## 9. Analyzers Reference

| Analyzer flag | Metric | Formula |
|---|---|---|
| `with_average_cluster_size` | Mean cluster size | sum(s^2 * n_s) / sum(s * n_s), excluding percolating clusters |
| `with_largest_cluster_size` | S_max | Size of the largest cluster per frame |
| `with_spanning_cluster_size` | S_span | Size of the largest non-percolating (finite) cluster |
| `with_gyration_radius` | R_g | RMS radius of gyration (from unwrapped coordinates), binned by cluster size |
| `with_correlation_length` | xi | Second moment of the gyration radius distribution |
| `with_percolation_probability` | Pi | Fraction of frames containing a 3D-spanning cluster |
| `with_order_parameter` | P_inf | Fraction of networking nodes in the percolating cluster, per dimension |
| `with_concentration` | c | Ratio of networking nodes to total nodes |
| `with_cluster_size_distribution` | n_s | Histogram of cluster sizes |

Enable all analyzers at once with `with_all=True` in `AnalysisSettings`.

Enable individual analyzers selectively:

```python
config_analysis = c.AnalysisSettings(
    with_percolation_probability=True,
    with_order_parameter=True,
    with_average_cluster_size=True,
)
```

---

## 10. Configuration Reference

### GeneralSettings

| Parameter | Type | Default | Description |
|---|---|---|---|
| `project_name` | str | "Project" | Used for the output subdirectory name |
| `export_directory` | str | "exports" | Root export directory |
| `file_location` | str | "" | Path to the trajectory file |
| `range_of_frames` | (int, int) | (0, -1) | Start and end frame indices; -1 means last frame |
| `apply_pbc` | bool | False | Apply periodic boundary conditions |
| `verbose` | bool | False | Print title, settings, and progress bars |
| `save_logs` | bool | False | Write logs.txt to the export directory |
| `save_performance` | bool | False | Write per-frame timing and memory metrics to JSON |

### LatticeSettings

| Parameter | Type | Default | Description |
|---|---|---|---|
| `apply_custom_lattice` | bool | False | Override lattice from file with a user-supplied matrix |
| `custom_lattice` | ndarray | zeros (3x3) | User-supplied 3x3 lattice matrix (Angstroms) |
| `get_lattice_from_file` | bool | False | Read lattice from a separate file |
| `apply_lattice_to_all_frames` | bool | True | Use the same lattice for every frame |

For trajectories in extended XYZ format (with `Lattice="..."` header) or LAMMPS dump format, set `apply_custom_lattice=False` and the lattice is read automatically per frame.

### ClusteringSettings (key parameters)

| Parameter | Type | Description |
|---|---|---|
| `criterion` | str | `"distance"` or `"bond"` |
| `node_types` | list of str | Atom type symbols to include |
| `node_masses` | list of float | Masses for each type (reduced units) |
| `connectivity` | list of str | 2-element (distance) or 3-element (bond) pattern |
| `cutoffs` | list of Cutoff | One `Cutoff(type1, type2, distance)` per relevant pair |
| `with_coordination_number` | bool | Enable coordination-number filtering |
| `coordination_mode` | str | Count neighbors of type: `"all_types"`, `"same_type"`, `"different_type"`, or a specific type symbol |
| `coordination_range` | [int, int] | Accepted coordination number range (inclusive) |
| `with_pairwise` | bool | Separate results per matching coordination pair |
| `with_mixing` | bool | Separate results per cross-coordination pair |
| `with_alternating` | bool | Separate results per adjacent coordination pair |
| `with_number_of_shared` | bool | Enable shared-neighbor filtering |
| `shared_threshold` | int | Minimum number of shared bridging neighbors |
| `shared_threshold_mode` | str | `"exact"` or `"at_least"` |
| `with_printed_unwrapped_clusters` | bool | Write cluster XYZ files |
| `print_mode` | str | `"all"`, `"connectivity"`, `"individual"`, or `"none"` |

### SettingsBuilder validation

`SettingsBuilder.build()` raises `ValueError` with an informative message if any of the following constraints are violated (non-exhaustive list):

- `criterion="bond"` with a connectivity list of length other than 3.
- `criterion="distance"` with a connectivity list of length other than 2.
- `with_coordination_number=True` with an invalid `coordination_mode` or an empty/reversed `coordination_range`.
- `with_pairwise`, `with_mixing`, or `with_alternating` set to `True` without `with_coordination_number=True`.
- `with_number_of_shared=True` without `with_coordination_number=True`.
- `node_types` and `node_masses` lists of different lengths.

---

## 11. Notes on Performance

- **Numba JIT compilation:** Hot-path functions in `src/nexus/utils/geometry.py` are decorated with `@njit`. Compilation occurs once on first execution and is cached in a `__pycache__` directory. Subsequent runs are not affected.

- **Memory:** Frame iteration is generator-based. Only one frame is held in memory at a time, regardless of trajectory length. The 1008-atom example uses less than 100 MB of RAM. The 1 056 000-atom trajectory (`example-SiO2-1056000at.xyz`) can be processed on a machine with 16 GB of RAM.

- **Performance tracking:** Set `save_performance=True` in `GeneralSettings` to record per-frame execution time, memory usage, and CPU usage. Results are written as a JSON file to the export directory.

- **Parallelism:** The main pipeline (`nexus.main`) is single-threaded. The validation script (`benchmarks/validation/validation-perco.py`) demonstrates how to parallelize multiple independent analyses using Python's `multiprocessing.Pool`.
