# NEXUS-CAT — Cluster Analysis Toolkit
## Code Verification Guide for Editors

**Package name:** nexus-cat  
**Version:** 0.1.1  
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
5. [Testing and Verification with Examples](#5-testing-and-verification-with-examples)

For strategies, analyzers, configuration options, and performance notes, see the full documentation: https://nexus-cat.readthedocs.io/en/latest/

---

## 1. Overview

`nexus-cat` is a Python package for network and cluster analysis of atomistic simulation trajectories. It is designed for researchers working with molecular dynamics (MD) or Monte Carlo (MC) data, with a particular focus on percolation theory applied to disordered systems (e.g., silicate glasses).

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

### Option A — Install from PyPI (version 0.1.1)

```
pip install nexus-cat==0.1.1
```

### Option B — Install from source (recommended for editors)

```bash
tar -xzvf nexus_cat-0.1.1-archive_for_Editors.tar.gz
cd nexus
pip install -e .
```

### Verify the installation

```python
import nexus
print(nexus.__version__)
```

Expected output: `0.1.1`

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
│   ├── inputs/
│   │   ├── README.md                           # Full dataset documentation
│   │   ├── quickstart-a_SiO2.py                # Amorphous silica analysis
│   │   ├── quickstart-a_H2O.py                 # Amorphous water (LD/HD/VHD)
│   │   ├── quickstart-a_Si.py                  # Amorphous silicon (LD/HD/VHD)
│   │   ├── example-a_SiO2-8064at-{0,5,10,...,35}GPa.xyz   # 8064 atoms, 10 frames
│   │   ├── example-a_H2O-8192mol-{0.01,0.5,1}kbar.xyz     # 8192 molecules, 10 frames
│   │   └── example-a_Si-100000at-{0,12}GPa.xyz            # 100000 atoms, 1 frame
│   └── outputs/                                    # Created at runtime
│       └── verification-JPerradin-20260417.tar.gz  # Contain validated results for verification
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


## 5. Testing and Verification with Examples

The bundled examples, quickstart scripts, dataset descriptions, step-by-step run instructions, and the expected output tree are fully documented in:

**→ [`examples/inputs/README.md`](examples/inputs/README.md)**

In short:

```bash
cd examples/inputs
python quickstart-a_SiO2.py   # amorphous silica (SiO4/SiO5/SiO6 + Stishovite)
python quickstart-a_H2O.py    # amorphous water   (LD / HD / VHD)
python quickstart-a_Si.py     # amorphous silicon (LD / HD / VHD)
```

Outputs are written to `examples/outputs/<trajectory-name>/`. Reference results for comparison are packaged in `examples/outputs/verification-JPerradin-20260417.tar.gz`. Extract in place with `tar -xzvf examples/outputs/verification-JPerradin-20260417.tar.gz -C examples/outputs/` and diff against your own `examples/outputs/<trajectory-name>/`.
---

## Further documentation

Strategies, analyzers, full configuration reference, and performance notes are maintained in the online documentation: https://nexus-cat.readthedocs.io/en/latest/
