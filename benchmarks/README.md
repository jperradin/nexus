# Benchmarks

This directory contains two independent benchmark suites used to characterize `nexus-cat`:

- `performance/` — wall-time and memory scaling of the four clustering strategies on synthetic cubic lattices.
- `validation/` — correctness check against the exact results of 3D site-percolation theory on a simple cubic lattice.

For any issue, question, or remark, please contact the corresponding author at julien.perradin@protonmail.com.

---

## Directory layout

```
benchmarks/
├── performance/
│   ├── benchmarks.py                           # Runs all four strategies on each dataset
│   ├── data/
│   │   ├── generate_simple_cubic_lattice.py
│   │   └── generate_intricated_simple_cubic_lattice.py
│   ├── performance-121.agr / .eps              # xmgrace plots (published figures)
│   └── performance-122.agr / .eps
│
└── validation/
    ├── validation-perco.py                     # Runs nexus on all site-percolation lattices
    ├── data/
    │   ├── generate.py                         # Generate site-percolation lattices (all p, all L)
    │   └── generate_single_frame.py
    ├── nexus_inputs_all                        # Input list: (path, label, L) for all p
    ├── nexus_inputs_critical_point             # Input list: only p = p_c = 0.3116
    ├── output/
    │   ├── 20/, 25/, ..., 50/                  # Pre-computed results, one subdir per lattice size L
    │   ├── data_collapse/                      # Finite-size scaling data
    │   └── scaling_laws/                       # Critical-exponent fits
    ├── snapshot.png                            # Sample percolating cluster
    └── validation.agr / .eps / .jpg            # xmgrace plot (published figure)
```

---

## 1. Performance benchmarks

The performance suite measures the cost of each clustering strategy (`distance`, `bond`, `coordination`, `shared`) on binary `A–B` cubic lattices of increasing size (N atoms, cubic side L).

### 1.1 Generate the benchmark datasets

From the repository root:

```bash
cd benchmarks/performance/data
python generate_simple_cubic_lattice.py             # writes universal_benchmark-<N>-<L>.xyz
python generate_intricated_simple_cubic_lattice.py  # writes an interpenetrated variant
```

Each script emits extended-XYZ files named `universal_benchmark-<N>-<L>.xyz` where `<N>` is the atom count and `<L>` is the cubic side (reduced units). Move or symlink the generated files to `benchmarks/data/` so that the hard-coded paths inside `benchmarks.py` resolve.

### 1.2 Run the benchmarks

```bash
# From the repository root
python benchmarks/performance/benchmarks.py
```

For each dataset listed in the `benchmarks` list at the bottom of `benchmarks.py`, the script runs the four clustering strategies sequentially with `save_performance=True`. Per-frame wall-time, RSS memory, and CPU usage are written as JSON files to:

```
benchmarks/output/<strategy>/universal_benchmark-<N>-<L>/performance*.json
```

To enable or disable a dataset, comment/un-comment the corresponding tuple in the `benchmarks` list.

### 1.3 Published figures

`performance-121.agr` / `performance-122.agr` are the xmgrace project files used to produce Figures 12.1 and 12.2 of the accompanying manuscript; the `.eps` files are the exports.

---

## 2. Validation against percolation theory

The validation suite benchmarks `nexus-cat` against the analytical predictions of 3D site-percolation on a simple cubic lattice (critical probability p_c ≈ 0.3116, universal critical exponents β, γ, ν).

### 2.1 Generate the percolation lattices

```bash
cd benchmarks/validation/data
python generate.py                  # all probabilities p ∈ [0.20, 0.40], sizes L ∈ {20, 25, ..., 50}
python generate_single_frame.py     # single-frame variant, useful for quick checks
```

`generate.py` writes:

- `./<L>/percolation_sites_0.3116.xyz` — 100 independent frames at p = p_c for each L.
- `./<L>/percolation_sites_<p>.xyz` — sweep over p for finite-size scaling (step 0.002).

Move the generated `<L>/` directories to `benchmarks/validation/data/` so that the paths in `nexus_inputs_*` resolve.

### 2.2 Run the analysis

Two input lists are provided:

- `nexus_inputs_critical_point` — only p = p_c, for computing percolation probability and order parameter at the transition.
- `nexus_inputs_all` — full sweep, required for finite-size scaling and critical-exponent fits. **Note:** by default these paths point to `benchmarks/validation_percolation/data/...`; edit the list to match `benchmarks/validation/data/...` before running.

```bash
# From the repository root
python benchmarks/validation/validation-perco.py
```

The script parallelises over 4 workers (`multiprocessing.Pool`); adjust `num_workers` inside the file if needed. Each trajectory is analysed with the Distance strategy (single type, cutoff = 1.1) and all analyzers enabled. Results are written to `benchmarks/validation/output/<L>/<label>/`.

### 2.3 Pre-computed reference results

`output/20/` … `output/50/` contain reference runs used in the manuscript. `output/data_collapse/` and `output/scaling_laws/` contain the post-processed finite-size-scaling curves and critical-exponent fits that reproduce the published figure `validation.jpg` / `validation.eps`.

### 2.4 Expected outcome

At p = p_c and in the thermodynamic limit the analysis should reproduce:

- Percolation probability Π(p_c) → ≈ 0.5 (crossing point for all L).
- Order parameter P_∞ ∼ (p − p_c)^β with β ≈ 0.4181.
- Correlation length ξ ∼ |p − p_c|^(−ν) with ν ≈ 0.8764.
- Mean cluster size χ ∼ |p − p_c|^(−γ) with γ ≈ 1.7951.

Data collapse of the rescaled quantities vs. (p − p_c) L^(1/ν) for L ∈ {20, …, 50} provides the visual validation (see `validation.jpg`).

---

## Notes

- All generation scripts use `numba` JIT compilation; the first run incurs a one-time compilation overhead.
- Hard-coded relative paths inside `benchmarks.py` and `validation-perco.py` assume execution from the repository root, not from within the `benchmarks/` subdirectory.
- Memory footprint for the largest performance dataset (`universal_benchmark-600000-100`) stays below 4 GB thanks to the generator-based frame iterator.
