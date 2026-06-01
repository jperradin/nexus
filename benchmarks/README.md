# Benchmarks

This directory contains two independent benchmark suites used to characterize `Nexus-CAT`:

- `performance/` — wall-time and memory scaling of the four clustering strategies on synthetic cubic lattices.
- `validation/` — correctness check against the exact results of 3D site-percolation theory on a simple cubic lattice.

Each suite ships a `*-reference` directory (and a matching `.tar.gz`) holding the pre-computed results used in the manuscript, so the figures can be reproduced without re-running the analysis. `performance/` also has its own `README.md` documenting the exact run environment (CPU, OS, package versions).

The full benchmark datasets (generated lattices and reference outputs) are archived on Zenodo: [https://doi.org/10.5281/zenodo.20487537](https://doi.org/10.5281/zenodo.20487537).

For any issue, question, or remark, please contact the corresponding author at julien.perradin@protonmail.com.

---

## Directory layout

```
benchmarks/
├── performance/
│   ├── README.md                              # Run environment + figure description
│   ├── benchmarks.py                          # Runs all four strategies on each dataset
│   ├── plot_performance.py                    # Reads output/, writes the two figures
│   ├── data/
│   │   └── generate_simple_cubic_lattice.py   # Generate binary cubic lattices
│   ├── output-reference/                      # Pre-computed per-strategy performance JSON
│   ├── export-reference/                      # Post-processed scaling curves + fits (.dat)
│   ├── performance_execution_time.png         # Published figure: time vs N (log-log)
│   ├── performance_memory.png                 # Published figure: memory vs N (log-log)
│   └── performance-122.agr / .eps             # xmgrace project + export
│
└── validation/
    ├── validation-perco.py                    # Runs nexus on all site-percolation lattices
    ├── data/
    │   ├── generate.py                        # Generate site-percolation lattices (all p, all L)
    │   └── generate_single_frame.py           # Single-frame variant (quick checks)
    ├── nexus_inputs_all                       # Input list: (path, label, L) for all p
    ├── nexus_inputs_critical_point            # Input list: only p = p_c = 0.3116
    ├── output/
    │   └── extract_exponents_and_plot.py      # Post-processing: scaling laws + data collapse
    ├── output-reference/                      # Pre-computed results
    │   ├── 20/, 25/, ..., 50/                 #   one subdir per lattice size L (export/ + p<...>/)
    │   ├── data_collapse/                     #   finite-size-scaling collapse curves
    │   ├── scaling_laws/                      #   critical-exponent fits
    │   └── validation_reproduced.png          #   matplotlib figure
    ├── validation.agr                         # xmgrace template (left untouched)
    └── validation_latest.agr / .eps / .png    # filled xmgrace project + exports
```

---

## 1. Performance benchmarks

The performance suite measures the cost of each clustering strategy (`distance`, `bond`, `coordination`, `shared`) on binary `A–B` cubic lattices of increasing size (N atoms, cubic side L). See `performance/README.md` for the exact machine and package versions used for the published numbers.

### 1.1 Generate the benchmark datasets

```bash
cd benchmarks/performance/data
python generate_simple_cubic_lattice.py   # writes benchmark-<N>-<L>.xyz
```

The script emits extended-XYZ files named `benchmark-<N>-<L>.xyz` (where `<N>` is the atom count and `<L>` is the cubic side in reduced units) into the current directory, plus a serialized `benchmark-<N>-<L>.npz` artifact holding the raw `networking_sites` / `bridge_sites` coordinate arrays. `benchmarks.py` reads the XYZ files from `./data/` relative to its own location, so no further move is needed.

Generation is **deterministic**: the RNG is seeded with a hardcoded `seed = 1_000_000 * L + N`, so re-running reproduces byte-identical datasets.

### 1.2 Run the benchmarks

```bash
cd benchmarks/performance
python benchmarks.py        # runs the four strategies, writes output/<strategy>_strategy/...
python plot_performance.py  # reads output/, writes both figures
```

For each dataset listed in the `benchmarks` list at the bottom of `benchmarks.py`, the script runs the four clustering strategies sequentially with `save_performance=True`. Per-frame wall-time, RSS memory, and CPU usage are written as JSON files to:

```
benchmarks/performance/output/<strategy>_strategy/benchmark-<N>-<L>/performance_*.json
```

To enable or disable a dataset, comment/un-comment the corresponding tuple in the `benchmarks` list.

### 1.3 Plots and published figures

`plot_performance.py` reads the per-run JSON files back and produces two log-log figures, each with two stacked panels (`p < pc` top, `p > pc` bottom):

- `performance_execution_time.png` — execution time vs N, with iminuit power-law fits `t = A * N**alpha`; the fitted exponent `alpha` (empirical algorithmic complexity) is shown per strategy in the legend.
- `performance_memory.png` — memory usage vs N, markers only, no fits.

The percolation regime is inferred from the site concentration `c = N / L**3`, so for a given lattice size the denser dataset is the `p > pc` one. `performance-122.agr` / `.eps` are the xmgrace project file and export used in the manuscript. Pre-computed JSON and the post-processed scaling curves live in `output-reference/` and `export-reference/`.

---

## 2. Validation against percolation theory

The validation suite benchmarks `nexus-cat` against the analytical predictions of 3D site-percolation on a simple cubic lattice (critical probability p_c ≈ 0.3116, universal critical exponents β, γ, ν).

### 2.1 Generate the percolation lattices

```bash
cd benchmarks/validation/data
python generate.py                  # all probabilities p ∈ [0.20, 0.40] (step 0.002), sizes L ∈ {20, 25, ..., 50}
python generate_single_frame.py     # single-frame variant, useful for quick checks
```

`generate.py` writes, directly into `benchmarks/validation/data/`:

- `./<L>/percolation_sites_0.3116.xyz` — 100 independent frames at p = p_c for each L.
- `./<L>/percolation_sites_<p>.xyz` — sweep over p (`p_index` 0..100, `p = 0.2 + 0.002 * p_index`) for finite-size scaling.
- `./<L>/percolation_sites_L<L>.npz` — one serialized artifact per L holding every generated frame as an `(M, 3)` int32 array of occupied-site coordinates (keys `pc_n###` and `p###_n###`), plus the `probabilities` grid.

**Reproducibility:** every frame is generated from a hardcoded, unique-per-frame seed `seed = 1_000_000 * L + 1000 * p_index + n` (the critical point uses the reserved `p_index = 500`). Because generation runs inside a `numba` nopython function, the seed is set with `np.random.seed()` *inside* that function — numba keeps its own RNG state, independent of NumPy's Python-level RNG, so seeding from plain Python would have no effect. The XYZ files are written in overwrite mode, so re-running reproduces identical datasets rather than appending.

### 2.2 Run the analysis

Two input lists are provided:

- `nexus_inputs_critical_point` — only p = p_c, for computing percolation probability and order parameter at the transition.
- `nexus_inputs_all` — full sweep, required for finite-size scaling and critical-exponent fits.

Both lists use paths relative to `benchmarks/validation/` (e.g. `./data/<L>/...`). To switch between them, edit the filename passed to `np.loadtxt` in `validation-perco.py` (default: `nexus_inputs_critical_point`).

```bash
cd benchmarks/validation
python validation-perco.py
```

The script parallelises over 4 workers (`multiprocessing.Pool`); adjust `num_workers` inside the file if needed. Each trajectory is analysed with the Distance strategy (single type `"1"`, cutoff = 1.1) and all analyzers enabled (`with_all=True`). Results are written to `benchmarks/validation/output/<L>/`, with the per-observable curves under `output/<L>/export/` as `probability-<observable>.dat` (plus `--errors` / `--std` companions).

### 2.3 Post-processing: critical exponents and figure

```bash
cd benchmarks/validation/output
python extract_exponents_and_plot.py
```

`extract_exponents_and_plot.py` consumes the per-size `export/` curves and runs two analyses for three observables (`correlation_length` ξ, `average_cluster_size` ⟨S⟩, `spanning_cluster_size` S_max):

1. **Scaling laws** — peak value of each observable vs L fitted to a power law `y = a * L**b`; written to `scaling_laws/`.
2. **Data collapse** — finite-size-scaling collapse with the Kawashima-Ito quality function (`fssa`) optimised with `iminuit` over `p_c`, `nu`, `zeta`; written to `data_collapse/`.

It then writes the matplotlib figure `validation_reproduced.png` and fills a copy of the xmgrace template, leaving `validation.agr` untouched and emitting `validation_latest.agr` (only data blocks, collapse-exponent annotations, and scaling-law legends are rewritten; layout/fonts/colours are preserved). This single script supersedes the former `extract_exponents.py` + `fill_agr.py` pair.

### 2.4 Pre-computed reference results

`output-reference/20/` … `output-reference/50/` contain the reference runs used in the manuscript (each with its `export/` curves and per-probability subdirs). `output-reference/data_collapse/` and `output-reference/scaling_laws/` hold the post-processed finite-size-scaling curves and critical-exponent fits, alongside `output-reference/validation_reproduced.png`.

### 2.5 Expected outcome

At p = p_c and in the thermodynamic limit the analysis should reproduce:

- Percolation probability Π(p_c) → ≈ 0.5 (crossing point for all L).
- Order parameter P_∞ ∼ (p − p_c)^β with β ≈ 0.4181.
- Correlation length ξ ∼ |p − p_c|^(−ν) with ν ≈ 0.8764.
- Mean cluster size χ ∼ |p − p_c|^(−γ) with γ ≈ 1.7951.

Data collapse of the rescaled quantities vs. (p − p_c) L^(1/ν) for L ∈ {20, …, 50} provides the visual validation (see `validation_reproduced.png` / `validation_latest.png`).

---

## Notes

- All generation scripts use `numba` JIT compilation; the first run incurs a one-time compilation overhead.
- All randomized generators are deterministically seeded (see the per-suite reproducibility notes above), and each lattice size `L` is also dumped to a serialized `.npz` artifact so the exact dataset can be reloaded with `np.load(...)` without regeneration. Seed-based reproduction assumes the same `numpy`/`numba` versions; the `.npz` artifacts are version-independent.
- Every script uses paths relative to its own directory — always `cd` into the subfolder shown in each command block before running.
- The `*-reference` directories (and their `.tar.gz` archives) are the frozen manuscript results; the live `output/` directories are overwritten by re-runs.
- Memory footprint for the largest performance dataset (`benchmark-600000-100`) stays below 4 GB thanks to the generator-based frame iterator.
