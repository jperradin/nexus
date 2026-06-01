# Performance benchmarks

Wall-time and memory scaling of the four `nexus-cat` clustering strategies
(`distance`, `bond`, `coordination`, `shared`) on synthetic binary cubic
lattices, plotted as a function of the number of sites N.

## Run environment

Benchmarks executed on **2026-05-28**.

| Component | Value |
|-----------|-------|
| CPU | 12th Gen Intel(R) Core(TM) i5-12600KF (16 logical CPUs) |
| Memory | 16 GB (16279632 kB) |
| OS | Ubuntu 24.04.4 LTS (Linux 6.6.114.1-microsoft-standard-WSL2, x86_64) |
| Python | 3.12.3 |

Software versions:

| Package | Version |
|---------|---------|
| numpy | 2.3.5 |
| scipy | 1.17.0 |
| numba | 0.63.1 |
| matplotlib | 3.10.9 |
| iminuit | 2.32.0 |

## Figures

- `performance_execution_time.png` — execution time vs N (log-log), with
  iminuit power-law fits `t = A * N**alpha`; the fitted exponent `alpha`
  (empirical algorithmic complexity) is shown per strategy in the legend.
- `performance_memory.png` — memory usage vs N (log-log), markers only, no fits.

Both figures have two stacked panels: `p < pc` (top) and `p > pc` (bottom).
The percolation regime is inferred from the site concentration `c = N / L**3`
(`c < 0.3116` -> `p < pc`, `c >= 0.3116` -> `p > pc`), so for a given lattice size the
denser dataset is the `p > pc` one.

## Regenerate

```bash
cd benchmarks/performance/data
python generate_simple_cubic_lattice.py   # writes benchmark-<N>-<L>.xyz 

cd ..
python benchmarks.py        # runs the four strategies, writes output/<strategy>_strategy/...
python plot_performance.py  # reads output/, writes both figures above
```
