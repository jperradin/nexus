"""Plot performance scaling vs number of sites (log-log).

Reads the per-run performance JSON files written by ``benchmarks.py`` under
``output/<strategy>_strategy/benchmark-<N>-<L>/performance_*.json`` and produces
two figures, each with two stacked panels (p < pc top, p > pc bottom):

    performance_execution_time.png  -- execution time vs N, with power-law fits
    performance_memory.png          -- memory usage vs N, no fits

Both axes are logarithmic. For the execution-time figure a power law
t = A * N**alpha appears as a straight line whose slope ``alpha`` is the
empirical algorithmic complexity exponent, fitted with iminuit (MINUIT least
squares in log-log space) and reported in each legend entry.

The percolation regime is inferred from the site concentration c = N / L**3:
c < 0.5 is the dilute (p < pc) regime, c >= 0.5 is the dense (p > pc) regime,
i.e. for a given lattice size the dataset with more nodes is the p > pc one.
"""

import json
import re
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FuncFormatter
from iminuit import Minuit
from iminuit.cost import LeastSquares

HERE = Path(__file__).resolve().parent
OUTPUT_DIR = HERE / "output"
EXPORT_DIR = HERE / "export"
FIG_TIME_PATH = HERE / "performance_execution_time.png"
FIG_MEM_PATH = HERE / "performance_memory.png"

BENCH_RE = re.compile(r"benchmark-(\d+)-(\d+)")

# Fixed order/styling so the legend is stable across runs.
STRATEGY_STYLE = {
    "distance": ("tab:blue", "o"),
    "bond": ("tab:orange", "s"),
    "coordination": ("tab:green", "^"),
    "shared": ("tab:red", "D"),
}

REGIMES = ["p < pc", "p > pc"]  # row order (top, bottom)


def regime_of(num_sites: int, lattice_size: int) -> str:
    """Classify a dataset by site concentration c = N / L**3."""
    concentration = num_sites / lattice_size**3
    return "p < pc" if concentration < 0.5 else "p > pc"


def collect():
    """Return data[strategy][regime] -> sorted list of (N, time_s, mem_mb)."""
    # Aggregate (mean) over any repeated runs for the same (strategy, regime, N).
    acc = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))

    for strat_dir in sorted(OUTPUT_DIR.glob("*_strategy")):
        strategy = strat_dir.name.replace("_strategy", "")
        for bench_dir in strat_dir.glob("benchmark-*"):
            m = BENCH_RE.search(bench_dir.name)
            if not m:
                continue
            num_sites, lattice_size = int(m.group(1)), int(m.group(2))
            regime = regime_of(num_sites, lattice_size)

            for jf in bench_dir.glob("performance_*.json"):
                with open(jf) as f:
                    rec = json.load(f)
                t_ms = rec.get("execution_time_ms")
                mem = rec.get("memory_usage_mb")
                if t_ms is None or mem is None:
                    continue
                n_nodes = rec.get("metrics", {}).get("number_nodes", num_sites)
                acc[strategy][regime][n_nodes].append((t_ms / 1000.0, mem))

    data = {}
    for strategy, by_regime in acc.items():
        data[strategy] = {}
        for regime, by_n in by_regime.items():
            series = []
            for n_nodes, runs in sorted(by_n.items()):
                t_mean = sum(t for t, _ in runs) / len(runs)
                m_mean = sum(mm for _, mm in runs) / len(runs)
                series.append((n_nodes, t_mean, m_mean))
            data[strategy][regime] = series
    return data


def _line(x, c, alpha):
    """Straight line in log-log space: log10(y) = c + alpha * log10(N)."""
    return c + alpha * x


def fit_power_law(ns, ys):
    """Fit y = A * N**alpha with iminuit. Returns (alpha, alpha_err, c) or None."""
    ns = np.asarray(ns, dtype=float)
    ys = np.asarray(ys, dtype=float)
    mask = (ns > 0) & (ys > 0)
    ns, ys = ns[mask], ys[mask]
    if ns.size < 2:
        return None

    log_n = np.log10(ns)
    log_y = np.log10(ys)
    # Unweighted fit in log space (no per-point uncertainties available).
    least_squares = LeastSquares(log_n, log_y, np.ones_like(log_y), _line)
    minuit = Minuit(least_squares, c=0.0, alpha=1.0)
    minuit.migrad()
    return minuit.values["alpha"], minuit.errors["alpha"], minuit.values["c"]


def _slug(text: str) -> str:
    """Filesystem-safe token (e.g. 'p < pc' -> 'p_lt_pc')."""
    return (
        text.replace("<", "lt")
        .replace(">", "gt")
        .replace(" ", "_")
        .strip("_")
    )


def _write_dat(path, header, columns):
    """Write whitespace-separated columns to ``path`` with a comment header."""
    rows = np.column_stack([np.asarray(col, dtype=float) for col in columns])
    EXPORT_DIR.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        f.write(f"# {header}\n")
        for row in rows:
            f.write(" ".join(f"{v:.10g}" for v in row) + "\n")
    print(f"Saved {path}")


def _plot_panel(ax, data, regime, metric_index, ylabel, fit=True, metric_slug="metric"):
    """Plot one panel: y[metric_index] vs N (log-log).

    With ``fit=True`` markers are overlaid with a dashed power-law fit and the
    fitted exponent is shown in the legend. With ``fit=False`` the points are
    simply connected by a solid line (no fit).

    Each curve is also exported to a ``.dat`` file under ``EXPORT_DIR``: one raw
    file per (metric, regime, strategy), plus a matching ``*-fit.dat`` holding
    the fitted power-law curve whenever a fit was performed.
    """
    regime_slug = _slug(regime)
    for strategy, (color, marker) in STRATEGY_STYLE.items():
        series = data.get(strategy, {}).get(regime, [])
        if not series:
            continue
        ns = [s[0] for s in series]
        ys = [s[metric_index] for s in series]

        stem = f"{metric_slug}-{regime_slug}-{strategy}"
        _write_dat(
            EXPORT_DIR / f"{stem}.dat",
            f"strategy={strategy} regime='{regime}' raw data\n# N {metric_slug}",
            (ns, ys),
        )

        if fit:
            result = fit_power_law(ns, ys)
            if result is not None:
                alpha, alpha_err, c = result
                label = rf"{strategy} ($\alpha$={alpha:.2f}$\pm${alpha_err:.2f})"
                x_fit = np.array([min(ns), max(ns)], dtype=float)
                y_fit = 10 ** (c + alpha * np.log10(x_fit))
                ax.plot(
                    x_fit, y_fit, color=color, linestyle="--", linewidth=1, alpha=0.7
                )
                x_dense = np.logspace(
                    np.log10(min(ns)), np.log10(max(ns)), 100
                )
                y_dense = 10 ** (c + alpha * np.log10(x_dense))
                _write_dat(
                    EXPORT_DIR / f"{stem}-fit.dat",
                    f"strategy={strategy} regime='{regime}' "
                    f"power-law fit y=A*N**alpha "
                    f"alpha={alpha:.6g} alpha_err={alpha_err:.6g} c={c:.6g} "
                    f"(A=10**c={10 ** c:.6g})\n# N {metric_slug}_fit",
                    (x_dense, y_dense),
                )
            else:
                label = strategy
            ax.plot(ns, ys, color=color, marker=marker, linestyle="none", label=label)
        else:
            ax.plot(ns, ys, color=color, marker=marker, linestyle="-", label=strategy)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_ylabel(ylabel)
    ax.grid(True, which="both", alpha=0.3)


def _make_figure(
    data,
    metric_index,
    ylabel,
    suptitle,
    path,
    fit,
    legend_title,
    metric_slug,
    plain_yticks=False,
):
    """Build one figure: two stacked panels (p<pc top, p>pc bottom) vs N."""
    fig, axes = plt.subplots(2, 1, figsize=(8, 9), sharex=True)
    plain = FuncFormatter(lambda y, _: f"{y:.0f}")
    for row, regime in enumerate(REGIMES):
        ax = axes[row]
        _plot_panel(
            ax, data, regime, metric_index, ylabel, fit=fit, metric_slug=metric_slug
        )
        ax.set_title(regime, fontweight="bold")
        ax.legend(title=legend_title, fontsize=8)
        if plain_yticks:
            # Show plain integers (200, 300, ...) instead of 2x10^2 on log axis.
            ax.yaxis.set_major_formatter(plain)
            ax.yaxis.set_minor_formatter(plain)
    axes[1].set_xlabel("Number of sites N")
    fig.suptitle(suptitle, fontweight="bold")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    print(f"Saved {path}")


def main():
    data = collect()
    if not data:
        raise SystemExit(f"No performance JSON found under {OUTPUT_DIR}")

    _make_figure(
        data,
        metric_index=1,
        ylabel="Execution time (s)",
        suptitle="nexus-cat execution-time complexity vs number of sites (log-log)",
        path=FIG_TIME_PATH,
        fit=True,
        legend_title="strategy (slope = complexity)",
        metric_slug="execution_time_s",
    )
    _make_figure(
        data,
        metric_index=2,
        ylabel="Memory usage (MB)",
        suptitle="nexus-cat memory usage vs number of sites (log-log)",
        path=FIG_MEM_PATH,
        fit=False,
        legend_title="strategy",
        metric_slug="memory_mb",
        plain_yticks=True,
    )


if __name__ == "__main__":
    main()
