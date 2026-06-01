#!/usr/bin/env python3
"""Extract critical exponents for 3D site-percolation validation and fill the
Grace project with the latest results.

Reproduces ``benchmarks/validation/validation.png`` from the *current* export
data found in this directory. Two analyses are run for three observables:

    * correlation_length   (xi)      -> panels a, b, e
    * average_cluster_size (<S>)     -> panels c, f
    * spanning_cluster_size (S_max)  -> panels d, g

1. Scaling laws (panels b, c, d)
   Peak value of each observable is taken at every system size L and fitted to
   a power law  y = a * L**b. Results go to ``scaling_laws/``.

2. Data collapse (panels e, f, g)
   Finite-size-scaling collapse with the Kawashima-Ito quality function
   (``fssa``) optimised with ``iminuit``. The free parameters are the critical
   occupation probability ``p_c``, the correlation-length exponent ``nu`` and
   the observable exponent ``zeta`` (so that y = L**(-zeta/nu) * observable and
   x = L**(1/nu) * (p - p_c)). Results go to ``data_collapse/``.

Data layout (identical control-parameter grid for every size)::

    <L>/export/probability-<observable>.dat            # col0 = p, col1 = value
    <L>/export/probability-<observable>--errors.dat    # col0 = p, col1 = error

Outputs:
    * ``validation.png``  -- the matplotlib figure (next to this script)
    * ``validation.agr``      -- a filled copy of ``validation.agr``

For the Grace project the original ``.agr`` is left untouched; only the numeric
data blocks, the three collapse exponent annotations and the three scaling-law
legends are rewritten, everything else (layout, fonts, colours) is preserved
verbatim.

Graph -> panel mapping (read from the .agr axis labels):

    G0 -> e)  xi  data collapse
    G1 -> f)  <S> data collapse
    G2 -> g)  S_max data collapse
    G3 -> b)  xi  scaling law      (S0 = data, S1 = fit line)
    G4 -> c)  <S> scaling law
    G5 -> d)  S_max scaling law
    G6 -> a)  xi vs p (raw curves)
"""

import os
import shutil

import fssa
import numpy as np
import matplotlib.pyplot as plt
from iminuit import Minuit
from iminuit.cost import LeastSquares

# --------------------------------------------------------------------------- #
# Configuration
# --------------------------------------------------------------------------- #
HERE = os.path.dirname(os.path.abspath(__file__))

# Linear system sizes L (one export directory each).
SIZES = np.array([20, 25, 30, 35, 40, 45, 50])

# observable -> (file stem, latex symbol, power-law (a0, b0),
#                collapse zeta (init, lo, hi), reference exponent string)
OBSERVABLES = {
    "correlation_length": {
        "stem": "probability-correlation_length",
        "symbol": r"$\xi$",
        "pl_init": (1.0, 1.0),
        "zeta": (0.9, 0.4, 1.5),
        "ref": r"$L^{1.00(5)}$",
    },
    "average_cluster_size": {
        "stem": "probability-average_cluster_size",
        "symbol": r"$\langle S\rangle$",
        "pl_init": (1.0, 2.0),
        "zeta": (1.75, 1.0, 2.6),
        "ref": r"$L^{2.1(1)}$",
    },
    "spanning_cluster_size": {
        "stem": "probability-spanning_cluster_size",
        "symbol": r"$S_{max}$",
        "pl_init": (1.0, 2.5),
        "zeta": (1.95, 1.0, 2.8),
        "ref": r"$L^{2.54(9)}$",
    },
}

# Critical region used for the collapse fit (window around the xi peak).
P_C_GUESS = 0.3116          # simple-cubic site percolation threshold
P_C_THEORY = 0.3116         # theoretical p_c used by the fixed-p_c scaling method
CRITICAL_HALF_WIDTH = 0.04  # keep p in [p_c - w, p_c + w] for the collapse
PARABOLA_HALF_WIN = 3       # +/- points around argmax used for the parabola fit

NU_LIMITS = (0.70, 1.15)
RHO_C_LIMITS = (0.300, 0.325)

SCALING_DIR = os.path.join(HERE, "scaling_laws")
COLLAPSE_DIR = os.path.join(HERE, "data_collapse")

# Grace project: source template and filled output (parent validation dir).
SRC = os.path.join(os.path.dirname(HERE), "validation.agr")
DST = os.path.join(os.path.dirname(HERE), "validation_latest.agr")

# graph letter -> observable name (for the 7-set graphs)
COLLAPSE_GRAPH = {"G0": "correlation_length",
                  "G1": "average_cluster_size",
                  "G2": "spanning_cluster_size"}
SCALING_GRAPH = {"G3": "correlation_length",
                 "G4": "average_cluster_size",
                 "G5": "spanning_cluster_size"}


# --------------------------------------------------------------------------- #
# Data loading
# --------------------------------------------------------------------------- #
def load_observable(stem):
    """Return (p, values, errors) stacked over all sizes.

    ``p`` is the shared 1D grid; ``values`` and ``errors`` are (n_sizes, n_p).
    """
    p_ref = None
    values, errors = [], []
    for L in SIZES:
        base = os.path.join(HERE, str(L), "export")
        data = np.loadtxt(os.path.join(base, stem + ".dat"), comments="#")
        err = np.loadtxt(os.path.join(base, stem + "--errors.dat"), comments="#")
        p = data[:, 0]
        if p_ref is None:
            p_ref = p
        elif not np.allclose(p, p_ref):
            raise ValueError(f"control-parameter grid differs for L={L}")
        values.append(data[:, 1])
        # clamp tiny/zero errors so the quality function stays finite
        errors.append(np.maximum(err[:, 1], 1e-3))
    return p_ref, np.array(values), np.array(errors)


# --------------------------------------------------------------------------- #
# 1. Scaling laws  (peak value vs L, power-law fit)
# --------------------------------------------------------------------------- #
def power_law(L, a, b):
    return a * L ** b


def parabola_peak(p, y, dy, half=PARABOLA_HALF_WIN):
    """Refine the peak of a noisy curve by a local quadratic fit.

    A parabola is fitted to the +/-``half`` grid points around ``argmax`` and
    its vertex is taken as the peak. Falls back to the raw ``argmax`` bin if the
    parabola is not concave or its vertex falls outside the fit window (e.g. a
    monotonic / edge curve). The error is read at the nearest grid point.
    """
    idx = int(np.argmax(y))
    lo = max(idx - half, 0)
    hi = min(idx + half + 1, len(p))
    coef = np.polyfit(p[lo:hi], y[lo:hi], 2)        # coef = [a2, a1, a0]
    p_star, y_star = p[idx], y[idx]                 # default: raw bin
    if coef[0] < 0:                                 # concave -> real maximum
        vertex = -coef[1] / (2 * coef[0])
        if p[lo] <= vertex <= p[hi - 1]:
            p_star, y_star = vertex, np.polyval(coef, vertex)
    j = int(np.argmin(np.abs(p - p_star)))
    return p_star, y_star, dy[j]


def value_at_pc(p, y, dy, pc=P_C_THEORY):
    """Observable value at the grid point closest to the theoretical p_c."""
    j = int(np.argmin(np.abs(p - pc)))
    return p[j], y[j], dy[j]


def _fit_power(p_pts, y_pts, dy_pts, init):
    a0, b0 = init
    cost = LeastSquares(SIZES.astype(float), y_pts, dy_pts, power_law)
    m = Minuit(cost, a=a0, b=b0)
    m.limits["a"] = (0, None)
    m.migrad()
    m.hesse()
    return dict(p_pts=p_pts, y_pts=y_pts, dy_pts=dy_pts,
                a=m.values["a"], da=m.errors["a"],
                b=m.values["b"], db=m.errors["b"],
                rchi2=m.fval / max(len(SIZES) - 2, 1))


def fit_scaling_law(name, info):
    p, values, errors = load_observable(info["stem"])

    # two peak-picking strategies, one point per size each
    n = len(SIZES)
    par = dict(p=np.zeros(n), y=np.zeros(n), dy=np.zeros(n))
    pcm = dict(p=np.zeros(n), y=np.zeros(n), dy=np.zeros(n))
    for i in range(n):
        par["p"][i], par["y"][i], par["dy"][i] = parabola_peak(
            p, values[i], errors[i])
        pcm["p"][i], pcm["y"][i], pcm["dy"][i] = value_at_pc(
            p, values[i], errors[i])

    parabola = _fit_power(par["p"], par["y"], par["dy"], info["pl_init"])
    fixed_pc = _fit_power(pcm["p"], pcm["y"], pcm["dy"], info["pl_init"])

    os.makedirs(SCALING_DIR, exist_ok=True)
    out = os.path.join(SCALING_DIR, f"{name}_scaling.dat")
    with open(out, "w") as f:
        f.write(f"# scaling law for {name}:  y = a * L^b\n")
        f.write(f"# parabola peak : a = {parabola['a']:.6g} +/- {parabola['da']:.3g}"
                f"  b = {parabola['b']:.6g} +/- {parabola['db']:.3g}"
                f"  rchi2 = {parabola['rchi2']:.4g}\n")
        f.write(f"# fixed p_c={P_C_THEORY} : a = {fixed_pc['a']:.6g} +/- {fixed_pc['da']:.3g}"
                f"  b = {fixed_pc['b']:.6g} +/- {fixed_pc['db']:.3g}"
                f"  rchi2 = {fixed_pc['rchi2']:.4g}\n")
        f.write("# L  y_parab  dy_parab  p_parab  y_fit_parab"
                "  y_pc  dy_pc  y_fit_pc\n")
        for i in range(n):
            f.write(f"{SIZES[i]:d} {par['y'][i]:.6g} {par['dy'][i]:.6g} "
                    f"{par['p'][i]:.6g} {power_law(SIZES[i], parabola['a'], parabola['b']):.6g} "
                    f"{pcm['y'][i]:.6g} {pcm['dy'][i]:.6g} "
                    f"{power_law(SIZES[i], fixed_pc['a'], fixed_pc['b']):.6g}\n")

    print(f"[scaling] {name:22s} parabola b = {parabola['b']:.3f} +/- {parabola['db']:.3f}"
          f"   fixed-p_c b = {fixed_pc['b']:.3f} +/- {fixed_pc['db']:.3f}")
    # parabola is the default method; expose both, keep back-compat aliases
    return dict(p=p, values=values, errors=errors,
                parabola=parabola, fixed_pc=fixed_pc,
                p_peak=par["p"], y_peak=par["y"], dy_peak=par["dy"],
                a=parabola["a"], b=parabola["b"], db=parabola["db"])


# --------------------------------------------------------------------------- #
# 2. Data collapse  (Kawashima-Ito quality via fssa + iminuit)
# --------------------------------------------------------------------------- #
class KIQuality:
    """Kawashima-Ito collapse quality as an iminuit cost function."""

    def __init__(self, l, rho, a, da):
        self.l, self.rho, self.a, self.da = l, rho, a, da

    def __call__(self, rho_c, nu, zeta):
        scaled = fssa.scaledata(self.l, self.rho, self.a, self.da,
                                rho_c, nu, zeta)
        return fssa.quality(scaled.x, scaled.y, scaled.dy)


def fit_data_collapse(name, info, scaling):
    p = scaling["p"]
    values = scaling["values"]
    errors = scaling["errors"]

    # restrict to the critical window so the collapse is well conditioned
    mask = np.abs(p - P_C_GUESS) <= CRITICAL_HALF_WIDTH
    rho = p[mask]
    a = values[:, mask]
    da = errors[:, mask]

    z0, zlo, zhi = info["zeta"]
    cost = KIQuality(SIZES.astype(float), rho, a, da)
    m = Minuit(cost, rho_c=P_C_GUESS, nu=0.88, zeta=z0)
    m.limits = [RHO_C_LIMITS, NU_LIMITS, (zlo, zhi)]
    m.fixed["rho_c"] = True  # fix p_c to the theoretical value for the main result

    # iterate migrad until a valid minimum is reached
    for _ in range(50):
        m.migrad(ncall=100000, iterate=10, use_simplex=True)
        m.hesse()
        if m.valid:
            break
        m.values["nu"] = np.random.uniform(*NU_LIMITS)
        m.values["zeta"] = np.random.uniform(zlo, zhi)

    rho_c, nu, zeta = m.values["rho_c"], m.values["nu"], m.values["zeta"]
    err = m.errors
    scaled = fssa.scaledata(SIZES.astype(float), rho, a, da, rho_c, nu, zeta)

    os.makedirs(COLLAPSE_DIR, exist_ok=True)
    params = os.path.join(COLLAPSE_DIR, f"{name}_params.dat")
    with open(params, "w") as f:
        f.write(f"# data-collapse exponents for {name}\n")
        f.write(f"p_c   {rho_c:.6f} {err['rho_c']:.6f}\n")
        f.write(f"nu    {nu:.6f} {err['nu']:.6f}\n")
        f.write(f"zeta  {zeta:.6f} {err['zeta']:.6f}\n")
        f.write(f"zeta_over_nu  {zeta / nu:.6f}\n")
        f.write(f"quality  {m.fval:.6f}\n")
        f.write(f"valid  {m.valid}\n")

    # one collapsed curve per size
    for i, L in enumerate(SIZES):
        out = os.path.join(COLLAPSE_DIR, f"{name}_collapse_L{L}.dat")
        np.savetxt(out, np.column_stack([scaled.x[i], scaled.y[i], scaled.dy[i]]),
                   header="x = L^(1/nu)(p-p_c)   y = L^(-zeta/nu) obs   dy")

    print(f"[collapse] {name:21s} p_c={rho_c:.4f} nu={nu:.3f} "
          f"zeta={zeta:.3f} zeta/nu={zeta / nu:.3f} Q={m.fval:.3f} "
          f"valid={m.valid}")
    return dict(scaled=scaled, rho_c=rho_c, nu=nu, zeta=zeta,
                nu_err=err["nu"], quality=m.fval)


# --------------------------------------------------------------------------- #
# 3. Figure (reproduce validation.png layout)
# --------------------------------------------------------------------------- #
def make_figure(scaling_res, collapse_res):
    names = list(OBSERVABLES.keys())
    palette = plt.cm.viridis(np.linspace(0.15, 0.85, len(SIZES)))
    panel_color = {"correlation_length": "tab:blue",
                   "average_cluster_size": "tab:green",
                   "spanning_cluster_size": "tab:purple"}

    fig = plt.figure(figsize=(15, 7))
    gs = fig.add_gridspec(2, 4, width_ratios=[1.4, 1, 1, 1], hspace=0.32,
                          wspace=0.33)

    # panel a) raw correlation length vs p
    ax_a = fig.add_subplot(gs[:, 0])
    sc = scaling_res["correlation_length"]
    for i, L in enumerate(SIZES):
        ax_a.plot(sc["p"], sc["values"][i], "-", color=palette[i], lw=1,
                  label=str(L))
    ax_a.axvline(collapse_res["correlation_length"]["rho_c"], color="k",
                 ls="--", lw=0.8, label=r"fitted $p_c$")
    ax_a.axvline(P_C_THEORY, color="r", ls=":", lw=1.2,
                 label=rf"theory $p_c$={P_C_THEORY}")
    ax_a.set_xlabel("Occ. probability $p$")
    ax_a.set_ylabel(r"Correlation length $\xi$ [$\AA$]")
    ax_a.legend(title="$L$", fontsize=7, ncol=2)
    ax_a.text(0.04, 0.95, "a)", transform=ax_a.transAxes, weight="bold")

    # panels b, c, d) scaling laws (log-log): parabola peak vs fixed-p_c
    Lf = np.linspace(SIZES.min(), SIZES.max(), 100)
    for j, name in enumerate(names):
        ax = fig.add_subplot(gs[0, j + 1])
        r = scaling_res[name]
        col = panel_color[name]

        par, pcm = r["parabola"], r["fixed_pc"]
        # parabola-refined peak (filled) + fit
        ax.errorbar(SIZES, par["y_pts"], yerr=par["dy_pts"], fmt="o", mfc=col,
                    mec="k", color=col, ms=6, label="parabola peak")
        ax.plot(Lf, power_law(Lf, par["a"], par["b"]), "-", color=col,
                label=rf"$L^{{{par['b']:.2f}({int(round(par['db'] * 100)):d})}}$")
        # value at theoretical p_c (open red) + fit
        ax.errorbar(SIZES, pcm["y_pts"], yerr=pcm["dy_pts"], fmt="s", mfc="w",
                    mec="r", color="r", ms=5, label=rf"at $p_c$={P_C_THEORY}")
        ax.plot(Lf, power_law(Lf, pcm["a"], pcm["b"]), "--", color="r",
                label=rf"$L^{{{pcm['b']:.2f}({int(round(pcm['db'] * 100)):d})}}$")

        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("$L$")
        ax.set_ylabel(OBSERVABLES[name]["symbol"] + "$(p_c)$")
        ax.legend(fontsize=6.5)
        ax.text(0.05, 0.9, "bcd"[j] + ")", transform=ax.transAxes,
                weight="bold")

    # panels e, f, g) data collapse
    for j, name in enumerate(names):
        ax = fig.add_subplot(gs[1, j + 1])
        c = collapse_res[name]
        scaled = c["scaled"]
        for i, L in enumerate(SIZES):
            ax.errorbar(scaled.x[i], scaled.y[i], yerr=scaled.dy[i], fmt="-o",
                        ms=2, lw=0.6, color=palette[i], label=str(L))
        ax.set_xlabel(r"$L^{1/\nu}(p-p_c)$")
        ax.set_ylabel(OBSERVABLES[name]["symbol"] + r"$\,L^{-\zeta/\nu}$")
        txt = (rf"$\nu$ = {c['nu']:.2f}" + "\n"
               + rf"$\zeta/\nu$ = {c['zeta'] / c['nu']:.2f}")
        ax.text(0.04, 0.96, txt, transform=ax.transAxes, va="top", fontsize=8)
        ax.text(0.84, 0.9, "efg"[j] + ")", transform=ax.transAxes,
                weight="bold")

    out = os.path.join(HERE, "validation.png")
    fig.savefig(out, dpi=200, bbox_inches="tight")
    print(f"\nfigure written: {out}")


# --------------------------------------------------------------------------- #
# 4. Grace (.agr) project filling
# --------------------------------------------------------------------------- #
def fmt_row(values):
    return " ".join(f"{v:.6g}" for v in values)


def uncert(b, db):
    """Render b +/- db as the compact ``b(d)`` notation, e.g. 2.43(7)."""
    return f"{b:.2f}({int(round(db * 100)):d})"


def build_blocks(scaling, collapse):
    """Return {(graph, set_index): [data row strings]}."""
    blocks = {}
    sizes = SIZES

    # collapse panels (e, f, g): one set per size, scaled (x, y, dy)
    for g, name in COLLAPSE_GRAPH.items():
        scaled = collapse[name]["scaled"]
        for i in range(len(sizes)):
            rows = [fmt_row((x, y, dy)) for x, y, dy
                    in zip(scaled.x[i], scaled.y[i], scaled.dy[i])
                    if np.isfinite(x) and np.isfinite(y)]
            blocks[(g, i)] = rows

    # scaling panels (b, c, d): S0 = peak data, S1 = fitted power law
    for g, name in SCALING_GRAPH.items():
        r = scaling[name]
        blocks[(g, 0)] = [fmt_row((L, y, dy)) for L, y, dy
                          in zip(sizes, r["y_peak"], r["dy_peak"])]
        Lf = np.linspace(sizes.min(), sizes.max(), 7)
        blocks[(g, 1)] = [fmt_row((L, power_law(L, r["a"], r["b"])))
                          for L in Lf]

    # panel a (G6): raw xi vs p, one set per size
    xi = scaling["correlation_length"]
    for i in range(len(sizes)):
        blocks[("G6", i)] = [fmt_row((p, y, dy)) for p, y, dy
                             in zip(xi["p"], xi["values"][i], xi["errors"][i])]
    return blocks


def rewrite_data(lines, blocks):
    """Replace every ``@target Gx.Sy`` data block with new rows."""
    out = []
    i = 0
    n = len(lines)
    while i < n:
        line = lines[i]
        out.append(line)
        if line.startswith("@target "):
            graph, sset = line.split()[1].split(".")          # e.g. G0, S3
            key = (graph, int(sset[1:]))
            i += 1
            # keep the @type line as-is
            if i < n and lines[i].startswith("@type"):
                out.append(lines[i])
                i += 1
            if key in blocks:                                  # swap data rows
                out.extend(blocks[key])
                while i < n and lines[i].strip() != "&":
                    i += 1
            # fall through to copy the closing '&' (and any kept old rows)
            continue
        i += 1
    return out


def patch_annotations(text, scaling, collapse):
    """Update collapse exponent strings and scaling-law legends in place."""
    c_xi = collapse["correlation_length"]
    c_s = collapse["average_cluster_size"]
    c_sm = collapse["spanning_cluster_size"]

    replacements = [
        # panel g) S_max: nu / 1s(=zeta) / 1sn(=zeta/nu)
        (r'@    string def "\xn \c;\C 0.88\f{}\n\x1/s \c;\C 2.22\f{}\n\x1/sn \c;\C 2.53"',
         f'@    string def "\\xn \\c;\\C {c_sm["nu"]:.2f}\\f{{}}\\n'
         f'\\x1/s \\c;\\C {c_sm["zeta"]:.2f}\\f{{}}\\n'
         f'\\x1/sn \\c;\\C {c_sm["zeta"] / c_sm["nu"]:.2f}"'),
        # panel f) <S>: nu / gamma(=zeta) / gamma/nu
        (r'@    string def "\xn \c;\C 0.85\f{}\n\xg \c;\C 1.75\f{}\n\xg/n \c;\C 2.06"',
         f'@    string def "\\xn \\c;\\C {c_s["nu"]:.2f}\\f{{}}\\n'
         f'\\xg \\c;\\C {c_s["zeta"]:.2f}\\f{{}}\\n'
         f'\\xg/n \\c;\\C {c_s["zeta"] / c_s["nu"]:.2f}"'),
        # panel e) xi: nu / nu'(=zeta) / nu'/nu
        ("@    string def \"\\xn \\c;\\C 0.88\\f{}\\n\\xn\\f{}'\\x \\c;\\C 0.90\\f{}\\n\\xn\\f{}'\\x/n \\c;\\C 1.02\"",
         f'@    string def "\\xn \\c;\\C {c_xi["nu"]:.2f}\\f{{}}\\n'
         f"\\xn\\f{{}}'\\x \\c;\\C {c_xi['zeta']:.2f}\\f{{}}\\n"
         f"\\xn\\f{{}}'\\x/n \\c;\\C {c_xi['zeta'] / c_xi['nu']:.2f}\""),
        # scaling-law fit legends
        ('@    s1 legend  "L\\S1.00(5)"',
         f'@    s1 legend  "L\\S{uncert(scaling["correlation_length"]["b"], scaling["correlation_length"]["db"])}"'),
        ('@    s1 legend  "L\\S2.1(1)"',
         f'@    s1 legend  "L\\S{uncert(scaling["average_cluster_size"]["b"], scaling["average_cluster_size"]["db"])}"'),
        ('@    s1 legend  "L\\S2.54(9)"',
         f'@    s1 legend  "L\\S{uncert(scaling["spanning_cluster_size"]["b"], scaling["spanning_cluster_size"]["db"])}"'),
    ]

    for old, new in replacements:
        if old not in text:
            print(f"  WARNING: annotation not found, skipped:\n    {old!r}")
            continue
        text = text.replace(old, new)
    return text


def fill_agr(scaling, collapse):
    """Write ``validation_latest.agr`` from the template + fresh results."""
    shutil.copyfile(SRC, DST)
    with open(DST) as f:
        lines = [ln.rstrip("\n") for ln in f]

    blocks = build_blocks(scaling, collapse)
    lines = rewrite_data(lines, blocks)
    text = "\n".join(lines) + "\n"
    text = patch_annotations(text, scaling, collapse)

    with open(DST, "w") as f:
        f.write(text)

    print(f"\nfilled project written: {DST}")


# --------------------------------------------------------------------------- #
def main():
    np.random.seed(0)  # collapse fit converges first try; keep it reproducible

    print("=== scaling laws ===")
    scaling_res = {name: fit_scaling_law(name, info)
                   for name, info in OBSERVABLES.items()}

    print("\n=== data collapse ===")
    collapse_res = {name: fit_data_collapse(name, info, scaling_res[name])
                    for name, info in OBSERVABLES.items()}

    make_figure(scaling_res, collapse_res)
    fill_agr(scaling_res, collapse_res)


if __name__ == "__main__":
    main()
