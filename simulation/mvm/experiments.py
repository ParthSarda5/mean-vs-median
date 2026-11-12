"""The experiments of the simulation study.  Each function writes a tidy CSV to results/data."""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats

from . import config as C
from .distributions import Dist
from .engine import cell_rng, run_cells
from .estimators import ALL, ESTIMATORS, MAIN, running_median


def _scaled(x: int, quick: bool, floor: int = 200) -> int:
    return max(floor, x // C.QUICK_FACTOR) if quick else x


def _tukey(eps):
    return Dist("tukey", eps=eps, scale=C.TUKEY_SCALE)


def _shift(eps):
    return Dist("shift", eps=eps, shift=C.SHIFT_LOC)


def main_dists():
    ds = [Dist("normal")] + [Dist("t", nu=v) for v in C.T_NUS] + [Dist("laplace")]
    ds += [_tukey(e) for e in C.EPS_GRID] + [_shift(e) for e in C.EPS_GRID]
    return ds


def exp_main(quick=False, workers=2):
    """Mean & median over the full (distribution x n) grid with B = 1e5 replications."""
    B = _scaled(C.B_MAIN, quick, 2000)
    ns = C.N_GRID
    cells = [(d, n) for d in main_dists() for n in ns]
    df = run_cells(cells, B, MAIN, workers, "main")
    df.to_csv(C.DATA_DIR / "main_grid.csv", index=False)
    return df


def exp_t_sweep(quick=False, workers=2):
    """Simulated Var(mean)/Var(median) over a fine grid of t degrees of freedom."""
    B = _scaled(C.B_SWEEP, quick, 2000)
    ns = C.T_SWEEP_N
    cells = [(Dist("t", nu=v), n) for v in C.T_SWEEP_NUS for n in ns]
    df = run_cells(cells, B, MAIN, workers, "t-sweep")
    df.to_csv(C.DATA_DIR / "t_sweep.csv", index=False)
    return df


def exp_comparators(quick=False, workers=2):
    """Five estimators on selected distributions plus the contamination sweeps (n = 100)."""
    B = _scaled(C.B_COMP, quick, 2000)
    sel = [Dist("normal"), Dist("t", nu=5), Dist("t", nu=3), Dist("laplace"), Dist("t", nu=1),
           _tukey(0.10), _shift(0.10)]
    cells = {(d.name, n): (d, n) for d in sel for n in C.N_COMP}
    for e in C.EPS_SWEEP:
        for d in (_tukey(e), _shift(e)):
            cells[(d.name, C.N_SWEEP_CONTAM)] = (d, C.N_SWEEP_CONTAM)
    df = run_cells(list(cells.values()), B, ALL, workers, "comparators")
    df.to_csv(C.DATA_DIR / "comparators.csv", index=False)
    return df


def exp_cauchy_paths(quick=False):
    """Running mean and running median along a few independent Cauchy streams."""
    N = _scaled(C.CAUCHY_PATH_LEN, quick, 2000)
    checkpoints = np.unique(np.r_[np.arange(1, 101), np.geomspace(100, N, 400).astype(int)])
    rows = []
    for p in range(C.CAUCHY_N_PATHS):
        x = cell_rng("cauchy-path", N, salt=p).standard_cauchy(N)
        rmean = np.cumsum(x) / np.arange(1, N + 1)
        rmed = running_median(x)
        for n in checkpoints:
            rows.append(dict(path=p, n=int(n), running_mean=rmean[n - 1], running_median=rmed[n - 1]))
    df = pd.DataFrame(rows)
    df.to_csv(C.DATA_DIR / "cauchy_paths.csv", index=False)
    return df


def exp_breakdown(quick=False):
    """Replace k of n = 101 clean N(0,1) points by 1e6; track the median |estimate| over replications."""
    n, reps = C.BREAKDOWN_N, _scaled(C.BREAKDOWN_REPS, quick, 100)
    X0 = cell_rng("breakdown", n).standard_normal((reps, n))
    rows = []
    for k in range(0, int(0.6 * n) + 1):
        X = X0.copy()
        X[:, :k] = C.BREAKDOWN_VALUE
        for name, f in ESTIMATORS.items():
            est = f(X)
            rows.append(dict(k=k, frac=k / n, estimator=name, median_abs=float(np.median(np.abs(est)))))
    df = pd.DataFrame(rows)
    df.to_csv(C.DATA_DIR / "breakdown.csv", index=False)
    return df


def exp_sensitivity(quick=False, save=True):
    """Sensitivity curves SC_n(y) = n [T(x_1..x_{n-1}, y) - T(x_1..x_{n-1})] with x_i = Phi^{-1}((i-.5)/(n-1))."""
    n = C.SC_N
    x = stats.norm.ppf((np.arange(1, n) - 0.5) / (n - 1))
    ys = np.linspace(-C.SC_YMAX, C.SC_YMAX, C.SC_NPOINTS)
    X = np.tile(x, (len(ys), 1))
    X = np.column_stack([X, ys])
    base = x[None, :]
    out = {"y": ys}
    for name, f in ESTIMATORS.items():
        out[name] = n * (f(X) - f(base)[0])
    df = pd.DataFrame(out)
    if save:
        df.to_csv(C.DATA_DIR / "sensitivity.csv", index=False)
    return df


ALL_EXPERIMENTS = {
    "main": exp_main,
    "tsweep": exp_t_sweep,
    "comparators": exp_comparators,
    "cauchy": exp_cauchy_paths,
    "breakdown": exp_breakdown,
    "sensitivity": exp_sensitivity,
}
