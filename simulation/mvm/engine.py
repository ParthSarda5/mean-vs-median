"""Monte Carlo engine: simulate one (distribution, n) cell and summarise every estimator.

Each cell draws B samples in N_BATCHES batches from its own RNG stream (derived from the
master seed and a CRC of the cell label, so results do not depend on the order cells run in).
Monte Carlo standard errors come from the spread of the per-batch statistics.
"""
from __future__ import annotations

import sys
import time
import zlib
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd

from . import config as C
from .distributions import Dist
from .estimators import ESTIMATORS


def cell_rng(label: str, n: int, salt: int = 0) -> np.random.Generator:
    ss = np.random.SeedSequence([C.SEED, zlib.crc32(f"{label}|{n}".encode()), salt])
    return np.random.default_rng(ss)


def simulate_cell(dist: Dist, n: int, B: int, est_names, n_batches: int = C.N_BATCHES) -> list[dict]:
    """Return one summary row per estimator for B Monte Carlo samples of size n."""
    n_batches = min(n_batches, B)
    bs = B // n_batches
    B = bs * n_batches
    rng = cell_rng(dist.name, n)
    T = {k: np.empty(B) for k in est_names}
    for b in range(n_batches):
        X = dist.sample(rng, (bs, n))
        for k in est_names:
            T[k][b * bs:(b + 1) * bs] = ESTIMATORS[k](X)

    theta = dist.theta
    sq = np.sqrt(n_batches)
    var_all = {k: T[k].var(ddof=1) for k in est_names}
    var_b = {k: T[k].reshape(n_batches, bs).var(axis=1, ddof=1) for k in est_names}
    rows = []
    for k in est_names:
        t = T[k].reshape(n_batches, bs)
        bias_b = t.mean(axis=1) - theta
        mse_b = ((t - theta) ** 2).mean(axis=1)
        ae = np.abs(T[k] - theta)
        q = np.quantile(ae, C.QUANTILES)
        row = dict(
            dist=dist.name, kind=dist.kind, nu=dist.nu, eps=dist.eps, n=n, B=B, estimator=k,
            bias=bias_b.mean(), se_bias=bias_b.std(ddof=1) / sq,
            var=var_all[k], se_var=var_b[k].std(ddof=1) / sq,
            mse=mse_b.mean(), se_mse=mse_b.std(ddof=1) / sq,
            n_var=n * var_all[k], n_mse=n * mse_b.mean(),
            q50_abs=q[0], q90_abs=q[1], q99_abs=q[2],
        )
        if "mean" in T:   # efficiency of every estimator relative to the sample mean
            ratio_b = var_b["mean"] / var_b[k]
            row["var_ratio_vs_mean"] = var_all["mean"] / var_all[k]
            row["se_var_ratio_vs_mean"] = ratio_b.std(ddof=1) / sq
        rows.append(row)
    return rows


def _work(args):
    dist, n, B, est_names = args
    t0 = time.time()
    rows = simulate_cell(dist, n, B, est_names)
    return rows, dist.name, n, time.time() - t0


def run_cells(cells, B: int, est_names, workers: int = 2, tag: str = "") -> pd.DataFrame:
    """cells: iterable of (Dist, n).  Runs them (in parallel) and returns one tidy DataFrame."""
    jobs = [(d, n, B, tuple(est_names)) for d, n in cells]
    out, t0 = [], time.time()
    with ProcessPoolExecutor(max_workers=workers) as ex:
        for i, (rows, name, n, dt) in enumerate(ex.map(_work, jobs), 1):
            out.extend(rows)
            print(f"[{tag}] {i:3d}/{len(jobs)} {name:18s} n={n:5d}  {dt:6.1f}s  (elapsed {time.time()-t0:6.0f}s)",
                  file=sys.stderr, flush=True)
    return pd.DataFrame(out)
