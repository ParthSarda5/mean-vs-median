"""Vectorised location estimators.  Each maps an array of shape (B, n) to B estimates."""
from __future__ import annotations

import numpy as np
from scipy import stats

HUBER_K = 1.345          # tuning constant: ~95% asymptotic efficiency at the normal
TRIM_PROP = 0.10         # proportion trimmed from each tail


def mean_est(X):
    return X.mean(axis=1)


def median_est(X):
    return np.median(X, axis=1)


def trimmed_mean(X, prop: float = TRIM_PROP):
    return stats.trim_mean(X, prop, axis=1)


def huber_location(X, k: float = HUBER_K, maxit: int = 60, tol: float = 1e-9):
    """Huber M-estimator of location with MAD scale, by iteratively reweighted averaging.

    Starts at the median; each step reweights w_i = min(1, k / |r_i|), r_i = (x_i - mu)/s,
    with the (fixed) scale s = 1.4826 * MAD.  Rows with MAD = 0 return the median.
    """
    mu = np.median(X, axis=1)
    s = 1.4826 * np.median(np.abs(X - mu[:, None]), axis=1)
    ok = s > 1e-12
    s = np.where(ok, s, 1.0)
    for _ in range(maxit):
        r = np.abs(X - mu[:, None]) / s[:, None]
        w = np.minimum(1.0, k / np.maximum(r, 1e-300))
        mu_new = (w * X).sum(axis=1) / w.sum(axis=1)
        mu_new = np.where(ok, mu_new, mu)
        done = np.max(np.abs(mu_new - mu) / s) < tol
        mu = mu_new
        if done:
            break
    return mu


def hodges_lehmann(X, max_elems: int = 20_000_000):
    """Hodges-Lehmann estimator: median of all Walsh averages (x_i + x_j)/2, i <= j."""
    B, n = X.shape
    i, j = np.triu_indices(n)
    out = np.empty(B)
    step = max(1, max_elems // len(i))
    for s in range(0, B, step):
        blk = X[s:s + step]
        out[s:s + step] = np.median(0.5 * (blk[:, i] + blk[:, j]), axis=1)
    return out


ESTIMATORS = {
    "mean": mean_est,
    "median": median_est,
    "trim10": trimmed_mean,
    "huber": huber_location,
    "hl": hodges_lehmann,
}

LABELS = {
    "mean": "Sample mean",
    "median": "Sample median",
    "trim10": "10% trimmed mean",
    "huber": "Huber ($k{=}1.345$)",
    "hl": "Hodges–Lehmann",
}
SHORT = {"mean": "Mean", "median": "Median", "trim10": "Trim-10%", "huber": "Huber", "hl": "HL"}
MAIN = ("mean", "median")
ALL = tuple(ESTIMATORS)


def running_median(x):
    """Median of x[:1], x[:2], ..., x[:N] in O(N log N) with two heaps."""
    import heapq
    lo, hi = [], []          # max-heap (negated) of the lower half, min-heap of the upper half
    out = np.empty(len(x))
    for t, v in enumerate(x):
        if not lo or v <= -lo[0]:
            heapq.heappush(lo, -v)
        else:
            heapq.heappush(hi, v)
        if len(lo) > len(hi) + 1:
            heapq.heappush(hi, -heapq.heappop(lo))
        elif len(hi) > len(lo):
            heapq.heappush(lo, -heapq.heappop(hi))
        out[t] = -lo[0] if len(lo) > len(hi) else 0.5 * (-lo[0] + hi[0])
    return out
