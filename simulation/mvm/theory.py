"""Closed-form / numerically exact asymptotic results used to check the simulation."""
from __future__ import annotations

import numpy as np
from scipy import integrate, optimize, special, stats

from .distributions import Dist


# ---- asymptotic variances (of sqrt(n) * (estimator - limit)) ----------------------------------
def avar_mean(d: Dist) -> float:
    return d.var


def avar_median(d: Dist) -> float:
    return 1.0 / (4.0 * d.f_at_median ** 2)


def are_median_vs_mean(d: Dist) -> float:
    """ARE = Var(mean) / Var(median) = 4 sigma^2 f(m)^2; values > 1 favour the median."""
    return d.var * 4.0 * d.f_at_median ** 2 if np.isfinite(d.var) else float("inf")


def crlb_n(d: Dist) -> float:
    """n * Cramer-Rao lower bound for the location parameter (nan if not a location model)."""
    I = d.fisher_location
    return 1.0 / I if np.isfinite(I) else float("nan")


def asymptotic_mse(d: Dist, estimator: str, n: int) -> float:
    """First-order MSE of the mean or median at sample size n, about theta = 0."""
    if estimator == "mean":
        return avar_mean(d) / n + (d.mean_limit - d.theta) ** 2
    if estimator == "median":
        return avar_median(d) / n + (d.median_limit - d.theta) ** 2
    raise ValueError(estimator)


# ---- exact finite-n variance of the median at N(0,1) (odd n) ----------------------------------
def exact_median_nvar_normal(n: int) -> float:
    """Exact n * Var(M_n) for an odd-size N(0,1) sample, by integrating the order-statistic density
    g(x) = n!/(k!)^2 Phi(x)^k (1-Phi(x))^k phi(x),  k = (n-1)/2."""
    if n % 2 == 0:
        raise ValueError("closed form is for odd n (even n needs the joint law of two order statistics)")
    k = (n - 1) // 2
    logc = special.gammaln(n + 1) - 2 * special.gammaln(k + 1)
    g = lambda x: np.exp(logc + k * stats.norm.logcdf(x) + k * stats.norm.logsf(x)) * stats.norm.pdf(x)
    val, _ = integrate.quad(lambda x: x * x * g(x), -12, 12, points=[0.0], limit=400)
    return float(n * val)


# ---- Student t: efficiency curve and crossover ------------------------------------------------
def are_t(nu: float) -> float:
    """ARE(median vs mean) for Student t_nu (nu > 2): nu/(nu-2) * 4 f(0)^2."""
    f0 = np.exp(special.gammaln((nu + 1) / 2) - special.gammaln(nu / 2)) / np.sqrt(nu * np.pi)
    return nu / (nu - 2) * 4 * f0 ** 2


def t_crossover() -> float:
    """Degrees of freedom nu* at which the median and the mean are equally efficient."""
    return float(optimize.brentq(lambda v: are_t(v) - 1.0, 3.0, 10.0))


# ---- comparator estimators at the normal model (ARE vs. the mean) -----------------------------
def are_hl_normal() -> float:
    return 3.0 / np.pi


def are_trimmed_normal(alpha: float) -> float:
    """ARE of the alpha-trimmed mean vs. the mean at N(0,1):  1 / avar,
    avar = [ int_{-c}^{c} x^2 phi(x) dx + 2 alpha c^2 ] / (1 - 2 alpha)^2,  c = z_{1-alpha}."""
    c = stats.norm.ppf(1 - alpha)
    inner, _ = integrate.quad(lambda x: x * x * stats.norm.pdf(x), -c, c)
    return (1 - 2 * alpha) ** 2 / (inner + 2 * alpha * c * c)


def are_huber_normal(k: float) -> float:
    """ARE of the Huber M-estimator (known scale) vs. the mean at N(0,1): (E psi')^2 / E psi^2."""
    e_dpsi = 2 * stats.norm.cdf(k) - 1
    e_psi2, _ = integrate.quad(lambda x: np.clip(x, -k, k) ** 2 * stats.norm.pdf(x), -np.inf, np.inf)
    return e_dpsi ** 2 / e_psi2


def normal_n_var(estimator: str, trim: float = 0.10, k: float = 1.345) -> float:
    """n * Var of each estimator at N(0,1), asymptotically."""
    if estimator == "mean":
        return 1.0
    if estimator == "median":
        return np.pi / 2
    if estimator == "hl":
        return np.pi / 3
    if estimator == "trim10":
        return 1.0 / are_trimmed_normal(trim)
    if estimator == "huber":
        return 1.0 / are_huber_normal(k)
    raise ValueError(estimator)


# ---- influence functions at N(0,1) (limit of the sensitivity curves) ---------------------------
def influence_normal(estimator: str, y, trim: float = 0.10, k: float = 1.345):
    y = np.asarray(y, dtype=float)
    if estimator == "mean":
        return y
    if estimator == "median":
        return np.sign(y) / (2 * stats.norm.pdf(0.0))
    if estimator == "hl":
        return (stats.norm.cdf(y) - 0.5) / (1.0 / (2.0 * np.sqrt(np.pi)))
    if estimator == "trim10":
        c = stats.norm.ppf(1 - trim)
        return np.clip(y, -c, c) / (1 - 2 * trim)
    if estimator == "huber":
        return np.clip(y, -k, k) / (2 * stats.norm.cdf(k) - 1)
    raise ValueError(estimator)
