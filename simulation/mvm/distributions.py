"""Sampling distributions with their exact theoretical ingredients.

Every distribution is centred at theta = 0, the quantity we want to estimate.  For the
asymmetric 'shift' contamination the *target* is still the clean centre 0, so both the
mean and the median are biased there and the comparison is on MSE.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import integrate, optimize, special

SQRT2PI = np.sqrt(2.0 * np.pi)


def _phi(x):
    return np.exp(-0.5 * np.square(x)) / SQRT2PI


@dataclass(frozen=True)
class Dist:
    kind: str                 # 'normal' | 't' | 'laplace' | 'tukey' | 'shift'
    nu: float = 0.0           # t degrees of freedom
    eps: float = 0.0          # contamination fraction
    scale: float = 3.0        # 'tukey': sd of contaminating normal
    shift: float = 10.0       # 'shift': mean of contaminating normal (unit variance)

    # ---- labels --------------------------------------------------------------------------
    @property
    def name(self) -> str:
        if self.kind == "normal":
            return "N(0,1)"
        if self.kind == "t":
            return "Cauchy" if self.nu == 1 else f"t{self.nu:g}"
        if self.kind == "laplace":
            return "Laplace"
        if self.kind == "tukey":
            return f"Tukey eps={self.eps:g}"
        if self.kind == "shift":
            return f"Shift eps={self.eps:g}"
        raise ValueError(self.kind)

    # ---- sampling ------------------------------------------------------------------------
    def sample(self, rng: np.random.Generator, size) -> np.ndarray:
        if self.kind == "normal":
            return rng.standard_normal(size)
        if self.kind == "t":
            return rng.standard_t(self.nu, size)
        if self.kind == "laplace":
            return rng.laplace(0.0, 1.0, size)
        z = rng.standard_normal(size)
        if self.eps == 0.0:
            return z
        mask = rng.random(size) < self.eps
        if self.kind == "tukey":
            return np.where(mask, self.scale * z, z)
        if self.kind == "shift":
            return z + self.shift * mask
        raise ValueError(self.kind)

    # ---- exact density / cdf -------------------------------------------------------------
    def pdf(self, x):
        x = np.asarray(x, dtype=float)
        if self.kind == "normal":
            return _phi(x)
        if self.kind == "t":
            nu = self.nu
            c = np.exp(special.gammaln((nu + 1) / 2) - special.gammaln(nu / 2)) / np.sqrt(nu * np.pi)
            return c * (1 + x**2 / nu) ** (-(nu + 1) / 2)
        if self.kind == "laplace":
            return 0.5 * np.exp(-np.abs(x))
        if self.kind == "tukey":
            return (1 - self.eps) * _phi(x) + self.eps * _phi(x / self.scale) / self.scale
        if self.kind == "shift":
            return (1 - self.eps) * _phi(x) + self.eps * _phi(x - self.shift)
        raise ValueError(self.kind)

    def cdf(self, x):
        x = np.asarray(x, dtype=float)
        if self.kind == "normal":
            return special.ndtr(x)
        if self.kind == "t":
            return special.stdtr(self.nu, x)
        if self.kind == "laplace":
            return np.where(x < 0, 0.5 * np.exp(x), 1 - 0.5 * np.exp(-x))
        if self.kind == "tukey":
            return (1 - self.eps) * special.ndtr(x) + self.eps * special.ndtr(x / self.scale)
        if self.kind == "shift":
            return (1 - self.eps) * special.ndtr(x) + self.eps * special.ndtr(x - self.shift)
        raise ValueError(self.kind)

    def _dpdf(self, x):
        """Derivative of the density (only needed for the Fisher information of the mixtures)."""
        if self.kind == "tukey":
            s = self.scale
            return -(1 - self.eps) * x * _phi(x) - self.eps * (x / s**2) * _phi(x / s) / s
        raise NotImplementedError

    # ---- theoretical ingredients ---------------------------------------------------------
    @property
    def theta(self) -> float:
        """The parameter we aim to estimate: the centre of the (clean) distribution."""
        return 0.0

    @property
    def symmetric(self) -> bool:
        return self.kind != "shift"

    @property
    def var(self) -> float:
        """Variance sigma^2 (inf when it does not exist)."""
        if self.kind == "normal":
            return 1.0
        if self.kind == "t":
            return self.nu / (self.nu - 2) if self.nu > 2 else float("inf")
        if self.kind == "laplace":
            return 2.0
        if self.kind == "tukey":
            return 1.0 + (self.scale**2 - 1.0) * self.eps
        if self.kind == "shift":
            return 1.0 + self.eps * (1 - self.eps) * self.shift**2
        raise ValueError(self.kind)

    @property
    def mean_limit(self) -> float:
        """Value the sample mean converges to (when it converges)."""
        if self.kind == "t" and self.nu <= 1:
            return float("nan")
        return self.eps * self.shift if self.kind == "shift" else 0.0

    @property
    def median_limit(self) -> float:
        """Value the sample median converges to: the population median."""
        if self.kind == "shift" and self.eps > 0:
            return float(optimize.brentq(lambda m: self.cdf(m) - 0.5, -1.0, self.shift))
        return 0.0

    @property
    def f_at_median(self) -> float:
        return float(self.pdf(self.median_limit))

    @property
    def fisher_location(self) -> float:
        """Fisher information for the location parameter (nan if not a location family)."""
        if self.kind == "normal":
            return 1.0
        if self.kind == "t":
            return (self.nu + 1) / (self.nu + 3)
        if self.kind == "laplace":
            return 1.0
        if self.kind == "tukey":
            if self.eps == 0.0:
                return 1.0
            f = lambda x: self._dpdf(x) ** 2 / self.pdf(x)
            val, _ = integrate.quad(f, -30, 30, points=[0.0], limit=400)
            return float(val)
        return float("nan")
