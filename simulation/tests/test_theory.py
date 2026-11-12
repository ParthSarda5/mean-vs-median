import numpy as np
import pytest

from mvm import theory as T
from mvm.distributions import Dist


def test_normal_are_is_two_over_pi():
    assert T.are_median_vs_mean(Dist("normal")) == pytest.approx(2 / np.pi)


def test_t_are_values_and_crossover():
    assert T.are_t(3) == pytest.approx(1.6211, abs=1e-3)
    assert T.are_t(5) == pytest.approx(0.9607, abs=1e-3)
    assert T.t_crossover() == pytest.approx(4.679, abs=1e-3)
    assert T.are_t(3) == pytest.approx(T.are_median_vs_mean(Dist("t", nu=3)))


def test_laplace_median_is_mle():
    d = Dist("laplace")
    assert T.are_median_vs_mean(d) == pytest.approx(2.0)
    assert T.avar_median(d) == pytest.approx(T.crlb_n(d))          # median attains the CRLB


def test_cauchy_median_efficiency_and_crlb():
    d = Dist("t", nu=1)
    assert T.crlb_n(d) == pytest.approx(2.0)
    assert T.avar_median(d) == pytest.approx(np.pi**2 / 4)
    assert T.crlb_n(d) / T.avar_median(d) == pytest.approx(8 / np.pi**2)
    assert np.isinf(T.are_median_vs_mean(d))


def test_normal_mean_attains_crlb():
    d = Dist("normal")
    assert T.avar_mean(d) == pytest.approx(T.crlb_n(d))


def test_fisher_information_of_t_matches_numerical_integral():
    from scipy import integrate
    d = Dist("t", nu=5)
    h = 1e-5
    score2 = lambda x: ((np.log(d.pdf(x + h)) - np.log(d.pdf(x - h))) / (2 * h)) ** 2 * d.pdf(x)
    val, _ = integrate.quad(score2, -200, 200, limit=400)
    assert val == pytest.approx(d.fisher_location, rel=1e-4)


def test_tukey_fisher_exceeds_one_over_variance_bound():
    d = Dist("tukey", eps=0.1, scale=3.0)
    assert d.fisher_location >= 1.0 / d.var          # Cramer-Rao: the mean cannot beat 1/I


def test_comparator_efficiencies_at_normal():
    assert T.are_hl_normal() == pytest.approx(0.9549, abs=1e-3)
    assert T.are_trimmed_normal(0.10) == pytest.approx(0.943, abs=1e-3)
    assert T.are_huber_normal(1.345) == pytest.approx(0.95, abs=1e-3)


def test_shift_contamination_limits():
    d = Dist("shift", eps=0.10, shift=10.0)
    assert d.mean_limit == pytest.approx(1.0)
    assert d.cdf(d.median_limit) == pytest.approx(0.5)
    assert 0 < d.median_limit < 0.5                   # median moves far less than the mean
    assert d.var == pytest.approx(1 + 0.1 * 0.9 * 100)


def test_variance_of_mixtures():
    assert Dist("tukey", eps=0.1, scale=3.0).var == pytest.approx(1 + 8 * 0.1)
    assert np.isinf(Dist("t", nu=2).var) and np.isinf(Dist("t", nu=1).var)


def test_influence_functions_bounded_vs_unbounded():
    y = np.array([-100.0, 100.0])
    assert np.all(np.abs(T.influence_normal("mean", y)) == 100.0)
    for est in ("median", "huber", "trim10", "hl"):
        assert np.all(np.abs(T.influence_normal(est, y)) < 3.0)
    assert np.abs(T.influence_normal("median", 5.0)) == pytest.approx(1.2533, abs=1e-3)


def test_exact_median_variance_odd_n_and_the_even_odd_effect():
    """n*Var(median) at N(0,1): exact for odd n, simulated for even n.  Even n averages two central
    order statistics and is markedly lower than both odd neighbours (Hodges-Lehmann 1967)."""
    from mvm.engine import simulate_cell
    exact = {n: T.exact_median_nvar_normal(n) for n in (9, 11, 25)}
    assert exact[9] == pytest.approx(1.4949, abs=1e-3)
    assert exact[11] == pytest.approx(1.5088, abs=1e-3)
    assert exact[9] < exact[11] < exact[25] < np.pi / 2                       # increases towards pi/2
    r11 = next(r for r in simulate_cell(Dist("normal"), 11, 100_000, ("median",)))
    assert r11["n_var"] == pytest.approx(exact[11], abs=4 * 11 * r11["se_var"])    # simulation agrees with exact
    r10 = next(r for r in simulate_cell(Dist("normal"), 10, 100_000, ("median",)))
    assert r10["n_var"] < exact[9] - 0.05 and r10["n_var"] < exact[11] - 0.05  # the even-n dip
