"""Statistical sanity checks of the Monte Carlo output (fixed seeds => deterministic)."""
import numpy as np
import pytest

from mvm import config as C
from mvm import theory as T
from mvm.distributions import Dist
from mvm.engine import simulate_cell
from mvm.estimators import ESTIMATORS


def _row(rows, est):
    return next(r for r in rows if r["estimator"] == est)


def test_samplers_have_the_right_moments():
    rng = np.random.default_rng(0)
    for d in (Dist("normal"), Dist("laplace"), Dist("t", nu=10), Dist("tukey", eps=0.2, scale=3.0)):
        x = d.sample(rng, 400_000)
        assert x.var() == pytest.approx(d.var, rel=0.03)
    x = Dist("shift", eps=0.1, shift=10.0).sample(rng, 400_000)
    assert x.mean() == pytest.approx(1.0, abs=0.02)
    assert np.median(x) == pytest.approx(Dist("shift", eps=0.1).median_limit, abs=0.01)


def test_cdf_pdf_consistent():
    for d in (Dist("normal"), Dist("t", nu=3), Dist("laplace"), Dist("tukey", eps=0.1), Dist("shift", eps=0.1)):
        xs = np.linspace(-5, 5, 11)
        h = 1e-5
        num = (d.cdf(xs + h) - d.cdf(xs - h)) / (2 * h)
        assert np.allclose(num, d.pdf(xs), atol=1e-6)


def test_normal_variances_match_theory():
    rows = simulate_cell(Dist("normal"), 400, 40_000, ("mean", "median"))
    m, md = _row(rows, "mean"), _row(rows, "median")
    assert m["n_var"] == pytest.approx(1.0, abs=4 * 400 * m["se_var"])
    assert md["n_var"] == pytest.approx(np.pi / 2, abs=4 * 400 * md["se_var"] + 0.02)
    assert md["var_ratio_vs_mean"] == pytest.approx(2 / np.pi, abs=0.02)


def test_t3_median_beats_mean_t10_mean_beats_median():
    for nu, median_wins in ((3, True), (10, False)):
        rows = simulate_cell(Dist("t", nu=nu), 200, 40_000, ("mean", "median"))
        ratio = _row(rows, "median")["var_ratio_vs_mean"]      # Var(mean)/Var(median)
        assert (ratio > 1) == median_wins
        assert ratio == pytest.approx(T.are_t(nu), rel=0.06)


def test_cauchy_mean_does_not_improve_with_n_but_median_does():
    q = {}
    for n in (25, 400):
        rows = simulate_cell(Dist("t", nu=1), n, 20_000, ("mean", "median"))
        q[n] = (_row(rows, "mean")["q50_abs"], _row(rows, "median")["q50_abs"])
    assert q[400][0] == pytest.approx(q[25][0], rel=0.15)       # |mean| has the same law for every n
    assert q[400][1] < 0.4 * q[25][1]                           # median error shrinks ~ n^{-1/2}


def test_shift_contamination_bias_matches_theory():
    d = Dist("shift", eps=0.1, shift=10.0)
    rows = simulate_cell(d, 200, 20_000, ("mean", "median"))
    assert _row(rows, "mean")["bias"] == pytest.approx(d.mean_limit, abs=0.05)
    assert _row(rows, "median")["bias"] == pytest.approx(d.median_limit, abs=0.03)


def test_engine_is_reproducible_and_order_independent():
    a = simulate_cell(Dist("normal"), 50, 2000, ("mean", "median"))
    b = simulate_cell(Dist("normal"), 50, 2000, ("mean", "median"))
    assert a == b


def test_breakdown_mean_at_one_outlier_median_at_half():
    n = 101
    X0 = np.random.default_rng(5).standard_normal((200, n))
    def med_abs(f, k):
        X = X0.copy(); X[:, :k] = C.BREAKDOWN_VALUE
        return np.median(np.abs(f(X)))
    assert med_abs(ESTIMATORS["mean"], 1) > 1e3
    assert med_abs(ESTIMATORS["median"], 50) < 5
    assert med_abs(ESTIMATORS["median"], 51) > 1e5
    assert med_abs(ESTIMATORS["trim10"], 10) < 5 and med_abs(ESTIMATORS["trim10"], 11) > 1e3
    assert med_abs(ESTIMATORS["hl"], 25) < 5 and med_abs(ESTIMATORS["hl"], 35) > 1e3   # breakdown ~ 29.3%


def test_sensitivity_curves_approach_influence_functions():
    from mvm.experiments import exp_sensitivity
    import pandas as pd
    df = exp_sensitivity(save=False)          # do not touch results/ from a test
    for est in ("mean", "median", "huber", "trim10", "hl"):
        for y in (-5.0, 5.0):
            sc = float(np.interp(y, df["y"], df[est]))
            assert sc == pytest.approx(float(T.influence_normal(est, y)), abs=0.35), (est, y, sc)
