import numpy as np
import pytest

from mvm.estimators import (ESTIMATORS, hodges_lehmann, huber_location, mean_est, median_est,
                            running_median, trimmed_mean)


def test_basic_values():
    X = np.array([[1.0, 2.0, 3.0, 4.0, 100.0]])
    assert mean_est(X)[0] == pytest.approx(22.0)
    assert median_est(X)[0] == 3.0
    assert trimmed_mean(X, 0.2)[0] == pytest.approx(3.0)            # drops 1 and 100 -> mean(2,3,4)
    assert hodges_lehmann(np.array([[1.0, 2.0, 3.0]]))[0] == 2.0    # Walsh averages incl. i = j


def test_hodges_lehmann_walsh_by_brute_force():
    x = np.random.default_rng(0).normal(size=9)
    walsh = [(x[i] + x[j]) / 2 for i in range(9) for j in range(i, 9)]
    assert hodges_lehmann(x[None, :])[0] == pytest.approx(np.median(walsh))


@pytest.mark.parametrize("name", list(ESTIMATORS))
def test_location_and_scale_equivariance(name):
    f = ESTIMATORS[name]
    X = np.random.default_rng(1).standard_normal((5, 31))
    assert np.allclose(f(X + 7.5), f(X) + 7.5, atol=1e-7)
    assert np.allclose(f(3.0 * X), 3.0 * f(X), atol=1e-7)


@pytest.mark.parametrize("name", list(ESTIMATORS))
def test_symmetric_sample_gives_centre(name):
    x = np.r_[-np.arange(1, 21), np.arange(1, 21), 0.0]
    assert ESTIMATORS[name](x[None, :])[0] == pytest.approx(0.0, abs=1e-6)


def test_huber_is_between_mean_and_median_on_skewed_data():
    x = np.r_[np.random.default_rng(3).normal(size=50), 6.0, 7.0, 8.0][None, :]
    lo, hi = np.median(x), x.mean()
    assert lo <= huber_location(x)[0] <= hi


def test_running_median_matches_numpy():
    x = np.random.default_rng(4).standard_cauchy(300)
    rm = running_median(x)
    assert all(np.isclose(rm[t], np.median(x[:t + 1])) for t in range(300))
