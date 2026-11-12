import numpy as np
import pytest

from geoscope.forward import LayeredEarth, WennerOperator, simulate_wenner
from geoscope.inverse import Sounding, InverseProblem
from geoscope.posterior import sample_posterior, diagnostics
from geoscope.analysis import read_sounding


def test_reused_operator_updates_all_parameters():
    a = np.geomspace(1, 200, 30)
    operator = WennerOperator(a, 3)
    for values in ([100, 25, 800, 5, 20], [150, 40, 600, 10, 15], [100, 25, 800, 5, 20]):
        expected = simulate_wenner(LayeredEarth(values[:3], values[3:]), a)
        np.testing.assert_allclose(operator(np.log(values)), expected, rtol=1e-8)


def test_recover_known_three_layer_model():
    a = np.geomspace(1, 200, 40)
    true = np.array([100, 25, 800, 5, 20])
    data = simulate_wenner(LayeredEarth(true[:3], true[3:]), a)
    problem = InverseProblem(Sounding(a, data, np.full(a.size, .03)))
    best, fits = problem.fit(starts=5)
    np.testing.assert_allclose(np.exp(best.x), true, rtol=.005)
    assert len(fits) == 5
    assert np.sqrt(np.mean(best.fun**2)) < 1e-4
    assert problem.log_probability(problem.lower - 1) == -np.inf


def test_sampler_against_analytic_halfspace_posterior():
    a = np.geomspace(1, 200, 20)
    sigma = .1
    data = Sounding(a, np.full(20, 100.), np.full(20, sigma))
    problem = InverseProblem(data, n_layers=1)
    fit, _ = problem.fit(starts=2)
    samples, report = sample_posterior(problem, fit, draws=4000, warmup=1000)
    # Uniform log prior and Gaussian log errors yield a Normal posterior for
    # log(rho), with mean log(100), SD sigma/sqrt(N); bounds are remote here.
    assert abs(samples.mean() - np.log(100)) < .004
    assert .85 < samples.std() / (sigma / np.sqrt(20)) < 1.15
    assert report["status"] == "diagnostics_passed"
    assert all(.05 < x < .8 for x in report["acceptance_rates"])


def test_diagnostics_detect_separated_chains():
    rng = np.random.default_rng(2)
    bad = rng.normal(size=(4, 500, 1))
    bad[2:] += 5
    assert diagnostics(bad, ["rho_1"])["status"] == "inconclusive_sampling"


@pytest.mark.parametrize("a,y,sigma", [([1], [100], [0]), ([1,2], [100], [.03]), ([1], [np.nan], [.03])])
def test_invalid_measurement_errors(a, y, sigma):
    with pytest.raises(ValueError):
        Sounding(a, y, sigma)


def test_csv_requires_correct_geometry_fields(tmp_path):
    path = tmp_path / "bad.csv"
    path.write_text("depth,rho\n1,100\n")
    with pytest.raises(ValueError, match="CSV must contain"):
        read_sounding(path)


def test_sampler_repeatability():
    problem = InverseProblem(Sounding([1,2,3], [100,100,100], [.1,.1,.1]), 1)
    fit, _ = problem.fit(starts=1)
    first, _ = sample_posterior(problem, fit, draws=30, warmup=50, seed=5)
    second, _ = sample_posterior(problem, fit, draws=30, warmup=50, seed=5)
    np.testing.assert_array_equal(first, second)

