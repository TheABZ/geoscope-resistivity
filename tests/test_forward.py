"""Physics benchmarks independent of layered-earth numerical inversion."""

import numpy as np
import pytest

from geoscope.forward import LayeredEarth, simulate_wenner, synthetic_observations
from geoscope.demo import DEFAULT_CONFIG, run


def test_uniform_halfspace():
    """A uniform half-space must return its resistivity at every spacing."""
    a = np.geomspace(0.1, 1000, 20)
    np.testing.assert_allclose(simulate_wenner(LayeredEarth((123.0,), ()), a), 123, rtol=1e-5)


def test_invisible_interfaces():
    earth = LayeredEarth((100, 100, 100), (5, 20))
    np.testing.assert_allclose(simulate_wenner(earth, np.geomspace(1, 200, 30)), 100, rtol=1e-5)


@pytest.mark.parametrize("rho2", [25.0, 800.0])
def test_two_layer_against_image_series(rho2):
    """Compare Hankel solver with the independent two-layer image expansion.

    rho_a/rho1 = 1 + 4 sum k^n [a/sqrt(a²+(2nh)²)
                                - a/sqrt(4a²+(2nh)²)]
    k = (rho2-rho1)/(rho2+rho1). Both conductive and resistive cases tested.
    """
    rho1, h = 100.0, 5.0
    a = np.geomspace(0.2, 500, 35)
    n = np.arange(1, 1001)[:, None]
    k = (rho2 - rho1) / (rho2 + rho1)
    expected = rho1 * (1 + 4 * np.sum(k**n * (
        a / np.sqrt(a**2 + (2*n*h)**2)
        - a / np.sqrt(4*a**2 + (2*n*h)**2)), axis=0))
    actual = simulate_wenner(LayeredEarth((rho1, rho2), (h,)), a)
    np.testing.assert_allclose(actual, expected, rtol=2e-4)


def test_resistivity_scaling():
    a = np.geomspace(1, 200, 15)
    base = simulate_wenner(LayeredEarth((100, 25, 800), (5, 20)), a)
    scaled = simulate_wenner(LayeredEarth((300, 75, 2400), (5, 20)), a)
    np.testing.assert_allclose(scaled, 3 * base, rtol=1e-7)


@pytest.mark.parametrize("rho,h", [([], []), ([0], []), ([np.nan], []), ([100, 50], []), ([100, 50], [-1]), ([[100]], [])])
def test_invalid_earth(rho, h):
    with pytest.raises(ValueError):
        LayeredEarth(rho, h)


@pytest.mark.parametrize("a", [[], [0], [-2], [np.nan], [[1, 2]]])
def test_invalid_spacing(a):
    with pytest.raises(ValueError):
        simulate_wenner(LayeredEarth((100,), ()), a)


def test_seed_and_noise_definition():
    clean = np.full(10000, 100.0)
    first, sigma = synthetic_observations(clean, 0.03, seed=7)
    second, _ = synthetic_observations(clean, 0.03, seed=7)
    np.testing.assert_array_equal(first, second)
    assert np.all(first > 0)
    np.testing.assert_allclose(sigma, 0.03)
    assert abs(np.std(np.log(first / clean)) - 0.03) < 0.001
    np.testing.assert_array_equal(synthetic_observations(clean, 0)[0], clean)
    with pytest.raises(ValueError):
        synthetic_observations(clean, -0.1)


def test_full_experiment(tmp_path):
    import json
    run(DEFAULT_CONFIG, tmp_path)
    table = np.loadtxt(tmp_path / "synthetic_sounding.csv", delimiter=",", skiprows=1)
    assert table.shape == (40, 4)
    assert np.all(table > 0)
    assert (tmp_path / "first_sounding.svg").stat().st_size > 1000
    metadata = json.loads((tmp_path / "experiment.json").read_text())
    assert metadata["config"]["seed"] == 42
    assert "synthetic" in metadata["data_origin"]
