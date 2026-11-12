"""Surface Wenner-alpha soundings over horizontal, isotropic layers.

Coordinates are metres; resistivity is ohm metres. The final layer is an
infinite half-space. SimPEG supplies the layered-earth Hankel-transform solver.
"""

from dataclasses import dataclass

import numpy as np
from simpeg import maps
from simpeg.electromagnetics.static import resistivity as dc


@dataclass(frozen=True)
class LayeredEarth:
    """N resistivities and N-1 finite layer thicknesses, top to bottom."""

    resistivities: tuple[float, ...]
    thicknesses: tuple[float, ...]

    def __post_init__(self):
        for name in ("resistivities", "thicknesses"):
            values = np.asarray(getattr(self, name), dtype=float)
            if values.ndim != 1:
                raise ValueError(f"{name} must be one-dimensional")
            if not np.all(np.isfinite(values)) or np.any(values <= 0):
                raise ValueError(f"{name} must contain finite, positive values")
            object.__setattr__(self, name, tuple(values.tolist()))
        if not self.resistivities:
            raise ValueError("At least one resistivity is required")
        if len(self.thicknesses) != len(self.resistivities) - 1:
            raise ValueError("N resistivities require N-1 thicknesses")


def simulate_wenner(earth: LayeredEarth, spacings) -> np.ndarray:
    """Predict apparent resistivity for Wenner-alpha electrode spacings a.

    The electrode order is A, M, N, B at x = -1.5a, -0.5a, 0.5a, 1.5a.
    Each measurement uses a 1 A source and has geometric factor K = 2*pi*a.
    Larger spacings change depth sensitivity; a is not a depth coordinate.
    """
    spacings = np.asarray(spacings, dtype=float)
    if spacings.ndim != 1 or spacings.size == 0:
        raise ValueError("spacings must be a nonempty one-dimensional sequence")
    if not np.all(np.isfinite(spacings)) or np.any(spacings <= 0):
        raise ValueError("spacings must be finite and positive")
    sources = []
    for a in spacings:
        receiver = dc.receivers.Dipole(
            locations_m=np.array([[-0.5 * a, 0.0, 0.0]]),
            locations_n=np.array([[0.5 * a, 0.0, 0.0]]),
            data_type="apparent_resistivity",
        )
        sources.append(dc.sources.Dipole(
            receiver_list=[receiver],
            location_a=np.array([-1.5 * a, 0.0, 0.0]),
            location_b=np.array([1.5 * a, 0.0, 0.0]),
            current=1.0,
        ))
    simulation = dc.Simulation1DLayers(
        survey=dc.Survey(sources),
        rho=np.asarray(earth.resistivities),
        thicknesses=np.asarray(earth.thicknesses),
    )
    prediction = np.asarray(simulation.dpred()).reshape(-1)
    if not np.all(np.isfinite(prediction)) or np.any(prediction <= 0):
        raise RuntimeError("Forward solver returned invalid apparent resistivities")
    return prediction


class WennerOperator:
    """Reuse survey geometry while changing log resistivities and thicknesses."""

    def __init__(self, spacings, n_layers):
        if isinstance(n_layers, bool) or not isinstance(n_layers, int) or n_layers < 1:
            raise ValueError("n_layers must be a positive integer")
        a = np.asarray(spacings, dtype=float)
        if a.ndim != 1 or not a.size or not np.all(np.isfinite(a)) or np.any(a <= 0):
            raise ValueError("spacings must be a nonempty positive finite 1D sequence")
        sources = []
        for spacing in a:
            receiver = dc.receivers.Dipole(
                np.array([[-0.5 * spacing, 0, 0]]),
                np.array([[0.5 * spacing, 0, 0]]),
                data_type="apparent_resistivity",
            )
            sources.append(dc.sources.Dipole(
                [receiver], location_a=np.array([-1.5 * spacing, 0, 0]),
                location_b=np.array([1.5 * spacing, 0, 0]), current=1,
            ))
        self.n_layers = n_layers
        self.n_params = 2 * n_layers - 1
        if n_layers > 1:
            wires = maps.Wires(("rho", n_layers), ("h", n_layers - 1))
            kwargs = {"rhoMap": maps.ExpMap(nP=n_layers) * wires.rho}
            kwargs["thicknessesMap"] = maps.ExpMap(nP=n_layers - 1) * wires.h
        else:
            kwargs = {"rhoMap": maps.ExpMap(nP=1), "thicknesses": np.array([])}
        self.simulation = dc.Simulation1DLayers(survey=dc.Survey(sources), **kwargs)

    def __call__(self, log_parameters):
        m = np.asarray(log_parameters, dtype=float)
        if m.shape != (self.n_params,) or not np.all(np.isfinite(m)):
            raise ValueError("Invalid log-parameter vector")
        predicted = np.asarray(self.simulation.dpred(m)).reshape(-1)
        if not np.all(np.isfinite(predicted)) or np.any(predicted <= 0):
            raise ValueError("Invalid numerical forward response")
        return predicted


def synthetic_observations(clean, relative_noise=0.03, seed=42):
    """Gaussian noise in log space: log(observed/clean) ~ N(0, sigma_log).

    relative_noise is sigma_log, approximately fractional error when small.
    The resulting positive observations have median clean (not mean clean).
    Returns (observed, sigma_log), with one standard deviation per datum.
    """
    clean = np.asarray(clean, dtype=float)
    if clean.ndim != 1 or not clean.size or not np.all(np.isfinite(clean)) or np.any(clean <= 0):
        raise ValueError("clean must be a nonempty positive finite 1D sequence")
    if not np.isfinite(relative_noise) or relative_noise < 0:
        raise ValueError("relative_noise must be finite and nonnegative")
    rng = np.random.default_rng(seed)
    observed = clean * np.exp(rng.normal(0, relative_noise, clean.size))
    if not np.all(np.isfinite(observed)) or np.any(observed <= 0):
        raise ValueError("Noise level produced nonfinite or zero observations")
    return observed, np.full(clean.shape, relative_noise)
