"""Bounded log-space inversion and explicit data uncertainty validation."""

from dataclasses import dataclass
import numpy as np
from scipy.optimize import least_squares

from .forward import WennerOperator


@dataclass
class Sounding:
    spacings: np.ndarray
    observed: np.ndarray
    sigma_log: np.ndarray

    def __post_init__(self):
        for key in ("spacings", "observed", "sigma_log"):
            value = np.array(getattr(self, key), dtype=float, copy=True)
            if value.ndim != 1 or not value.size or not np.all(np.isfinite(value)) or np.any(value <= 0):
                raise ValueError(f"{key} must be a positive finite nonempty 1D array")
            setattr(self, key, value)
        if not (self.spacings.shape == self.observed.shape == self.sigma_log.shape):
            raise ValueError("All sounding arrays must have identical lengths")


class InverseProblem:
    """Known sigma_log, fixed layer count, uniform priors in bounded log space."""

    def __init__(self, sounding, n_layers=3, rho_bounds=(10., 3000.), h_bounds=(0.5, 100.)):
        self.data = sounding
        self.operator = WennerOperator(sounding.spacings, n_layers)
        self.names = [f"rho_{i+1}" for i in range(n_layers)] + [f"h_{i+1}" for i in range(n_layers-1)]
        for bounds in (rho_bounds, h_bounds):
            if len(bounds) != 2 or not np.all(np.isfinite(bounds)) or not 0 < bounds[0] < bounds[1]:
                raise ValueError("Bounds must be finite, positive, and increasing")
        self.lower = np.log([rho_bounds[0]] * n_layers + [h_bounds[0]] * (n_layers-1))
        self.upper = np.log([rho_bounds[1]] * n_layers + [h_bounds[1]] * (n_layers-1))
        if sounding.observed.size < len(self.names):
            raise ValueError("Need at least as many measurements as model parameters")

    def residual(self, m):
        return (np.log(self.operator(m)) - np.log(self.data.observed)) / self.data.sigma_log

    def log_probability(self, m):
        if np.any(m < self.lower) or np.any(m > self.upper) or not np.all(np.isfinite(m)):
            return -np.inf
        try:
            residual = self.residual(m)
        except (ValueError, FloatingPointError):
            return -np.inf
        return -0.5 * np.dot(residual, residual)

    def fit(self, starts=8, seed=42):
        if isinstance(starts, bool) or not isinstance(starts, int) or starts < 1:
            raise ValueError("starts must be a positive integer")
        rng = np.random.default_rng(seed)
        fits = []
        for i in range(starts):
            m0 = (self.lower + self.upper) / 2 if i == 0 else rng.uniform(self.lower, self.upper)
            result = least_squares(self.residual, m0, bounds=(self.lower, self.upper),
                                   max_nfev=500, ftol=1e-9, xtol=1e-9, gtol=1e-9)
            fits.append(result)
        successful = [r for r in fits if r.success and np.isfinite(r.cost)]
        if not successful:
            raise RuntimeError("No inversion start converged; inspect data, bounds, and geometry")
        best = min(successful, key=lambda r: r.cost)
        return best, fits
