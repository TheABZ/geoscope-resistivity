"""GeoScope: reproducible electrical resistivity experiments."""

from .forward import LayeredEarth, simulate_wenner

__all__ = ["LayeredEarth", "simulate_wenner"]
__version__ = "0.2.0"
