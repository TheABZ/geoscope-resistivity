"""Run the first GeoScope experiment: geoscope --output results."""

import argparse
import json
from importlib.metadata import version
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .forward import LayeredEarth, simulate_wenner, synthetic_observations

DEFAULT_CONFIG = {
    "resistivities_ohm_m": [100.0, 25.0, 800.0],
    "thicknesses_m": [5.0, 20.0],
    "spacing_min_m": 1.0,
    "spacing_max_m": 200.0,
    "spacing_count": 40,
    "sigma_log": 0.03,
    "seed": 42,
}


def run(config, output):
    """Write a figure, CSV, and metadata for an explicit experiment config."""
    earth = LayeredEarth(config["resistivities_ohm_m"], config["thicknesses_m"])
    count = config["spacing_count"]
    lo, hi = config["spacing_min_m"], config["spacing_max_m"]
    if isinstance(count, bool) or not isinstance(count, int) or count < 2:
        raise ValueError("spacing_count must be an integer of at least 2")
    if not np.isfinite(lo) or not np.isfinite(hi) or not 0 < lo < hi:
        raise ValueError("Require 0 < spacing_min_m < spacing_max_m, both finite")
    spacings = np.geomspace(lo, hi, count)
    clean = simulate_wenner(earth, spacings)
    observed, sigma = synthetic_observations(clean, config["sigma_log"], config["seed"])
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    np.savetxt(
        output / "synthetic_sounding.csv",
        np.column_stack([spacings, clean, observed, sigma]), delimiter=",",
        header="a_m,rhoa_clean_ohm_m,rhoa_observed_ohm_m,sigma_log", comments="",
    )
    metadata = {
        "data_origin": "synthetic; not Eruwa field measurements",
        "array": "Wenner-alpha; A M N B; K=2*pi*a",
        "noise_model": "log(observed/clean) ~ Normal(0, sigma_log**2)",
        "config": config,
        "versions": {name: version(name) for name in ["geoscope-resistivity", "simpeg", "numpy", "matplotlib"]},
    }
    (output / "experiment.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    plot_experiment(earth, spacings, clean, observed, sigma, output / "first_sounding.svg")
    return clean


def plot_experiment(earth, spacings, clean, observed, sigma, target):
    """Plot true layer properties and apparent resistivity on separate axes."""
    with plt.rc_context({"font.family": "DejaVu Sans", "font.size": 11, "svg.fonttype": "none"}):
        fig, (left, right) = plt.subplots(1, 2, figsize=(12, 5.8), gridspec_kw={"width_ratios": [1, 1.65]})
        fig.suptitle("GeoScope | Beneath the sounding curve", fontsize=19, fontweight="bold", x=0.06, ha="left")
        boundaries = np.r_[0.0, np.cumsum(earth.thicknesses)]
        bottom = max(40.0, 1.6 * boundaries[-1])
        depths = np.r_[boundaries, bottom]
        colors = ["#0b7285", "#e9a23b", "#485696"]
        for i, rho in enumerate(earth.resistivities):
            left.axhspan(depths[i], depths[i + 1], color=colors[i % len(colors)], alpha=0.09)
            left.plot([rho, rho], depths[i:i + 2], color=colors[i % len(colors)], lw=3)
            if i + 1 < len(earth.resistivities):
                left.plot([rho, earth.resistivities[i + 1]], [depths[i + 1]] * 2, color="#56616b", lw=1.3)
            left.text(0.04, (depths[i] + depths[i + 1]) / 2, f"Layer {i + 1}: {rho:g} Ω m", transform=left.get_yaxis_transform(), fontsize=10)
        left.set(xscale="log", ylim=(bottom, 0), xlabel="True resistivity (Ω m)", ylabel="Depth (m)", title="Known synthetic Earth")
        left.set_xlim(min(earth.resistivities) / 3, max(earth.resistivities) * 3)
        left.text(0.03, 0.02, "Bottom layer continues to infinite depth", transform=left.transAxes, fontsize=8)
        right.loglog(spacings, clean, color="#0b7285", lw=2.5, label="Noise-free forward prediction")
        error = np.vstack([observed - observed * np.exp(-sigma), observed * np.exp(sigma) - observed])
        right.errorbar(spacings, observed, yerr=error, fmt="o", ms=3.5, color="#b66a13", alpha=0.8, capsize=2, label="Synthetic observations (±1σ in log space)")
        right.set(xlabel="Wenner electrode spacing a (m)", ylabel="Apparent resistivity (Ω m)", title="Surface response | Wenner-alpha")
        right.legend(loc="best", fontsize=9)
        for ax in (left, right):
            ax.grid(alpha=0.18, which="both")
            ax.spines[["top", "right"]].set_visible(False)
        fig.text(0.06, 0.025, "SYNTHETIC DATA  •  Horizontal, isotropic layers  •  Electrode spacing is not depth", fontsize=9, color="#56616b")
        fig.tight_layout(rect=(0, 0.055, 1, 0.91))
        # Fixed metadata makes the SVG reproducible apart from library rendering changes.
        with plt.rc_context({"svg.hashsalt": "geoscope-v0.1"}):
            fig.savefig(target, metadata={"Date": None}, facecolor="white")
        plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description="GeoScope synthetic Wenner sounding")
    parser.add_argument("--config", type=Path, help="JSON experiment configuration")
    parser.add_argument("--output", type=Path, default=Path("results"))
    args = parser.parse_args()
    config = DEFAULT_CONFIG.copy()
    if args.config:
        supplied = json.loads(args.config.read_text(encoding="utf-8"))
        if not isinstance(supplied, dict):
            parser.error("Configuration must be a JSON object")
        unknown = set(supplied) - set(config)
        if unknown:
            parser.error(f"Unknown configuration keys: {sorted(unknown)}")
        config.update(supplied)
    try:
        clean = run(config, args.output)
    except (ValueError, TypeError) as exc:
        parser.error(str(exc))
    print(f"Simulated {clean.size} Wenner measurements.")
    print(f"Saved figure, CSV, and experiment metadata to {args.output.resolve()}")


if __name__ == "__main__":
    main()
