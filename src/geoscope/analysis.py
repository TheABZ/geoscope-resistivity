"""Complete reproducible inversion, sampling, and sensitivity experiment."""

import argparse
import csv
import json
from importlib.metadata import version
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .demo import DEFAULT_CONFIG, run
from .forward import LayeredEarth, synthetic_observations, simulate_wenner
from .inverse import Sounding, InverseProblem
from .posterior import sample_posterior


def read_sounding(path):
    table = np.genfromtxt(path, names=True, delimiter=",")
    required = {"a_m", "rhoa_observed_ohm_m", "sigma_log"}
    if not required.issubset(table.dtype.names or ()):
        raise ValueError(f"CSV must contain {sorted(required)}")
    return Sounding(*(np.atleast_1d(table[key]) for key in
                     ("a_m", "rhoa_observed_ohm_m", "sigma_log")))


def write_rows(path, rows):
    with Path(path).open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def response_bands(problem, log_samples, seed=44):
    flat = log_samples.reshape(-1, log_samples.shape[-1])
    rng = np.random.default_rng(seed)
    selected = flat[rng.choice(len(flat), min(1000, len(flat)), replace=False)]
    responses = np.array([problem.operator(m) for m in selected])
    # Separate latent-response uncertainty from future noisy measurement uncertainty.
    future = responses * np.exp(rng.normal(size=responses.shape) * problem.data.sigma_log)
    return np.quantile(responses, [0.025, 0.5, 0.975], axis=0), np.quantile(future, [0.025, 0.975], axis=0)


def make_figures(problem, fit, samples, report, output):
    response, future = response_bands(problem, samples)
    a = problem.data.spacings
    order = np.argsort(a)
    with plt.rc_context({"font.size": 10, "svg.fonttype": "none", "svg.hashsalt": "geoscope-v02"}):
        fig, axes = plt.subplots(1, 2, figsize=(12, 5), layout="constrained")
        left, right = axes
        left.fill_between(a[order], future[0, order], future[1, order], color="#f4bc69", alpha=.35, label="95% future-observation interval")
        left.fill_between(a[order], response[0, order], response[2, order], color="#0b7285", alpha=.3, label="95% latent-response interval")
        left.loglog(a, problem.data.observed, ".", color="#b66a13", label="Observations")
        left.loglog(a[order], problem.operator(fit.x)[order], color="#0b7285", label="Best least-squares fit")
        left.set(xlabel="Wenner spacing a (m)", ylabel="Apparent resistivity (Ω m)", title="Data fit and posterior response")
        left.legend(fontsize=8)
        # Compare parameters on a dimensionless log ratio axis to avoid mixing units.
        rows = report["parameters"]
        center = np.array([r["median"] for r in rows])
        low = np.array([r["q025"] for r in rows]) / center
        high = np.array([r["q975"] for r in rows]) / center
        right.errorbar(np.ones(len(rows)), np.arange(len(rows)), xerr=[1-low, high-1], fmt="o", capsize=4, color="#0b7285")
        right.set_yticks(np.arange(len(rows)), [r["parameter"] for r in rows])
        right.set(xscale="log", xlabel="Parameter / posterior median", title="95% marginal credible intervals")
        for axis in axes:
            axis.grid(alpha=.18)
        fig.suptitle(f"GeoScope | {report['status']} | Conditional on chosen layers and priors")
        fig.savefig(output / "inversion_uncertainty.svg", metadata={"Date": None})
        plt.close(fig)
        fig, axes = plt.subplots(len(problem.names), 1, figsize=(10, 1.6 * len(problem.names)), sharex=True, layout="constrained", squeeze=False)
        stride = max(1, samples.shape[1] // 350)
        for i, name in enumerate(problem.names):
            axis = axes[i, 0]
            for chain in range(samples.shape[0]):
                axis.plot(np.arange(0, samples.shape[1], stride), np.exp(samples[chain, ::stride, i]), lw=.6, alpha=.7, label=f"Chain {chain+1}")
            axis.set_ylabel(name)
        axes[0, 0].legend(ncol=4, fontsize=8)
        axes[-1, 0].set_xlabel("Retained iteration (warmup removed; display thinned only)")
        fig.suptitle("Independent posterior chains | rho in Ω m; h in m")
        fig.savefig(output / "posterior_traces.svg", metadata={"Date": None})
        plt.close(fig)
    return [{"a_m": float(a[i]), "latent_q025": float(response[0, i]),
             "latent_median": float(response[1, i]), "latent_q975": float(response[2, i]),
             "future_q025": float(future[0, i]), "future_q975": float(future[1, i])} for i in range(len(a))]


def sensitivity_experiment(seed=42):
    """Paired-noise synthetic benchmark; no general reliability claim."""
    earth = LayeredEarth((100, 25, 800), (5, 20))
    rows = []
    for max_spacing in (30., 200.):
        a = np.geomspace(1, max_spacing, 40)
        clean = simulate_wenner(earth, a)
        for noise in (.01, .03, .1):
            observed, sigma = synthetic_observations(clean, noise, seed)
            problem = InverseProblem(Sounding(a, observed, sigma))
            best, _ = problem.fit(starts=5, seed=seed)
            rows.append({"max_spacing_m": max_spacing, "sigma_log": noise,
                         "weighted_rms": float(np.sqrt(np.mean(best.fun**2))),
                         "jacobian_condition_number": float(np.linalg.cond(best.jac)),
                         **{name: float(value) for name, value in zip(problem.names, np.exp(best.x))}})
    return rows


def analyze(data, output, layers=3, starts=8, draws=4000, warmup=2000, seed=42, origin="user_supplied"):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    problem = InverseProblem(data, n_layers=layers)
    fit, fits = problem.fit(starts=starts, seed=seed)
    samples, report = sample_posterior(problem, fit, draws=draws, warmup=warmup, seed=seed+1)
    report.update({"data_origin": origin, "model_layers": layers,
                   "prior": "Independent uniform log(rho) in log([10,3000] Ω m); log(h) in log([0.5,100] m)",
                   "likelihood": "Independent Gaussian log-residuals with known sigma_log",
                   "best_fit": {name: float(v) for name, v in zip(problem.names, np.exp(fit.x))},
                   "weighted_rms": float(np.sqrt(np.mean(fit.fun**2))),
                   "versions": {name: version(name) for name in ("geoscope-resistivity", "simpeg", "numpy", "scipy", "arviz")}})
    # Portable numerical outputs; all retained samples saved without thinning.
    np.savez_compressed(output / "posterior_samples.npz", log_parameters=samples,
                        parameter_names=np.array(problem.names))
    write_rows(output / "posterior_summary.csv", report["parameters"])
    write_rows(output / "multistart_fits.csv", [{"start": i+1, "success": bool(r.success),
               "chi_squared": float(2*r.cost), "nfev": int(r.nfev),
               **{name: float(v) for name, v in zip(problem.names, np.exp(r.x))}} for i, r in enumerate(fits)])
    write_rows(output / "response_intervals.csv", make_figures(problem, fit, samples, report, output))
    # Strict JSON: nonfinite diagnostic values are null rather than invalid NaN tokens.
    def sanitize(value):
        if isinstance(value, float) and not np.isfinite(value):
            return None
        if isinstance(value, dict):
            return {k: sanitize(v) for k, v in value.items()}
        if isinstance(value, list):
            return [sanitize(v) for v in value]
        return value
    (output / "analysis.json").write_text(json.dumps(sanitize(report), indent=2, allow_nan=False) + "\n", encoding="utf-8")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, help="Wenner CSV; otherwise use the synthetic example")
    parser.add_argument("--output", type=Path, default=Path("results/analysis"))
    parser.add_argument("--layers", type=int, default=3)
    parser.add_argument("--starts", type=int, default=8)
    parser.add_argument("--draws", type=int, default=4000)
    parser.add_argument("--warmup", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--sensitivity", action="store_true", help="Also run six synthetic noise/aperture experiments")
    args = parser.parse_args()
    try:
        if args.data:
            data = read_sounding(args.data)
            origin = "user_supplied; not independently verified as field data"
        else:
            run(DEFAULT_CONFIG, args.output / "input")
            data = read_sounding(args.output / "input/synthetic_sounding.csv")
            origin = "synthetic"
        report = analyze(data, args.output, args.layers, args.starts, args.draws, args.warmup, args.seed, origin)
        if args.sensitivity:
            write_rows(args.output / "sensitivity.csv", sensitivity_experiment(args.seed))
    except (ValueError, RuntimeError) as exc:
        parser.error(str(exc))
    print(f"Status: {report['status']}; weighted RMS: {report['weighted_rms']:.3f}")
    print(f"Saved analysis to {args.output.resolve()}")


if __name__ == "__main__":
    main()
