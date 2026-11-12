"""Independent adaptive-warmup random-walk Metropolis chains.

All covariance adaptation stops before saved posterior draws. ArviZ provides
rank-normalized split R-hat plus bulk and tail effective sample sizes.
"""

import arviz as az
import numpy as np


def diagnostics(log_samples, names):
    physical = np.exp(log_samples)
    data = az.from_dict(posterior={name: physical[:, :, i] for i, name in enumerate(names)})
    rhat = az.rhat(data, method="rank")
    bulk = az.ess(data, method="bulk")
    tail = az.ess(data, method="tail")
    rows = []
    for i, name in enumerate(names):
        low, median, high = np.quantile(physical[:, :, i], [0.025, 0.5, 0.975])
        rows.append({"parameter": name, "q025": float(low), "median": float(median),
                     "q975": float(high), "rhat": float(rhat[name]),
                     "ess_bulk": float(bulk[name]), "ess_tail": float(tail[name])})
    passed = all(np.isfinite(r["rhat"]) and r["rhat"] <= 1.01
                 and r["ess_bulk"] >= 400 and r["ess_tail"] >= 400 for r in rows)
    return {"status": "diagnostics_passed" if passed else "inconclusive_sampling",
            "criteria": "All parameters: rank R-hat <=1.01; bulk and tail ESS >=400",
            "parameters": rows,
            "limitation": "Passing diagnostics does not prove all posterior modes were explored."}


def sample_posterior(problem, fit, draws=4000, warmup=2000, chains=4, seed=43, progress=None):
    for value, minimum, label in ((draws, 20, "draws"), (warmup, 50, "warmup"), (chains, 2, "chains")):
        if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
            raise ValueError(f"{label} must be an integer >= {minimum}")
    p = len(problem.names)
    # Local information initializes scale only; the posterior target is always
    # evaluated with the full nonlinear forward model and stated bounded prior.
    eigenvalues, eigenvectors = np.linalg.eigh(fit.jac.T @ fit.jac)
    covariance = (eigenvectors * (1 / np.maximum(eigenvalues, 1e-3))) @ eigenvectors.T
    covariance += np.eye(p) * 1e-9
    samples = np.empty((chains, draws, p))
    acceptance = []
    for chain, sequence in enumerate(np.random.SeedSequence(seed).spawn(chains)):
        rng = np.random.default_rng(sequence)
        current = np.clip(rng.multivariate_normal(fit.x, 4 * covariance),
                          problem.lower + 1e-8, problem.upper - 1e-8)
        logp = problem.log_probability(current)
        if not np.isfinite(logp):
            raise RuntimeError("Invalid initial chain state")
        cov = covariance.copy()
        history = []
        accepted = 0
        for step in range(warmup + draws):
            proposal = current + rng.multivariate_normal(np.zeros(p), cov * (2.38**2 / p))
            proposed_logp = problem.log_probability(proposal)
            if np.log(rng.uniform()) < proposed_logp - logp:
                current, logp = proposal, proposed_logp
                if step >= warmup:
                    accepted += 1
            if step < warmup:
                history.append(current.copy())
                if step >= 299 and (step + 1) % 100 == 0:
                    # Diminishing influence of initialization; use all warmup states.
                    empirical = np.cov(np.asarray(history).T)
                    empirical = np.atleast_2d(empirical)
                    cov = 0.95 * empirical + 0.05 * covariance + np.eye(p) * 1e-8
            else:
                samples[chain, step - warmup] = current
        acceptance.append(accepted / draws)
        if progress:
            progress(chain+1, chains)
    report = diagnostics(samples, problem.names)
    report.update({"sampler": "Independent adaptive-warmup random-walk Metropolis; adaptation stops before retained draws",
                   "acceptance_rates": acceptance, "seed": seed, "warmup": warmup,
                   "draws_per_chain": draws, "chains": chains,
                   "initialization": "Dispersed near best multistart fit using local covariance"})
    return samples, report
