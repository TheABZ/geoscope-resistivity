# Methods and numerical validation

## Forward problem

The DC conduction equation away from point current sources is
`div(sigma * grad(phi)) = 0`, where `sigma = 1/rho` and `phi` is electric
potential. Point sources provide the current boundary/source terms. The surface
is insulating away from electrodes, and the medium extends to infinite depth.

GeoScope constructs a SimPEG `Survey` of dipole sources and dipole receivers,
requests apparent-resistivity data, and uses `Simulation1DLayers` with explicit
resistivity and thickness arrays. SimPEG 0.25.0 supplies its default layered-earth
digital-filter Hankel-transform calculation. No mesh, inversion regularisation,
or parameter estimation is performed in the forward-only v0.1 workflow. The
v0.2 inverse workflow below uses the same solver with exponential parameter
maps, reusing only survey geometry between predictions.

For A, M, N, B at `(-1.5a, -0.5a, 0.5a, 1.5a)`, respectively, the surface
half-space geometric factor is `2πa`. Each source uses 1 A. A uniform half-space
must therefore yield `rho_a = rho` independently of spacing.

## Independent two-layer benchmark

For upper-layer resistivity rho1, thickness h, lower half-space resistivity rho2,
define `k = (rho2-rho1)/(rho2+rho1)`. The surface point-source image expansion
gives the Wenner apparent resistivity:

```math
\rho_a = \rho_1\left[1+4\sum_{n=1}^{\infty} k^n\left(
\frac{a}{\sqrt{a^2+(2nh)^2}}-
\frac{a}{\sqrt{4a^2+(2nh)^2}}\right)\right].
```

The test evaluates 1,000 image terms for rho1=100 Ω m, h=5 m and rho2=25 or
800 Ω m. These choices have |k|<0.78, making series truncation negligible
relative to the numerical acceptance tolerance. Spacings range from 0.2 to
500 m. The relative tolerance is 0.02%. This benchmark checks array geometry,
units, contrast behaviour, and the numerical solver independently of its code.
It does not validate all possible extreme models or establish real-world accuracy.

Other tests check a uniform half-space, identical adjacent layers, scaling all
resistivities by the same factor, noise generation, input validation, and exports.

## Inversion and Bayesian inference (v0.2)

The parameter vector is `m = [log(rho_1), ..., log(rho_N), log(h_1), ...,
log(h_(N-1))]`. Residuals are `(log(prediction)-log(observation))/sigma_log`.
SciPy's bounded `least_squares` uses eight starts by default, up to 500 function
evaluations per start, and ftol/xtol/gtol of 1e-9. Only converged finite fits
are eligible as the best fit. There is no smoothing penalty. Bounds are
rho in [10,3000] Ω m and h in [0.5,100] m.

The posterior density in m is proportional to exp(-0.5*sum(residual²)) inside
these bounds and zero outside, corresponding to independent uniform priors in
log coordinates. Measurement uncertainties are fixed, so likelihood constants
do not affect sampling. No Jacobian term is needed for a prior uniform in m;
the corresponding physical-parameter prior is proportional to 1/(rho*h).

Four independent Metropolis chains use SeedSequence children. Initialization is
dispersed near the best fit using four times the covariance from a regularized
inverse J.T@J; eigenvalues are floored at 1e-3. Proposals scale covariance by
2.38²/p. After 300 warmup steps, covariance is refreshed every 100 warmup steps
using 95% empirical covariance, 5% initial covariance, and small diagonal jitter.
Adaptation stops at the end of warmup; no adaptation occurs in retained draws.
Default warmup is 2000 and retained length 4000 per chain. Saved NPZ arrays keep
every retained sample in log coordinates; trace figures thin only the display.

ArviZ provides rank-normalized split R-hat (with folding), and bulk/tail ESS.
All R-hat values must be finite and ≤1.01, and both ESS estimates ≥400, to label
sampling diagnostics passed. Nonconvergence is retained as an explicit result.
Other posterior modes, incorrect layer count, uncertain noise, and violations
of 1D geometry can invalidate an interpretation even when diagnostics pass.

Response intervals use up to 1000 randomly selected draws. Noise-free curves
give pointwise latent-response intervals; an additional draw from the assumed
noise model gives pointwise future-observation intervals. The latter include
both parameter and assumed measurement-noise variation. Sampling validation
compares a uniform-Earth log-resistivity posterior with its analytic Normal
mean and variance under a broad bounded log-uniform prior.

## Reproducibility

The committed example includes configuration, seed, package versions, tabulated
clean/noisy responses, and SVG. Floating-point responses can vary slightly by
platform. The random-number generator is NumPy's default generator with seed 42.
The noise model is specified in the README and metadata. Regeneration overwrites
the three named outputs within the selected output directory.

## Sources

- [SimPEG source and citation guidance](https://github.com/simpeg/simpeg).
- [SimPEG user tutorials](https://simpeg.xyz/user-tutorials/).
- [SciPy least_squares](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.least_squares.html).
- [ArviZ rank R-hat](https://python.arviz.org/en/v0.23.4/api/generated/arviz.rhat.html).
- [ArviZ effective sample size](https://python.arviz.org/en/v0.23.4/api/generated/arviz.ess.html).
- [Streamlit app testing](https://docs.streamlit.io/develop/api-reference/app-testing).
- [Versioned layered-earth solver source](https://github.com/simpeg/simpeg/blob/v0.25.0/simpeg/electromagnetics/static/resistivity/simulation_1d.py).
- Cockett et al. (2015), *SimPEG: An open source framework for simulation and
  gradient based parameter estimation in geophysical applications*, Computers
  & Geosciences 85, 142–154. [DOI](https://doi.org/10.1016/j.cageo.2015.09.015).

GeoScope combines experiment setup, multistart optimisation, a Metropolis
sampling workflow, diagnostic handling, exports, dashboard, and learning
materials. The forward numerical method is SimPEG's; optimisation is SciPy's
and convergence diagnostics are ArviZ's. No new geophysical algorithm is claimed.
