# Reproducible synthetic demonstration

Command: `python -m geoscope.analysis --sensitivity` with the project's installed
virtual environment. Numerical versions are recorded in `analysis.json`.

The reference experiment has true parameters [100, 25, 800] Ω m and [5, 20] m,
40 Wenner spacings from 1–200 m, log-noise standard deviation 0.03, and seed 42.
All measurements in this demonstration are simulated.

## Parameter recovery and conditional uncertainty

| Parameter | Truth | Posterior median | 95% marginal interval | Unit |
| --- | ---: | ---: | --- | --- |
| rho_1 | 100 | 99.20 | 97.31–101.19 | Ω m |
| rho_2 | 25 | 23.00 | 18.31–26.43 | Ω m |
| rho_3 | 800 | 725.40 | 607.58–915.56 | Ω m |
| h_1 | 5 | 5.23 | 4.84–5.67 | m |
| h_2 | 20 | 17.86 | 13.76–21.19 | m |

Weighted residual RMS: approximately 0.790. All five reference parameters fall
within their marginal intervals in this single run; that is not a coverage
study. Across the retained samples, maximum rank R-hat is about 1.006, minimum
bulk ESS about 624, and minimum tail ESS about 417. These satisfy the declared
thresholds. Chains start near the best fit, so other modes may remain unexplored.

## Why survey extent matters in this example

The sensitivity experiment compares maximum spacings of 30 and 200 m at three
noise levels: 1%, 3%, and 10% in log space. The same random seed pairs the noise
patterns across cases; these six experiments are not independent statistical
replicates. All cases use the same truth and 40 measurements.

For all three 30 m cases, recovered rho_3 reaches the 3000 Ω m upper bound.
Their weighted RMS values are still about 0.8. The condition number of the
weighted residual Jacobian is above 4300, versus approximately 24–46 for the
200 m cases. This illustrates poor conditioning for this particular short
survey. It does not provide a general minimum spacing rule or an uncertainty
interval by itself. Inspect the exported values before generalizing.

## Validation limits

The forward solver passes independent uniform-half-space and analytical
two-layer tests. The inversion recovers an exact synthetic three-layer model
without noise. The sampler is checked against an analytically known one-layer
posterior and deliberately separated chains. Dashboard tests exercise
simulation, changed inputs, sampling, and clearing stale results. Field-data
validation, model comparison, and repeated-simulation coverage studies remain
future research work.
