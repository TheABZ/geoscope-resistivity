# Lesson 2 — Recovering the Earth and measuring uncertainty

## 1. Reverse the problem

In Lesson 1 we chose underground layers and predicted measurements. Now we
pretend the layer values are unknown. The inversion starts with guesses and
changes them to match the observations. The clean synthetic model is not
provided to the optimiser.

Run `python -m geoscope.analysis --sensitivity` with your virtual environment's
Python, as shown in the README. Read `results/analysis/multistart_fits.csv`.
Each row is an independent starting guess and its final fit. Compare the
chi-squared values and layer parameters, not just the success flag.

## 2. A good fit does not prove a unique Earth

The objective is the sum of squared residuals in log space, weighted by the
measurement errors. A weighted RMS near one means residuals have a scale
similar to the assumed errors. A smaller number is not proof of a correct
geological interpretation, and a value below one can occur by chance.

Read `sensitivity.csv`. The short-extent example drives the third-layer
resistivity to the bound at 3000 Ω m, even though the invented truth is 800 Ω m.
Its sounding still fits well. This illustrates insufficient sensitivity to a
deep property for this specific experiment.

## 3. What the posterior means

Bayesian sampling explores models under the stated prior and likelihood.
Here the prior gives equal density per log unit within fixed positive bounds;
it does not assign equal density per Ω m. The likelihood assumes independent
Gaussian errors in log apparent resistivity with known standard deviations.

`posterior_summary.csv` reports medians and the 2.5th/97.5th percentiles for
each parameter. These are marginal credible intervals conditional on these
assumptions. They do not include uncertainty about the chosen layer count or
whether a 1D model describes a real site.

## 4. Check the sampler before trusting its intervals

Four chains should explore similar distributions. Open `posterior_traces.svg`
and inspect the traces. R-hat compares within-chain and between-chain variation;
ESS estimates how much independent information the autocorrelated draws contain.
Every parameter must have R-hat ≤1.01 and bulk/tail ESS ≥400 for the status to
pass. If not, increase warmup/draws and investigate other starting regions,
priors, and model adequacy. Passing diagnostics still cannot prove all modes
were found. These chains start near the best least-squares solution.

## 5. Two different uncertainty bands

The latent-response band shows variation in noise-free predictions across
sampled Earth models. The future-observation band also adds measurement noise.
Neither is a direct underground depth interval. Curve bands are pointwise
95% intervals, not a simultaneous 95% envelope for the whole curve.

## 6. Use the dashboard

Start the dashboard using the README command. Change one parameter at a time,
then compare truth and recovered values. Reduce maximum spacing to 30 m.
Run uncertainty analysis and inspect the trace for rho_3. This 30 m example is
a known difficult case for the current sampler. If diagnostics fail, retain the
outputs as a failed sampling case; do not interpret the parameter intervals or
keep extending runs without a validated sampling strategy. Sampling failure
alone does not prove weak geological identifiability. Changed inputs clear
the previous samples. In v0.2.1, failed runs show diagnostic metrics and withhold
the estimated parameters unless a debugging checkbox is explicitly selected.

Try explaining why deeper-layer uncertainty can increase without a dramatic
change in the fitted sounding curve.
