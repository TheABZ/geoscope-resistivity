# Suggested LinkedIn post

What can electrical measurements tell us about the ground beneath our feet—and
what can they leave uncertain?

I have been developing GeoScope, a Python project that connects my physics and
groundwater geophysics background with computational modelling and uncertainty
analysis.

Using SimPEG for the forward physics, the project simulates layered-earth
resistivity soundings, recovers model parameters from several starting guesses,
and explores Bayesian uncertainty with convergence diagnostics. It also includes
an interactive Streamlit dashboard and reproducible examples.

One result from the synthetic demonstration stood out: limiting electrode
spacing to 30 m produced a good fit while pushing the estimated deep-layer
resistivity to its upper bound. The lesson is useful: matching observations does
not necessarily mean we have constrained the underlying Earth well.

This is a developing research and learning project using synthetic data. I am
continuing to build my understanding of inverse problems, validation, and the
assumptions behind scientific software.

Code, examples, and methods:
https://github.com/TheABZ/geoscope-resistivity

I would welcome feedback from people working in geophysics, scientific Python,
and uncertainty quantification.

#Geophysics #Python #ScientificComputing #InverseProblems #OpenSource
