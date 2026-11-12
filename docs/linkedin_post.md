# LinkedIn draft

One of the most useful results from my latest Python project was a warning:
the uncertainty estimates were not ready to trust.

I have been developing GeoScope, a scientific computing project that builds on
my physics background and undergraduate work in electrical resistivity methods
for groundwater assessment.

Using SimPEG for the forward physics, GeoScope simulates layered-earth
soundings, estimates layer properties, and checks Bayesian sampling with R-hat
and effective sample size. It includes a Streamlit dashboard and reproducible
examples.

In the synthetic comparison, the 200 m survey run passed the chosen sampling
checks. The 30 m run did not. Extending that run from 4,000 to 8,000 retained
draws per chain still did not resolve the convergence failure.

That distinction matters. A table of estimated values can look convincing even
when the sampling behind it is unreliable. I have documented the unsuccessful
case alongside the passed baseline, and updated the dashboard to withhold
unvalidated parameter estimates by default.

This is a synthetic software demonstration, not a field-validated groundwater
tool. The short-survey sampling problem remains open, and I am continuing to
learn how to test both the physics and the inference carefully.

Code, figures, methods, and verified results:
https://github.com/TheABZ/geoscope-resistivity

I would welcome feedback from researchers working in geophysics, scientific
Python, and uncertainty analysis.

#Geophysics #Python #ScientificComputing #InverseProblems #OpenSource
