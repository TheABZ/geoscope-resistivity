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
or parameter estimation is performed in v0.1.

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

## Reproducibility

The committed example includes configuration, seed, package versions, tabulated
clean/noisy responses, and SVG. Floating-point responses can vary slightly by
platform. The random-number generator is NumPy's default generator with seed 42.
The noise model is specified in the README and metadata. Regeneration overwrites
the three named outputs within the selected output directory.

## Sources

- [SimPEG source and citation guidance](https://github.com/simpeg/simpeg).
- [SimPEG user tutorials](https://simpeg.xyz/user-tutorials/).
- [Versioned layered-earth solver source](https://github.com/simpeg/simpeg/blob/v0.25.0/simpeg/electromagnetics/static/resistivity/simulation_1d.py).
- Cockett et al. (2015), *SimPEG: An open source framework for simulation and
  gradient based parameter estimation in geophysical applications*, Computers
  & Geosciences 85, 142–154. [DOI](https://doi.org/10.1016/j.cageo.2015.09.015).

GeoScope's contribution at this stage is its experiment setup, exports,
independent checks, and learning materials; the forward numerical method is
SimPEG's.
