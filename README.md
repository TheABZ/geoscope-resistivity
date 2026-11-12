# GeoScope — Resistivity

**How much can electrical measurements tell us about the Earth beneath us?**

GeoScope is a Python project exploring layered-earth electrical resistivity,
inverse problems, and eventually uncertainty. It builds on Gabriel Chiazor
Unomah's physics and groundwater geophysics background at the University of Lagos.

**Current milestone: v0.1 — forward modelling.** Generate a three-layer synthetic
Earth, simulate a Wenner-alpha sounding with SimPEG, add reproducible measurement
noise, and export a figure and dataset. Inversion, Bayesian sampling, and the
interactive dashboard are planned; they are not implemented in this version.

![Synthetic layered Earth and Wenner sounding](results/first_sounding.svg)

## Run it on Windows

Install Python 3.11 or 3.12 and Git first. Open **Command Prompt**, then run each
line in order:

```bat
git clone https://github.com/TheABZ/geoscope-resistivity.git
cd geoscope-resistivity
py -3.11 -m venv .venv
.venv\Scripts\python -m pip install -e ".[dev]"
.venv\Scripts\python -m geoscope.demo --config examples/three_layer.json --output results
```

Use `py -3.12` instead if Python 3.12 is installed. Explicit interpreter paths
avoid PowerShell activation-policy issues. In PowerShell, prefix the interpreter
path with `.\` (for example `.\.venv\Scripts\python`).

Expected output:

```text
Simulated 40 Wenner measurements.
Saved figure, CSV, and experiment metadata to ...results
```

Open `results/first_sounding.svg` in a browser. Open the CSV in Excel or Python.
The generated files are:

| File | Contents |
| --- | --- |
| `results/first_sounding.svg` | True layer properties beside the simulated sounding |
| `results/synthetic_sounding.csv` | Spacing, clean response, noisy response, log-space standard deviation |
| `results/experiment.json` | Model, noise assumptions, random seed, software versions |

## Run it on Linux or macOS

```sh
git clone https://github.com/TheABZ/geoscope-resistivity.git
cd geoscope-resistivity
python3 -m venv .venv
.venv/bin/python -m pip install -e '.[dev]'
.venv/bin/python -m geoscope.demo --config examples/three_layer.json --output results
.venv/bin/python -m pytest -q
```

After installation, `geoscope --output results` is also available when the virtual
environment is activated. `python -m geoscope.demo` works without activation when
using the environment's Python as above. Python >=3.10 is declared; 3.11/3.12 are
the recommended starting environments. See `requirements-tested.txt` for the
exact dependency versions used for the committed example.

## The experiment

| Layer | Resistivity (Ω m) | Thickness (m) |
| --- | ---: | ---: |
| 1 | 100 | 5 |
| 2 | 25 | 20 |
| 3 | 800 | Infinite half-space |

The survey uses 40 logarithmically spaced Wenner spacings from 1 to 200 m.
All four electrodes lie on a flat surface in the order A–M–N–B, equally spaced
by `a`. The geometric factor is `K = 2πa`; apparent resistivity is
`rho_a = K * ΔV / I`. SimPEG computes the response of the layered medium.

Noise is Gaussian in log space:
`log(rho_observed / rho_clean) ~ Normal(0, sigma_log²)`, with `sigma_log=0.03`
and seed 42. This is about 3% fractional noise for small errors; the noise-free
prediction is the distribution's median. Figure error bars show a one-standard-
deviation multiplicative noise scale around each observation. **They are not
uncertainty intervals for underground layer properties.**

Edit `examples/three_layer.json` to change resistivities, thicknesses, electrode
spacing, or noise. Keep N resistivities and N−1 thicknesses, all positive. Run
with `--output results/my_experiment` to preserve the original example.

## Scientific scope and checks

- Every included measurement is **synthetic**, not an Eruwa field measurement.
- Layers are horizontal, laterally infinite, isotropic, and homogeneous within
  each layer. The surface is flat and electrodes are ideal point electrodes.
- Electrode spacing controls sensitivity; it is not a measurement of depth.
- Low resistivity alone does not identify water. Clay, salinity, lithology, and
  other factors matter. This tool does not estimate borehole yield or recommend
  drilling locations.
- Tests compare the numerical forward solver against a homogeneous half-space
  and an independent analytical two-layer image series, then check invisible
  interfaces, resistivity scaling, noise reproducibility, invalid inputs, and
  the complete export workflow.

Windows test command: `.venv\Scripts\python -m pytest -q`.

## Development path

- [x] Forward model, synthetic noise, exports, physics benchmarks
- [ ] Bounded inversion of resistivity and thickness in log coordinates
- [ ] Multiple starting models and survey/noise sensitivity experiments
- [ ] Bayesian sampling with explicit priors and convergence diagnostics
- [ ] Interactive dashboard
- [ ] Field-data case study after data and metadata review

Read [Lesson 1](docs/lesson_01.md) for an introduction and a small experiment.
Read [the methods](docs/methods.md) for numerical assumptions and references.

## Attribution and licence

Project author: **Gabriel Chiazor Unomah**. MIT licence; see [LICENSE](LICENSE).
The forward solver is provided by **SimPEG**, with its own licence and citation
requirements. GeoScope is a learning and research software project, not a claim
of a new inversion algorithm or a peer-reviewed scientific finding.
