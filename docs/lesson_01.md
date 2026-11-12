# Lesson 1 — From underground properties to surface measurements

The aim today is to understand **forward modelling**: start with an underground
model, then calculate what the measuring instrument should observe.

## 1. Describe the Earth

The example has a 5 m thick upper layer, a 20 m thick middle layer, and a lower
half-space. The boundaries are therefore at 5 m and 25 m. Their resistivities
are 100, 25, and 800 Ω m. Resistivity describes how strongly a material opposes
electrical current flow. These numbers are an invented teaching example.

## 2. Describe the measurement

A and B inject and withdraw electrical current; M and N measure potential
difference. In the Wenner-alpha array all adjacent electrodes are separated by
the same distance a. For each a, SimPEG predicts the response of all the layers.

The instrument response is converted to **apparent resistivity**. This is the
resistivity a uniform Earth would need to produce the same response for that
electrode arrangement. It is not a direct reading of one underground layer.

## 3. Read the figure

The left panel shows the model we supplied, including exact boundary depths.
The right panel shows the surface sounding. Its horizontal axis is electrode
spacing, not depth. Its vertical axis is apparent resistivity, not a layer's
true resistivity. Larger spacings generally change sensitivity toward deeper
structures, but there is no one-to-one spacing-to-depth conversion here.

The smooth line is the numerical response without noise. The points are
synthetic noisy observations, generated with a fixed seed so you can reproduce
the exercise. Their bars describe the prescribed measurement-noise scale.

## 4. Try one change

1. Run the default example using the README instructions.
2. Copy `examples/three_layer.json` to `examples/uniform.json`.
3. In the copy, set `resistivities_ohm_m` to `[100, 100, 100]`.
4. Run the same command with `--config examples/uniform.json --output results/uniform`.
5. Compare the two figures.

Expected result: the noise-free uniform-Earth curve is horizontal at 100 Ω m.
Layer boundaries have no electrical contrast, so the sounding cannot see them.
Noisy points still scatter around that line.

## 5. The question that leads to inversion

In the field, we know electrode positions and measurements but not the true
layer properties. Inversion tries to work backwards. Our next milestone will
estimate model parameters and examine whether different models can explain the
same observations. A good fit alone will not prove the model is unique.

Before moving on, explain in your own words: why can a sounding curve differ
from the resistivity profile, and why do electrically identical layers produce
no boundary signature?
