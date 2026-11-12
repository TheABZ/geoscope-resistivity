"""Local dashboard: python -m streamlit run streamlit_app.py."""

import io
import json

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import streamlit as st

from geoscope.forward import LayeredEarth, simulate_wenner, synthetic_observations
from geoscope.inverse import Sounding, InverseProblem
from geoscope.posterior import sample_posterior

st.set_page_config(page_title="GeoScope | Resistivity", page_icon="🌍", layout="wide")
st.title("GeoScope")
st.subheader("Explore what electrical soundings reveal beneath the surface")
st.caption("Synthetic experiments and Wenner-alpha CSV analysis • Gabriel Chiazor Unomah")

with st.sidebar:
    st.header("Experiment")
    mode = st.radio("Data source", ["Synthetic three-layer Earth", "Upload Wenner CSV"])
    seed = st.number_input("Random seed", min_value=0, max_value=1000000, value=42, step=1)
    settings = []
    if mode.startswith("Synthetic"):
        for i, default in enumerate((100., 25., 800.)):
            settings.append(st.number_input(f"Layer {i+1} resistivity (Ω m)", min_value=10., max_value=3000., value=default))
        for i, default in enumerate((5., 20.)):
            settings.append(st.number_input(f"Layer {i+1} thickness (m)", min_value=.5, max_value=100., value=default))
        maximum = st.slider("Maximum electrode spacing (m)", 20, 500, 200)
        noise = st.slider("Log-noise standard deviation (%)", 1, 15, 3) / 100
        uploaded = None
        layers = 3
    else:
        uploaded = st.file_uploader("Wenner CSV", type="csv")
        layers = st.selectbox("Assumed layer count", [1, 2, 3], index=2)
        st.caption("Columns: a_m, rhoa_observed_ohm_m, sigma_log. All values must be positive. sigma_log is a log-space standard deviation.")


@st.cache_data(show_spinner=False)
def fit_data(a, observed, sigma, layers, seed):
    problem = InverseProblem(Sounding(a, observed, sigma), layers)
    best, fits = problem.fit(starts=6, seed=seed)
    return best, fits, problem.operator(best.x)


if mode.startswith("Synthetic"):
    a = np.geomspace(1, maximum, 40)
    clean = simulate_wenner(LayeredEarth(settings[:3], settings[3:]), a)
    observed, sigma = synthetic_observations(clean, noise, seed)
    st.info("These are simulated measurements. The supplied layer properties are known only because we created this experiment.")
else:
    if uploaded is None:
        st.info("Upload a Wenner-alpha sounding CSV to begin. Other array geometries need a different forward model.")
        st.stop()
    try:
        table = pd.read_csv(uploaded)
        a, observed, sigma = (table[c].to_numpy(dtype=float) for c in ("a_m", "rhoa_observed_ohm_m", "sigma_log"))
        Sounding(a, observed, sigma)
    except (KeyError, ValueError, TypeError) as exc:
        st.error(f"Invalid sounding CSV: {exc}")
        st.stop()
    clean = None
    st.info("Uploaded measurements have not been independently verified. Confirm Wenner geometry, units, and error estimates before interpreting results.")

try:
    with st.spinner("Fitting layer properties from six starting models…"):
        fit, fits, prediction = fit_data(a, observed, sigma, layers, seed)
except (ValueError, RuntimeError) as exc:
    st.error(str(exc))
    st.stop()
problem = InverseProblem(Sounding(a, observed, sigma), layers)
left, right = st.columns([1.6, 1])
with left:
    fig, axis = plt.subplots(figsize=(8, 4.5))
    order = np.argsort(a)
    axis.loglog(a, observed, "o", ms=4, color="#b66a13", label="Observations")
    if clean is not None:
        axis.loglog(a, clean, "--", color="#999999", label="True noise-free response")
    axis.loglog(a[order], prediction[order], color="#0b7285", label="Recovered response")
    axis.set(xlabel="Wenner spacing a (m)", ylabel="Apparent resistivity (Ω m)")
    axis.legend(); axis.grid(alpha=.2)
    st.pyplot(fig)
    plt.close(fig)
with right:
    st.metric("Weighted residual RMS", f"{np.sqrt(np.mean(fit.fun**2)):.3f}")
    rows = [{"Parameter": name, "Recovered": v, "Unit": "Ω m" if name.startswith("rho") else "m"}
            for name, v in zip(problem.names, np.exp(fit.x))]
    if clean is not None:
        for row, true in zip(rows, settings):
            row["True"] = true
    st.dataframe(pd.DataFrame(rows), hide_index=True)
st.caption("A small residual does not prove that the model is unique. Electrode spacing is not depth. Fixed layer count and parameter bounds affect the result.")
csv_data = pd.DataFrame({"a_m": a, "rhoa_observed_ohm_m": observed, "sigma_log": sigma}).to_csv(index=False)
st.download_button("Download sounding CSV", csv_data, "sounding.csv", "text/csv")
with st.expander("Compare starting models"):
    st.dataframe(pd.DataFrame([{"Start": i+1, "Converged": r.success, "Chi squared": 2*r.cost,
                               **dict(zip(problem.names, np.exp(r.x)))} for i, r in enumerate(fits)]), hide_index=True)

st.header("How uncertain are the recovered properties?")
st.write("Sample four independent chains under explicit log-uniform priors: resistivity 10–3000 Ω m and thickness 0.5–100 m. Noise estimates and layer count are held fixed.")
draws = st.selectbox("Retained draws per chain", [2000, 4000, 8000], index=1)
# Results are tied to exact input bytes, layer count, seed, and sampling settings.
fingerprint = (csv_data, layers, seed, draws)
if st.session_state.get("posterior_input") != fingerprint:
    st.session_state.pop("posterior_result", None)
if st.button("Run uncertainty analysis", type="primary"):
    with st.spinner("Sampling four chains and checking convergence…"):
        samples, report = sample_posterior(problem, fit, draws=draws, warmup=2000, seed=seed+1)
    st.session_state.posterior_result = (samples, report)
    st.session_state.posterior_input = fingerprint
if "posterior_result" in st.session_state:
    samples, report = st.session_state.posterior_result
    if report["status"] == "diagnostics_passed":
        st.success("Sampling diagnostics passed. Interpretation remains conditional on the chosen model and priors.")
    else:
        st.warning("Sampling is inconclusive. Treat intervals as provisional; increase draws and inspect chain mixing.")
    st.dataframe(pd.DataFrame(report["parameters"]), hide_index=True)
    st.caption("q025–q975 are 95% marginal credible intervals. R-hat ≤1.01 and bulk/tail ESS ≥400 are required for every parameter.")
    param = st.selectbox("Inspect parameter traces", problem.names)
    i = problem.names.index(param)
    traces = pd.DataFrame(np.exp(samples[:, :, i]).T, columns=[f"Chain {j+1}" for j in range(4)])
    st.line_chart(traces.iloc[::max(1, draws//500)])
    st.download_button("Download diagnostics JSON", json.dumps(report, indent=2), "diagnostics.json", "application/json")
    stream = io.BytesIO()
    np.savez_compressed(stream, log_parameters=samples, parameter_names=np.array(problem.names))
    st.download_button("Download posterior samples", stream.getvalue(), "posterior_samples.npz", "application/octet-stream")
st.caption("Resistivity alone cannot establish groundwater presence or borehole yield. This is research and learning software.")
