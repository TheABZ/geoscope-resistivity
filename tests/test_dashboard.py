import pytest
from pathlib import Path

pytest.importorskip("streamlit")
from streamlit.testing.v1 import AppTest
APP = Path(__file__).resolve().parents[1] / "streamlit_app.py"


def test_dashboard_simulates_and_refits():
    app = AppTest.from_file(APP, default_timeout=60).run()
    assert not app.exception
    assert app.metric[0].label == "Weighted residual RMS"
    old = app.dataframe[0].value["True"].iloc[0]
    app.sidebar.number_input[1].set_value(150.).run()
    assert not app.exception
    assert app.dataframe[0].value["True"].iloc[0] != old


def test_dashboard_sampling_and_stale_state():
    app = AppTest.from_file(APP, default_timeout=90).run()
    app.selectbox[0].select(2000)
    app.button[0].click().run(timeout=90)
    assert not app.exception
    assert "posterior_result" in app.session_state
    assert len(app.dataframe) >= 3
    app.sidebar.slider[0].set_value(100).run()
    assert not app.exception
    assert "posterior_result" not in app.session_state


def test_upload_mode_waits_for_data():
    app = AppTest.from_file(APP, default_timeout=60).run()
    app.sidebar.radio[0].set_value("Upload Wenner CSV").run()
    assert not app.exception
    assert "Upload a Wenner" in app.info[0].value
