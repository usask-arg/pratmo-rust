"""
Integration tests for the pratmo Python wrapper.

Both CTM and DIURN modes work with ``PratmoModel.with_defaults()`` using the
embedded science data (fort01.x, fort02.x, and the spectral/atmosphere files).
"""

import numpy as np
import pytest

from pratmo import (
    BoxSnapshot,
    CtmBoxSpec,
    CtmConfig,
    DiurnBoxSpec,
    DiurnConfig,
    JValues,
    LongLivedMixingRatios,
    PratmoModel,
)


@pytest.fixture(scope="module")
def model():
    return PratmoModel.with_defaults()


@pytest.fixture(scope="module")
def ctm_out_1box(model):
    """Single-box CTM run used by multiple tests."""
    cfg = CtmConfig(
        latitude_deg=60.0,
        julian_day=75,
        integration_days=10,
        boxes=[CtmBoxSpec(altitude_level=20)],
    )
    return model.run_ctm(cfg)


@pytest.fixture(scope="module")
def ctm_out_4box(model):
    """Four-box CTM run for array-shape tests."""
    boxes = [CtmBoxSpec(altitude_level=lvl) for lvl in [10, 15, 20, 25]]
    cfg = CtmConfig(
        latitude_deg=45.0,
        julian_day=180,
        integration_days=5,
        boxes=boxes,
    )
    return model.run_ctm(cfg)


# ── CTM basic ──────────────────────────────────────────────────────────────────

def test_ctm_basic(ctm_out_1box):
    out = ctm_out_1box
    assert len(out.boxes) == 1
    snap = out.boxes[0]
    assert snap.altitude_km > 0
    assert snap.pressure_mb > 0
    assert snap.temperature_k > 100
    assert snap.air_density_cm3 > 0


def test_ctm_box_index(ctm_out_1box):
    assert ctm_out_1box.boxes[0].box_index == 0


def test_ctm_multi_box(ctm_out_4box):
    out = ctm_out_4box
    assert len(out.boxes) == 4
    # Altitude should increase with box index (lower level → lower altitude)
    alts = [b.altitude_km for b in out.boxes]
    assert all(a > 0 for a in alts)


# ── Species profile (1-D numpy arrays from CTM) ────────────────────────────────

def test_species_profile_shape(ctm_out_4box):
    profile = ctm_out_4box.species_profile("o3")
    assert isinstance(profile, np.ndarray)
    assert profile.ndim == 1
    assert len(profile) == 4
    assert np.all(profile >= 0)


def test_species_profile_oh(ctm_out_4box):
    oh = ctm_out_4box.species_profile("oh")
    assert isinstance(oh, np.ndarray)
    assert len(oh) == 4


def test_species_profile_all_names(ctm_out_4box):
    """Verify every valid implicit species name can be queried."""
    valid = [
        "no", "no2", "no3", "n2o5", "hno3", "h", "oh", "ho2", "h2o2",
        "o", "o3", "bro", "br", "hbr", "hno2", "hcl", "cl", "cl2", "clo",
        "clono2", "hno4", "hocl", "brono2", "hobr", "h2co", "ch3o2",
        "ch3o2h", "oclo", "cl2o2", "brcl",
    ]
    assert len(valid) == 30
    for name in valid:
        arr = ctm_out_4box.species_profile(name)
        assert isinstance(arr, np.ndarray), f"species '{name}' did not return ndarray"


def test_species_profile_invalid(ctm_out_1box):
    with pytest.raises(ValueError, match="Unknown"):
        ctm_out_1box.species_profile("not_a_species")


def test_jvalue_profile_shape(ctm_out_4box):
    jp = ctm_out_4box.jvalue_profile("no2")
    assert isinstance(jp, np.ndarray)
    assert jp.ndim == 1
    assert len(jp) == 4
    assert np.all(jp >= 0)


def test_jvalue_profile_all_names(ctm_out_4box):
    """Verify every valid J-value name can be queried."""
    valid = [
        "no", "o2", "o3", "o3_o1d", "h2co_a", "h2co_b", "h2o2", "rooh",
        "no2", "no3_x", "no3_l", "n2o5", "hno2", "hno3", "hno4", "clono2",
        "cl2", "hocl", "oclo", "cl2o2", "clo", "bro", "brono2", "hobr",
        "n2o", "cfc11", "cfc12", "cfc113", "cfc114", "cfc115", "ccl4",
        "ch3cl", "ch3ccl3", "ch3br", "h1211", "h1301", "h2402", "hcfc22",
        "hcfc123", "hcfc141b", "chbr3", "ch3i", "cf3i", "ocs",
    ]
    assert len(valid) == 44
    for name in valid:
        arr = ctm_out_4box.jvalue_profile(name)
        assert isinstance(arr, np.ndarray), f"J-value '{name}' did not return ndarray"


def test_jvalue_profile_invalid(ctm_out_1box):
    with pytest.raises(ValueError, match="Unknown"):
        ctm_out_1box.jvalue_profile("not_a_jvalue")


# ── BoxSnapshot fields ─────────────────────────────────────────────────────────

def test_box_snapshot_implicit(ctm_out_1box):
    imp = ctm_out_1box.boxes[0].implicit
    assert imp.o3 >= 0
    assert isinstance(imp.oh, float)


def test_box_snapshot_long_lived(ctm_out_1box):
    ll = ctm_out_1box.boxes[0].long_lived
    assert ll.o3 >= 0
    assert isinstance(ll.ch4, float)


def test_box_snapshot_jvalues(ctm_out_1box):
    jv = ctm_out_1box.boxes[0].jvalues
    assert jv.no2 >= 0
    assert isinstance(jv.o3_o1d, float)


# ── to_dict methods ────────────────────────────────────────────────────────────

def test_implicit_species_to_dict(ctm_out_1box):
    d = ctm_out_1box.boxes[0].implicit.to_dict()
    assert len(d) == 30
    assert "o3" in d
    assert "oh" in d
    assert all(isinstance(v, float) for v in d.values())


def test_jvalues_to_dict(ctm_out_1box):
    d = ctm_out_1box.boxes[0].jvalues.to_dict()
    assert len(d) == 44
    assert "no2" in d
    assert "o3_o1d" in d


def test_long_lived_to_dict(ctm_out_1box):
    d = ctm_out_1box.boxes[0].long_lived.to_dict()
    assert len(d) == 18
    assert "o3" in d
    assert "ch4" in d


# ── LongLivedMixingRatios constructor and setters ─────────────────────────────

def test_long_lived_constructor():
    mr = LongLivedMixingRatios(o3=5e-6, ch4=1.8e-6, h2o=5e-3)
    assert mr.o3 == pytest.approx(5e-6)
    assert mr.ch4 == pytest.approx(1.8e-6)
    assert mr.h2o == pytest.approx(5e-3)
    assert mr.n2o == 0.0


def test_long_lived_default_zero():
    mr = LongLivedMixingRatios()
    assert mr.o3 == 0.0
    assert mr.brx == 0.0


def test_long_lived_setters():
    mr = LongLivedMixingRatios()
    mr.o3 = 4e-6
    mr.ch4 = 2e-6
    assert mr.o3 == pytest.approx(4e-6)
    assert mr.ch4 == pytest.approx(2e-6)


# ── DiurnBoxSpec / CtmBoxSpec constructors ────────────────────────────────────

def test_diurn_box_spec():
    b = DiurnBoxSpec(altitude_level=25, albedo=0.1, temp_offset_k=5.0)
    assert b.altitude_level == 25
    assert b.albedo == pytest.approx(0.1)
    assert b.temp_offset_k == pytest.approx(5.0)


def test_ctm_box_spec_defaults():
    b = CtmBoxSpec(altitude_level=10)
    assert b.altitude_level == 10
    assert b.albedo == 0.0
    assert b.temp_offset_k == 0.0


# ── DiurnConfig / CtmConfig constructors ──────────────────────────────────────

def test_diurn_config_defaults():
    cfg = DiurnConfig()
    assert cfg.latitude_deg == pytest.approx(0.0)
    assert cfg.julian_day == 120
    assert cfg.integration_days == 20
    assert cfg.bromine is False
    assert cfg.solar_flux_scale == pytest.approx(1.0)
    assert cfg.initial_mixing_ratios is None


def test_ctm_config_defaults():
    cfg = CtmConfig()
    assert cfg.latitude_deg == pytest.approx(60.0)
    assert cfg.julian_day == 75
    assert cfg.integration_days == 40


def test_ctm_config_setters():
    cfg = CtmConfig()
    cfg.latitude_deg = 45.0
    cfg.julian_day = 180
    assert cfg.latitude_deg == pytest.approx(45.0)
    assert cfg.julian_day == 180


# ── Diagnostics ───────────────────────────────────────────────────────────────

def test_diagnostics(ctm_out_1box):
    diag = ctm_out_1box.diagnostics
    assert isinstance(diag.raxloop, float)
    assert isinstance(diag.radcount, float)


# ── DIURN runs ───────────────────────────────────────────────────────────────

def test_diurn_basic(model):
    """DIURN mode works with embedded data: fort02.x provides initial conditions."""
    cfg = DiurnConfig(
        latitude_deg=0.0,
        julian_day=120,
        integration_days=20,
        boxes=[DiurnBoxSpec(altitude_level=20)],
    )
    out = model.run_diurn(cfg)
    assert len(out.boxes) == 1
    assert len(out.time_series) == 1
    snap = out.boxes[0]
    assert snap.altitude_km > 0
    assert snap.implicit.o3 > 0
    assert snap.implicit.oh > 0


def test_diurn_time_series(model):
    cfg = DiurnConfig(
        latitude_deg=0.0,
        julian_day=120,
        integration_days=20,
        boxes=[DiurnBoxSpec(altitude_level=20)],
    )
    out = model.run_diurn(cfg)
    ts = out.time_series[0]
    assert len(ts.steps) == 34   # standard 24-hour DIURN grid
    assert all(isinstance(step.time_hhmm, int) for step in ts.steps)


def test_diurn_species_grid(model):
    cfg = DiurnConfig(
        latitude_deg=0.0,
        julian_day=120,
        integration_days=20,
        boxes=[DiurnBoxSpec(altitude_level=15), DiurnBoxSpec(altitude_level=20)],
    )
    out = model.run_diurn(cfg)
    o3 = out.species_grid("o3")
    assert isinstance(o3, np.ndarray)
    assert o3.shape == (2, 34)
    assert np.all(o3 >= 0)


# ── PratmoModel repr ──────────────────────────────────────────────────────────

def test_pratmo_model_repr(model):
    assert repr(model) == "PratmoModel()"
