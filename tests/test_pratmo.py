"""
Integration tests for the pratmo Python wrapper.

Both CTM and DIURN modes work with ``PratmoModel.with_defaults()`` using the
embedded science data (fort01.x, fort02.x, and the spectral/atmosphere files).
"""

from datetime import date
from pathlib import Path

import numpy as np
import pytest

from pratmo import (
    IMPLICIT_SPECIES_NAMES,
    JVALUE_NAMES,
    LONG_LIVED_NAMES,
    BoxSnapshot,
    Atmosphere,
    Box,
    ChemistryOptions,
    CtmBoxSpec,
    CtmConfig,
    CtmOptions,
    DiurnBoxSpec,
    DiurnConfig,
    DiurnalOptions,
    JValues,
    LongLivedMixingRatios,
    Model,
    PhotolysisOptions,
    PratmoModel,
    background_mixing_ratios,
    mixing_ratio_as,
    number_density,
    ppbv,
    ppmv,
    pressure,
    temperature,
)


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
CLIMATOLOGY_FIXTURE = REPOSITORY_ROOT / "pratmo-core/tests/fixtures/legacy_inputs"


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
    assert len(out) == 4
    assert len(out.boxes) == 4
    # Altitude should increase with box index (lower level → lower altitude)
    alts = [b.altitude_km for b in out.boxes]
    assert all(a > 0 for a in alts)


def test_ctm_coordinate_arrays(ctm_out_4box):
    assert ctm_out_4box.altitude_km.shape == (4,)
    assert ctm_out_4box.pressure_mb.shape == (4,)
    assert ctm_out_4box.temperature_k.shape == (4,)
    assert ctm_out_4box.air_density_cm3.shape == (4,)
    assert np.allclose(
        ctm_out_4box.altitude_km,
        [box.altitude_km for box in ctm_out_4box.boxes],
    )


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


def test_field_names_are_case_insensitive(ctm_out_4box):
    assert np.array_equal(
        ctm_out_4box.species_profile(" O3 "),
        ctm_out_4box.species_profile("o3"),
    )
    assert np.array_equal(
        ctm_out_4box.jvalue_profile("NO2"),
        ctm_out_4box.jvalue_profile("no2"),
    )


def test_species_profile_all_names(ctm_out_4box):
    """Verify every valid implicit species name can be queried."""
    assert len(IMPLICIT_SPECIES_NAMES) == 40
    for name in IMPLICIT_SPECIES_NAMES:
        arr = ctm_out_4box.species_profile(name)
        assert isinstance(arr, np.ndarray), f"species '{name}' did not return ndarray"
        assert arr[0] == pytest.approx(
            getattr(ctm_out_4box.boxes[0].implicit, name)
        ), f"species '{name}' returned the wrong field"


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
    assert len(JVALUE_NAMES) == 52
    for name in JVALUE_NAMES:
        arr = ctm_out_4box.jvalue_profile(name)
        assert isinstance(arr, np.ndarray), f"J-value '{name}' did not return ndarray"
        assert arr[0] == pytest.approx(
            getattr(ctm_out_4box.boxes[0].jvalues, name)
        ), f"J-value '{name}' returned the wrong field"


def test_jvalue_profile_invalid(ctm_out_1box):
    with pytest.raises(ValueError, match="Unknown"):
        ctm_out_1box.jvalue_profile("not_a_jvalue")


def test_long_lived_profile_all_names(ctm_out_4box):
    assert len(LONG_LIVED_NAMES) == 19
    for name in LONG_LIVED_NAMES:
        profile = ctm_out_4box.long_lived_profile(name)
        assert profile.shape == (4,)
        assert profile[0] == pytest.approx(
            getattr(ctm_out_4box.boxes[0].long_lived, name)
        )


def test_long_lived_profile_invalid(ctm_out_1box):
    with pytest.raises(ValueError, match="Unknown long-lived species"):
        ctm_out_1box.long_lived_profile("not_a_family")


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
    implicit = ctm_out_1box.boxes[0].implicit
    d = implicit.to_dict()
    assert set(d) == set(IMPLICIT_SPECIES_NAMES)
    for name in IMPLICIT_SPECIES_NAMES:
        assert d[name] == getattr(implicit, name), f"wrong to_dict value for {name}"


def test_jvalues_to_dict(ctm_out_1box):
    jvalues = ctm_out_1box.boxes[0].jvalues
    d = jvalues.to_dict()
    assert set(d) == set(JVALUE_NAMES)
    for name in JVALUE_NAMES:
        assert d[name] == getattr(jvalues, name), f"wrong to_dict value for {name}"


def test_long_lived_to_dict(ctm_out_1box):
    d = ctm_out_1box.boxes[0].long_lived.to_dict()
    assert set(d) == set(LONG_LIVED_NAMES)


# ── LongLivedMixingRatios constructor and setters ─────────────────────────────

def test_long_lived_constructor():
    mr = LongLivedMixingRatios(o3=5e-6, ch4=1.8e-6, h2o=5e-6)
    assert mr.o3 == pytest.approx(5e-6)
    assert mr.ch4 == pytest.approx(1.8e-6)
    assert mr.h2o == pytest.approx(5e-6)
    assert mr.n2o == pytest.approx(300e-9)


def test_long_lived_defaults_are_representative_and_non_singular():
    mr = LongLivedMixingRatios()
    assert mr.o3 == pytest.approx(5e-6)
    assert mr.n2o == pytest.approx(300e-9)
    assert mr.brx == pytest.approx(20e-12)
    assert mr.iodx == pytest.approx(1e-12)


def test_long_lived_setters():
    mr = LongLivedMixingRatios()
    mr.o3 = 4e-6
    mr.ch4 = 2e-6
    assert mr.o3 == pytest.approx(4e-6)
    assert mr.ch4 == pytest.approx(2e-6)


# ── DiurnBoxSpec / CtmBoxSpec constructors ────────────────────────────────────

def test_diurn_box_spec():
    b = DiurnBoxSpec(
        altitude_level=25,
        altitude_km=26.5,
        aerosol_surface_area_um2_cm3=0.1,
        sea_salt_surface_area_um2_cm3=0.02,
        temp_offset_k=5.0,
    )
    assert b.altitude_level == 25
    assert b.altitude_km == pytest.approx(26.5)
    assert b.aerosol_surface_area_um2_cm3 == pytest.approx(0.1)
    assert b.sea_salt_surface_area_um2_cm3 == pytest.approx(0.02)
    assert b.temp_offset_k == pytest.approx(5.0)


def test_ctm_box_spec_defaults():
    b = CtmBoxSpec(altitude_level=10)
    assert b.altitude_level == 10
    assert b.aerosol_surface_area_um2_cm3 == 0.0
    assert b.sea_salt_surface_area_um2_cm3 == 0.0
    assert b.temp_offset_k == 0.0


# ── DiurnConfig / CtmConfig constructors ──────────────────────────────────────

def test_diurn_config_defaults():
    cfg = DiurnConfig()
    assert cfg.latitude_deg == pytest.approx(0.0)
    assert cfg.julian_day == 120
    assert cfg.integration_days == 20
    assert len(cfg.boxes) == 1
    assert cfg.boxes[0].altitude_level == 15
    assert cfg.bromine is True
    assert cfg.iodine is False
    assert cfg.cpp_compatibility is False
    assert cfg.elapsed_time_hours is None
    assert cfg.solar_flux_scale == pytest.approx(1.0)
    assert cfg.surface_albedo == pytest.approx(0.20)
    assert cfg.heterogeneous_chemistry is True
    assert cfg.radiative_aerosol is False
    assert cfg.initial_mixing_ratios is None


def test_pratmo_climatology_uses_cpp_pressure_height_grid():
    from pratmo import PratmoClimatology

    profile = PratmoClimatology(CLIMATOLOGY_FIXTURE).sample(
        -85.0, date(2008, 1, 16), np.array([0.0, 8.0, 80.0])
    )
    assert profile.temperature_k == pytest.approx([256.4, 228.42518988, 182.79799212])
    assert profile.o3[:2] == pytest.approx([0.0, 0.0573016796e-6])
    assert profile.noy[:2] == pytest.approx([0.235e-9, 0.759941626e-9])
    assert profile.n2o[:2] == pytest.approx([317.615e-9, 315.160443e-9])
    assert profile.aerosol_surface_area_um2_cm3 == pytest.approx(
        [0.9399, 0.9399, 2.09534636e-8]
    )


def test_pratmo_climatology_matches_cpp_osiris_case():
    from pratmo import PratmoClimatology

    profile = PratmoClimatology(CLIMATOLOGY_FIXTURE).sample(
        30.466005325317383, date(2008, 7, 2), np.array([26.5])
    )
    assert profile.noy[0] == pytest.approx(10.65843794e-9)
    assert profile.n2o[0] == pytest.approx(192.1951233e-9)
    assert profile.aerosol_surface_area_um2_cm3[0] == pytest.approx(0.25008213)


def test_pratmo_climatology_builds_correlated_long_lived_inputs():
    from pratmo import PratmoClimatology

    profile = PratmoClimatology(CLIMATOLOGY_FIXTURE).sample(
        5.0, date(2008, 7, 2), np.array([26.5])
    )
    ratios = profile.initial_mixing_ratios(o3=np.array([5.0e-6]))[0]
    assert ratios.o3 == pytest.approx(5.0e-6)
    assert ratios.n2o == pytest.approx(profile.n2o[0])
    assert ratios.noy == pytest.approx(profile.noy[0])
    assert ratios.ch4 > 0.0
    assert ratios.h2o == pytest.approx(7.0e-6 - 2.0 * ratios.ch4)
    assert ratios.co == pytest.approx(3.0e-8)
    assert ratios.ch3br == pytest.approx(10.0e-12)
    assert ratios.brx > 2.0e-12
    assert ratios.iodx == pytest.approx(1.0e-12)


def test_ctm_config_defaults():
    cfg = CtmConfig()
    assert cfg.latitude_deg == pytest.approx(60.0)
    assert cfg.julian_day == 75
    assert cfg.integration_days == 40
    assert [box.altitude_level for box in cfg.boxes] == [10, 15, 20, 25]
    assert cfg.bromine is True
    assert cfg.iodine is False


def test_ctm_config_setters():
    cfg = CtmConfig()
    cfg.latitude_deg = 45.0
    cfg.julian_day = 180
    assert cfg.latitude_deg == pytest.approx(45.0)
    assert cfg.julian_day == 180


def test_low_level_defaults_are_runnable(model):
    ctm = model.run_ctm(CtmConfig())
    diurn = model.run_diurn(DiurnConfig())
    assert len(ctm) == 4
    assert len(diurn) == 1


# ── Diagnostics ───────────────────────────────────────────────────────────────

def test_diagnostics(ctm_out_1box):
    diag = ctm_out_1box.diagnostics
    assert isinstance(diag.raxloop, float)
    assert isinstance(diag.radcount, float)
    assert isinstance(diag.newraf_nonconvergence_count, int)
    assert isinstance(diag.rafday_nonconvergence_count, int)
    assert isinstance(diag.rafday_max_final_relative_correction, float)
    assert isinstance(diag.rafday_max_correction_iterations, int)
    assert diag.rafday_max_final_relative_correction >= 0.0
    assert diag.rafday_max_correction_iterations >= 0


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
    assert len(out) == 1
    assert len(ts) == 34
    assert len(ts.steps) == 34   # standard 24-hour DIURN grid
    assert all(isinstance(step.time_hhmm, int) for step in ts.steps)
    elapsed = [step.elapsed_seconds for step in ts.steps]
    assert elapsed[0] == 0.0
    assert elapsed[-1] == pytest.approx(24.0 * 3600.0)
    assert all(a < b for a, b in zip(elapsed, elapsed[1:]))
    assert ts.steps[0].time_hhmm == ts.steps[-1].time_hhmm == 1200
    assert np.array_equal(out.elapsed_seconds, elapsed)
    assert np.array_equal(out.time_hhmm, [step.time_hhmm for step in ts.steps])
    assert out.altitude_km.shape == (1,)
    assert out.pressure_mb.shape == (1,)
    assert out.temperature_k.shape == (1,)
    assert out.air_density_cm3.shape == (1,)
    assert out.species_profile("o3").shape == (1,)
    assert out.long_lived_profile("noy").shape == (1,)
    assert out.jvalue_profile("no2").shape == (1,)


def test_diurn_exact_elapsed_time_grid(model):
    hours = np.arange(0.0, 24.5, 0.5)
    out = model.run_diurn(
        DiurnConfig(
            latitude_deg=30.0,
            julian_day=184,
            integration_days=1,
            boxes=[DiurnBoxSpec(altitude_level=20)],
            elapsed_time_hours=hours.tolist(),
        )
    )
    assert np.array_equal(out.elapsed_seconds / 3600.0, hours)


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


def test_pratmo_model_constructor_uses_embedded_data():
    assert repr(PratmoModel()) == "PratmoModel()"


@pytest.mark.parametrize(
    ("config", "message"),
    [
        (CtmConfig(boxes=[]), "at least one box"),
        (
            CtmConfig(boxes=[CtmBoxSpec(altitude_level=42)]),
            "altitude_level",
        ),
        (
            DiurnConfig(boxes=[DiurnBoxSpec(altitude_level=20)], julian_day=0),
            "julian_day",
        ),
    ],
)
def test_invalid_configs_raise_value_error(model, config, message):
    run = model.run_ctm if isinstance(config, CtmConfig) else model.run_diurn
    with pytest.raises(ValueError, match=message):
        run(config)


# ── High-level Python interface and quantities ───────────────────────────────

def test_quantity_helpers_convert_scalars_and_arrays():
    assert ppmv(5.0) == pytest.approx(5.0e-6)
    assert ppbv(300.0) == pytest.approx(300.0e-9)
    assert pressure(5000.0, "Pa") == pytest.approx(50.0)
    assert temperature(-55.0, "degC") == pytest.approx(218.15)
    assert number_density(1.0e18, "m-3") == pytest.approx(1.0e12)
    assert np.array_equal(ppmv([1.0, 2.0]), np.array([1.0e-6, 2.0e-6]))
    assert np.array_equal(
        mixing_ratio_as(np.array([1.0e-9, 2.0e-9]), "ppbv"),
        np.array([1.0, 2.0]),
    )


def test_quantity_helpers_reject_unknown_units():
    with pytest.raises(ValueError, match="Unknown pressure unit"):
        pressure(1.0, "psi")


def test_atmosphere_converts_and_validates_units():
    atmosphere = Atmosphere(
        pressure=[8000.0, 5000.0, 3000.0],
        pressure_unit="Pa",
        temperature=[-48.15, -53.15, -58.15],
        temperature_unit="degC",
        altitude=[18000.0, 21000.0, 24000.0],
        altitude_unit="m",
        ozone=[3.5, 5.0, 6.0],
        ozone_unit="ppmv",
    )
    assert np.array_equal(atmosphere.pressure_mb, [80.0, 50.0, 30.0])
    assert np.array_equal(atmosphere.altitude_km, [18.0, 21.0, 24.0])
    assert np.allclose(atmosphere.ozone, [3.5e-6, 5.0e-6, 6.0e-6])
    assert atmosphere.ozone_kind == "mixing_ratio"


def test_atmosphere_rejects_bad_grid():
    with pytest.raises(ValueError, match="pressure must decrease"):
        Atmosphere(
            pressure=[50.0, 80.0],
            temperature=[220.0, 225.0],
            ozone=[5.0, 3.5],
        )


def test_atmosphere_warns_for_likely_unconverted_ozone():
    with pytest.warns(UserWarning, match="set ozone_unit='ppmv'"):
        Atmosphere(
            pressure=[80.0, 50.0],
            temperature=[225.0, 220.0],
            ozone=[3.5, 5.0],
            ozone_unit="fraction",
        )


def test_background_mixing_ratios_accept_readable_overrides():
    ratios = background_mixing_ratios(o3=ppmv(6.0), noy=ppbv(12.0))
    assert ratios.o3 == pytest.approx(6.0e-6)
    assert ratios.noy == pytest.approx(12.0e-9)
    assert ratios.iodx > 0.0


def test_high_level_ctm_smoke():
    output = Model().ctm(
        latitude=45.0,
        day="2026-03-20",
        boxes=[Box.at_level(20)],
        options=CtmOptions(integration_days=3),
    )
    assert output.altitude_km.shape == (1,)
    assert output.species_profile("o3")[0] > 0.0
    assert output.jvalue_profile("no2")[0] > 0.0


def test_high_level_custom_atmosphere_and_plots_smoke():
    atmosphere = Atmosphere(
        pressure=[80.0, 50.0, 30.0],
        temperature=[225.0, 220.0, 215.0],
        altitude=[18.0, 21.0, 24.0],
        ozone=[3.5, 5.0, 6.0],
        aerosol_surface_area=[0.3, 0.2, 0.1],
    )
    output = Model().diurnal(
        atmosphere=atmosphere,
        options=DiurnalOptions(integration_days=3),
    )
    from pratmo.plotting import plot_atmosphere, plot_diurnal, plot_profile

    assert len(output) == 3
    assert len(plot_atmosphere(atmosphere).data) == 4
    assert len(plot_diurnal(output, ["oh", "no2"]).data) == 2
    assert len(plot_profile(output, ["o3", "no2"]).data) == 2


def test_high_level_no2_constrained_smoke():
    atmosphere = Atmosphere(
        pressure=[50.0],
        temperature=[220.0],
        altitude=[21.0],
        ozone=[5.0],
    )
    result = Model().diurnal_no2_constrained(
        atmosphere=atmosphere,
        observed_no2_cm3=[1.0e8],
        target_hhmm=630,
        iterations=1,
        options=DiurnalOptions(integration_days=3),
    )
    assert len(result.output) == 1
    assert np.isfinite(result.noy_scale[0])
    assert result.modeled_no2_cm3[0] >= 0.0


def test_high_level_option_warnings():
    with pytest.warns(UserWarning, match="unusually reflective"):
        PhotolysisOptions(surface_albedo=0.95)
    with pytest.warns(UserWarning, match="experimental"):
        Model._warn_chemistry(ChemistryOptions(iodine=True))


def test_high_level_ctm_warns_when_latitude_is_binned():
    with pytest.warns(UserWarning, match="50° latitude grid point"):
        Model().ctm(
            latitude=52.1,
            boxes=[Box.at_level(20)],
            options=CtmOptions(integration_days=3),
        )
