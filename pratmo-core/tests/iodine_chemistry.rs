// Iodine/sea-salt tests assert the normal Rust chemistry, including active
// heterogeneous recycling; parity mode deliberately zeros the legacy DIURN
// heterogeneous-rate block for Fortran differential testing.
#![cfg(not(feature = "fortran-parity"))]

/// Iodine chemistry integration tests.
///
/// Uses the high-level `PratmoModel` API with embedded defaults (which include
/// iodine chemistry) to verify physical self-consistency of the Iy family:
///   I, IO, HOI, IONO2, HI, OIO, I2, I2O2, I2O3, I2O4
///
/// Tests are written as sanity checks on stratospheric photochemistry rather
/// than bit-exact regression tests, since the iodine chemistry is new and
/// cross-section data is approximate.
use pratmo_core::{
    api::{DiurnBoxSpec, DiurnConfig, PratmoModel},
    reader::{FortranReader, ModelReader},
    state::ModelState,
};

/// Run a single-box DIURN at 20 km (~level 12) for 5 days, equatorial noon geometry.
fn run_20km() -> pratmo_core::api::DiurnOutput {
    run_20km_with_iodine(true)
}

fn run_20km_with_iodine(iodine: bool) -> pratmo_core::api::DiurnOutput {
    let model = PratmoModel::with_defaults();
    let cfg = DiurnConfig {
        latitude_deg: 0.0,
        julian_day: 120,
        integration_days: 5,
        boxes: vec![DiurnBoxSpec {
            altitude_level: 12,
            aerosol_surface_area_um2_cm3: 0.0,
            sea_salt_surface_area_um2_cm3: 0.0,
            temp_offset_k: 0.0,
        }],
        iodine,
        ..Default::default()
    };
    model.run_diurn(&cfg).expect("DIURN run failed")
}

fn run_16km_with_aerosol(aerosol_area: f64) -> pratmo_core::api::DiurnOutput {
    let model = PratmoModel::with_defaults();
    let cfg = DiurnConfig {
        latitude_deg: 0.0,
        julian_day: 120,
        integration_days: 2,
        boxes: vec![DiurnBoxSpec {
            altitude_level: 9,
            aerosol_surface_area_um2_cm3: 0.0,
            sea_salt_surface_area_um2_cm3: aerosol_area,
            temp_offset_k: 0.0,
        }],
        bromine: false,
        ..Default::default()
    };
    model.run_diurn(&cfg).expect("DIURN run failed")
}

/// Find the time step closest to solar noon (max OH proxy) and night (min OH).
fn noon_and_night_idx(steps: &[pratmo_core::api::DiurnTimeStep]) -> (usize, usize) {
    let (noon, _) = steps
        .iter()
        .enumerate()
        .max_by(|(_, a), (_, b)| a.implicit.oh.partial_cmp(&b.implicit.oh).unwrap())
        .unwrap();
    let (night, _) = steps
        .iter()
        .enumerate()
        .min_by(|(_, a), (_, b)| a.implicit.oh.partial_cmp(&b.implicit.oh).unwrap())
        .unwrap();
    (noon, night)
}

// ── Basic presence tests ──────────────────────────────────────────────────────

#[test]
fn test_iodine_species_positive() {
    let out = run_20km();
    let snap = &out.boxes[0];
    for (name, value) in [
        ("I", snap.implicit.i),
        ("IO", snap.implicit.io),
        ("HOI", snap.implicit.hoi),
        ("IONO2", snap.implicit.iono2),
        ("HI", snap.implicit.hi),
        ("OIO", snap.implicit.oio),
        ("I2", snap.implicit.i2),
        ("I2O2", snap.implicit.i2o2),
        ("I2O3", snap.implicit.i2o3),
        ("I2O4", snap.implicit.i2o4),
    ] {
        assert!(value.is_finite() && value > 0.0, "{name}={value:e}");
    }
}

// ── Total-iodine conservation ─────────────────────────────────────────────────

#[test]
fn test_iodine_conservation() {
    // The gas-phase family constraint must preserve iodine atom equivalents.
    let out = run_20km();
    let snap = &out.boxes[0];
    let n = snap.air_density_cm3;
    let total_iy = (snap.implicit.i
        + snap.implicit.io
        + snap.implicit.hoi
        + snap.implicit.iono2
        + snap.implicit.hi
        + snap.implicit.oio
        + 2.0 * (snap.implicit.i2 + snap.implicit.i2o2 + snap.implicit.i2o3 + snap.implicit.i2o4))
        / n;
    let fiodx = snap.long_lived.iodx;
    let ratio = total_iy / fiodx;
    assert!(
        (ratio - 1.0).abs() < 1.0e-8,
        "Total Iy/fiodx ratio = {ratio:.12} (expected unity): total={total_iy:.2e} fiodx={fiodx:.2e}"
    );
}

#[test]
fn test_all_iodine_jvalues_are_loaded() {
    let out = run_20km();
    let j = &out.boxes[0].jvalues;
    for (name, value) in [
        ("IO", j.io),
        ("HOI", j.hoi),
        ("IONO2", j.iono2),
        ("OIO", j.oio),
        ("I2", j.i2),
        ("I2O2", j.i2o2),
        ("I2O3", j.i2o3),
        ("I2O4", j.i2o4),
    ] {
        assert!(value.is_finite() && value > 0.0, "J({name})={value}");
    }
}

#[test]
fn test_evaluated_iodine_rates_parse_at_fixed_width() {
    let mut state = ModelState::new();
    FortranReader::embedded().read_all(&mut state).unwrap();
    for (index, expected) in [
        (227usize, 1.60e-11), // HI + OH
        (228, 2.59e-12),      // IO + ClO -> I + OClO
        (250, 3.60e-16),      // IO + O3
        (262, 2.12e-12),      // IO + ClO atomic/ICl proxy
        (263, 1.20e-11),      // IO + BrO -> Br + OIO
    ] {
        let parsed = state.rk[[0, index]];
        assert!(
            (parsed / expected - 1.0).abs() < 1.0e-12,
            "rate {index} parsed as {parsed}, expected {expected}"
        );
    }
}

#[test]
fn test_iodine_off_is_a_converged_control() {
    let with_iodine = run_20km_with_iodine(true);
    let without_iodine = run_20km_with_iodine(false);
    for (label, diagnostics) in [
        ("iodine on", &with_iodine.diagnostics),
        ("iodine off", &without_iodine.diagnostics),
    ] {
        assert_eq!(
            diagnostics.newraf_nonconvergence_count, 0,
            "{label} recorded NEWRAF failures"
        );
        assert_eq!(
            diagnostics.rafday_nonconvergence_count, 0,
            "{label} recorded RAFDAY failures"
        );
    }
    assert!(without_iodine.boxes[0].long_lived.iodx < 1.0e-18);
    assert!(
        without_iodine.boxes[0].implicit.i < 1.0e-6 * with_iodine.boxes[0].implicit.i,
        "iodine-off numerical floor is not negligible"
    );
    for (name, iodine_on, iodine_off) in [
        (
            "OH",
            with_iodine.boxes[0].implicit.oh,
            without_iodine.boxes[0].implicit.oh,
        ),
        (
            "NO2",
            with_iodine.boxes[0].implicit.no2,
            without_iodine.boxes[0].implicit.no2,
        ),
        (
            "HNO3",
            with_iodine.boxes[0].implicit.hno3,
            without_iodine.boxes[0].implicit.hno3,
        ),
    ] {
        let relative_change = (iodine_on / iodine_off - 1.0).abs();
        assert!(
            relative_change < 0.01,
            "1 ppt iodine unexpectedly changes {name} by {relative_change:.3e}"
        );
    }
}

#[test]
fn test_realistic_lower_stratosphere_bry_rafday_converges() {
    let model = PratmoModel::with_defaults();
    let boxes = vec![DiurnBoxSpec {
        altitude_level: 8,
        aerosol_surface_area_um2_cm3: 0.0,
        sea_salt_surface_area_um2_cm3: 0.0,
        temp_offset_k: 0.0,
    }];

    // Keep the embedded background gases, but use the tropical lower-
    // stratospheric O3/NOy/Bry values used by the Saiz-Lopez comparison.
    let baseline = model
        .run_diurn(&DiurnConfig {
            latitude_deg: 0.0,
            julian_day: 120,
            integration_days: 1,
            boxes: boxes.clone(),
            ..Default::default()
        })
        .expect("baseline DIURN run failed");
    let mut mixing = baseline.boxes[0].long_lived.clone();
    mixing.o3 = 110.0e-9;
    mixing.noy = 50.0e-12;
    mixing.brx = 6.0e-12;
    mixing.iodx = 1.0e-12;

    let out = model
        .run_diurn(&DiurnConfig {
            latitude_deg: 0.0,
            julian_day: 120,
            integration_days: 1,
            boxes: boxes.clone(),
            bromine: true,
            iodine: true,
            initial_mixing_ratios: Some(vec![mixing.clone()]),
            ..Default::default()
        })
        .expect("realistic Bry DIURN run failed");

    assert_eq!(out.diagnostics.newraf_nonconvergence_count, 0);
    assert_eq!(
        out.diagnostics.rafday_nonconvergence_count, 0,
        "HNO3/BrONO2 family coupling prevented RAFDAY convergence"
    );
    assert!(
        out.diagnostics.rafday_max_final_relative_correction < 3.0e-5,
        "RAFDAY reported convergence with a final relative correction of {:.3e}",
        out.diagnostics.rafday_max_final_relative_correction
    );
    assert!(
        (1..=10).contains(&out.diagnostics.rafday_max_correction_iterations),
        "unexpected RAFDAY correction count: {}",
        out.diagnostics.rafday_max_correction_iterations
    );

    // The failure is in the legacy NOy/Bry slow-family solve, not in iodine.
    let without_iodine = model
        .run_diurn(&DiurnConfig {
            latitude_deg: 0.0,
            julian_day: 120,
            integration_days: 1,
            boxes,
            bromine: true,
            iodine: false,
            initial_mixing_ratios: Some(vec![mixing]),
            ..Default::default()
        })
        .expect("iodine-off realistic Bry DIURN run failed");
    assert_eq!(
        without_iodine.diagnostics.rafday_nonconvergence_count, 0,
        "NOy/Bry RAFDAY convergence must not depend on iodine"
    );

    let snap = &out.boxes[0];
    assert!(snap.implicit.hno3.is_finite() && snap.implicit.hno3 > 0.0);
    assert!(snap.implicit.brono2.is_finite() && snap.implicit.brono2 > 0.0);
    let hno3_mixing_ratio = snap.implicit.hno3 / snap.air_density_cm3;
    let brono2_mixing_ratio = snap.implicit.brono2 / snap.air_density_cm3;
    assert!(
        (hno3_mixing_ratio / 6.665216563e-12 - 1.0).abs() < 0.01,
        "HNO3 regression: {hno3_mixing_ratio:e}"
    );
    assert!(
        (brono2_mixing_ratio / 1.319634955e-13 - 1.0).abs() < 0.01,
        "BrONO2 regression: {brono2_mixing_ratio:e}"
    );

    // FIXRAT must preserve both linked families after the RAFDAY correction.
    let noy = (snap.implicit.no
        + snap.implicit.no2
        + snap.implicit.no3
        + 2.0 * snap.implicit.n2o5
        + snap.implicit.hno3
        + snap.implicit.hno2
        + snap.implicit.hno4
        + snap.implicit.clono2
        + snap.implicit.brono2)
        / snap.air_density_cm3;
    let bry = (snap.implicit.bro
        + snap.implicit.br
        + snap.implicit.hbr
        + snap.implicit.brono2
        + snap.implicit.hobr
        + snap.implicit.brcl)
        / snap.air_density_cm3;
    assert!((noy / 50.0e-12 - 1.0).abs() < 1.0e-8, "NOy={noy:e}");
    assert!((bry / 6.0e-12 - 1.0).abs() < 1.0e-8, "Bry={bry:e}");
}

// ── Diurnal cycle shape ───────────────────────────────────────────────────────

#[test]
fn test_io_dominates_i_at_noon() {
    // In the lower stratosphere at noon, I + O3 → IO is fast and IO accumulates.
    // IO/I >> 1 is the expected ratio during the day.
    let out = run_20km();
    let steps = &out.time_series[0].steps;
    let (noon_idx, _) = noon_and_night_idx(steps);
    let noon = &steps[noon_idx].implicit;
    assert!(
        noon.io > noon.i,
        "IO ({:.2e}) should exceed I ({:.2e}) at noon",
        noon.io,
        noon.i
    );
}

#[test]
fn test_iono2_diurnal_variation() {
    // IONO2 is produced by IO + NO2 + M and destroyed by photolysis.
    // At 20 km, IO formation at noon can make IONO2 peak in the day; at higher
    // altitudes stronger UV photolysis shifts the peak to night. Either way,
    // IONO2 must show a clear diurnal variation (max/min ratio > 1.05).
    let out = run_20km();
    let steps = &out.time_series[0].steps;
    let iono2_max = steps
        .iter()
        .map(|s| s.implicit.iono2)
        .fold(0.0_f64, f64::max);
    let iono2_min = steps
        .iter()
        .map(|s| s.implicit.iono2)
        .fold(f64::MAX, f64::min);
    assert!(
        iono2_max > iono2_min * 1.05,
        "IONO2 should vary diurnally: max={iono2_max:.2e} min={iono2_min:.2e}"
    );
}

#[test]
fn test_hoi_peaks_in_daytime() {
    // HOI is produced by IO + HO2 (daytime reaction) and is photolysed to I + OH.
    // In this equatorial 20 km case, both HOI and OH peak at local noon.
    let out = run_20km();
    let steps = &out.time_series[0].steps;
    let (noon_idx, _) = noon_and_night_idx(steps);
    let hoi_peak_idx = steps
        .iter()
        .enumerate()
        .max_by(|(_, a), (_, b)| a.implicit.hoi.partial_cmp(&b.implicit.hoi).unwrap())
        .map(|(index, _)| index)
        .unwrap();
    let hoi_max = steps[hoi_peak_idx].implicit.hoi;
    let hoi_min = steps
        .iter()
        .map(|s| s.implicit.hoi)
        .fold(f64::MAX, f64::min);
    assert_eq!(
        hoi_peak_idx, noon_idx,
        "HOI peak moved away from the OH-defined local noon"
    );
    assert!(
        hoi_max > hoi_min * 1.01,
        "HOI should vary diurnally: max={hoi_max:.2e} min={hoi_min:.2e}"
    );
}

#[test]
fn test_ixoy_transient_fraction_regression() {
    let out = run_20km();
    let maximum_fraction = out.time_series[0]
        .steps
        .iter()
        .map(|step| {
            let s = &step.implicit;
            let ixoy = 2.0 * (s.i2o2 + s.i2o3 + s.i2o4);
            let iy = s.i + s.io + s.hoi + s.iono2 + s.hi + s.oio + 2.0 * s.i2 + ixoy;
            ixoy / iy
        })
        .fold(0.0_f64, f64::max);
    assert!(
        maximum_fraction > 5.0e-5 && maximum_fraction < 2.0e-4,
        "unexpected peak IxOy/Iy fraction: {maximum_fraction:.3e}"
    );
}

#[test]
fn test_sea_salt_heterogeneous_recycling_changes_partitioning() {
    let clean = run_16km_with_aerosol(0.0);
    // Use the explicit sea-salt field; sulfate/other aerosol is independently
    // covered by the SETUPR unit test and must not activate these two rates.
    let salty = run_16km_with_aerosol(0.1);

    let clean_snap = &clean.boxes[0].implicit;
    let salty_snap = &salty.boxes[0].implicit;

    assert!(
        salty_snap.i > clean_snap.i,
        "Sea-salt iodine recycling should increase atomic I: clean={:.2e} salty={:.2e}",
        clean_snap.i,
        salty_snap.i
    );
}

// ── Non-regression: iodine should not blow up the rest of chemistry ───────────

#[test]
fn test_o3_reasonable_at_20km() {
    // O3 at 20 km should be in the range 2–10 ppm.
    // Iodine at 1 ppt contributes <0.1% to O3 loss, so this is effectively
    // a sanity check that the run converged normally.
    let out = run_20km();
    let snap = &out.boxes[0];
    let o3_mr = snap.implicit.o3 / snap.air_density_cm3;
    assert!(
        o3_mr > 1e-6 && o3_mr < 1e-4,
        "O3 mixing ratio at 20 km = {o3_mr:.2e}, expected 1e-6 to 1e-4"
    );
}

#[test]
fn test_oh_reasonable_at_20km() {
    // Daily-mean OH at 20 km should be in the range 1e-14 to 1e-11 (mixing ratio).
    let out = run_20km();
    let snap = &out.boxes[0];
    let oh_mr = snap.implicit.oh / snap.air_density_cm3;
    assert!(
        oh_mr > 1e-15 && oh_mr < 1e-10,
        "OH mixing ratio at 20 km = {oh_mr:.2e}, expected 1e-15 to 1e-10"
    );
}
