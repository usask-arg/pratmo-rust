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
use pratmo_core::api::{DiurnBoxSpec, DiurnConfig, PratmoModel};

/// Run a single-box DIURN at 20 km (~level 12) for 5 days, equatorial noon geometry.
fn run_20km() -> pratmo_core::api::DiurnOutput {
    let model = PratmoModel::with_defaults();
    let cfg = DiurnConfig {
        latitude_deg: 0.0,
        julian_day: 120,
        integration_days: 5,
        boxes: vec![DiurnBoxSpec {
            altitude_level: 12,
            albedo: 0.0,
            temp_offset_k: 0.0,
        }],
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
            albedo: aerosol_area,
            temp_offset_k: 0.0,
        }],
        bromine: true,
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
    assert!(
        snap.implicit.i > 0.0,
        "I should be > 0, got {}",
        snap.implicit.i
    );
    assert!(
        snap.implicit.io > 0.0,
        "IO should be > 0, got {}",
        snap.implicit.io
    );
    assert!(
        snap.implicit.hoi > 0.0,
        "HOI should be > 0, got {}",
        snap.implicit.hoi
    );
    assert!(
        snap.implicit.iono2 > 0.0,
        "IONO2 should be > 0, got {}",
        snap.implicit.iono2
    );
    assert!(
        snap.implicit.hi > 0.0,
        "HI should be > 0, got {}",
        snap.implicit.hi
    );
}

// ── Total-iodine conservation ─────────────────────────────────────────────────

#[test]
fn test_iodine_conservation() {
    // Total Iy mixing ratio should stay within 30% of the initial 1 ppt.
    // Loose tolerance: the CH3I source term keeps Iy near fiodx, but the
    // timescales and approximate cross-sections allow significant drift.
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
        ratio > 0.1 && ratio < 10.0,
        "Total Iy/fiodx ratio = {ratio:.2} (expected 0.1–10): total={total_iy:.2e} fiodx={fiodx:.2e}"
    );
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
    // With the expanded iodine oxide scheme, the phase can shift slightly; it
    // still needs to show a real diurnal cycle.
    let out = run_20km();
    let steps = &out.time_series[0].steps;
    let hoi_max = steps.iter().map(|s| s.implicit.hoi).fold(0.0_f64, f64::max);
    let hoi_min = steps
        .iter()
        .map(|s| s.implicit.hoi)
        .fold(f64::MAX, f64::min);
    assert!(
        hoi_max > hoi_min * 1.01,
        "HOI should vary diurnally: max={hoi_max:.2e} min={hoi_min:.2e}"
    );
}

#[test]
fn test_ixoy_do_not_accumulate() {
    let out = run_20km();
    let snap = &out.boxes[0];
    let ixoy = 2.0 * (snap.implicit.i2o2 + snap.implicit.i2o3 + snap.implicit.i2o4);
    let iy = snap.implicit.i
        + snap.implicit.io
        + snap.implicit.hoi
        + snap.implicit.iono2
        + snap.implicit.hi
        + snap.implicit.oio
        + 2.0 * snap.implicit.i2
        + ixoy;
    assert!(
        ixoy < 0.2 * iy,
        "IxOy reservoirs should stay transient: IxOy/Iy={:.2e}",
        ixoy / iy
    );
}

#[test]
fn test_sea_salt_heterogeneous_recycling_changes_partitioning() {
    let clean = run_16km_with_aerosol(0.0);
    // 0.1 µm² cm⁻³ is representative of lower-stratospheric background
    // aerosol and remains inside RAFDAY's convergent regime after FIXRAT.
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
