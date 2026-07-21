//! Policy smoke tests for the opt-in `fortran-parity` feature.  The assertions
//! below intentionally check both sides of the policy: parity mode reproduces
//! the executable's GMU0/280-K, SSF, heterogeneous-rate, and R177 quirks,
//! while normal mode preserves the corrected/input-driven behavior.

mod common;

use common::input_dir;

use pratmo_core::{
    api::{DiurnBoxSpec, DiurnConfig, PratmoModel},
    chemistry::{chempl, setupr},
    diurnal::diurn,
    reader::{FortranReader, ModelReader},
    state::ModelState,
};

#[cfg(feature = "fortran-parity")]
use pratmo_core::ctm::ctmlfq;

#[test]
fn original_40_species_spectral_file_is_supported() {
    let dir = input_dir();
    let mut state = ModelState::new();
    let mut reader = FortranReader::new(dir);
    reader
        .read_spectral_data(&mut state)
        .expect("read_spectral_data failed");

    assert_eq!(
        state.njval, 44,
        "40 tabulated species plus four base J-values"
    );
}

#[test]
fn original_odf_isr_annotations_are_parsed() {
    let dir = input_dir();
    let mut state = ModelState::new();
    let mut reader = FortranReader::new(dir);
    reader
        .read_spectral_data(&mut state)
        .expect("read_spectral_data failed");

    assert!(
        state.isr.iter().any(|&band| band > 0),
        "ODF ISR bands must survive the trailing wavelength annotation"
    );
}

#[test]
fn read_all_uses_the_selected_legacy_policy() {
    let dir = input_dir();
    let mut state = ModelState::new();
    let mut reader = FortranReader::new(dir);
    reader.read_all(&mut state).expect("read_all failed");

    #[cfg(feature = "fortran-parity")]
    {
        assert_eq!(state.gmu0, -0.14, "bread.f hard-coded GMU0 override");
        assert!(
            state.t[..state.nc]
                .iter()
                .all(|&temperature| temperature == 280.0),
            "bread.f calibration block forces a flat 280 K profile"
        );
    }

    #[cfg(not(feature = "fortran-parity"))]
    {
        assert_eq!(state.gmu0, -0.12, "normal builds honor the fort01 value");
        assert!(
            state.t[..state.nc]
                .iter()
                .any(|&temperature| temperature != 280.0),
            "normal builds retain the atmospheric temperature profile"
        );
    }
}

#[test]
fn diurn_applies_the_selected_ssf_policy() {
    // An empty driver run is enough to exercise DIURN's policy setup without
    // coupling this test to a particular atmosphere or solver trajectory.
    let mut state = ModelState::new();
    diurn(&mut state).expect("empty DIURN policy run failed");

    #[cfg(feature = "fortran-parity")]
    assert!(state.ssf.iter().all(|&value| value == 0.0));

    #[cfg(not(feature = "fortran-parity"))]
    assert!(state.ssf.iter().all(|&value| value == 1.0));
}

#[test]
fn hno3_r177_policy_is_explicit() {
    // Set up a collision-free synthetic N1..N30 map.  With only R177 present,
    // CHEMPL isolates the exact parity switch without requiring a full run.
    let mut state = ModelState::new();
    for (slot, value) in state.n.iter_mut().take(30).enumerate() {
        *value = slot + 1;
    }
    state.ntot = 30;
    state.r.fill(0.0);
    state.r[176] = 1.0; // Fortran R177: BrONO2 + H2O -> HOBr + HNO3
    chempl(&mut state);

    let hno3_slot = state.n[4] - 1;
    #[cfg(feature = "fortran-parity")]
    assert_eq!(state.rp[hno3_slot], 0.0);

    #[cfg(not(feature = "fortran-parity"))]
    assert_eq!(state.rp[hno3_slot], 1.0);
}

#[test]
fn public_diurn_smoke_respects_photolysis_policy() {
    let model = PratmoModel::with_defaults();
    let config = DiurnConfig {
        latitude_deg: 0.0,
        julian_day: 120,
        integration_days: 1,
        boxes: vec![DiurnBoxSpec {
            altitude_level: 12,
            altitude_km: None,
            aerosol_surface_area_um2_cm3: 0.0,
            sea_salt_surface_area_um2_cm3: 0.0,
            temp_offset_k: 0.0,
        }],
        iodine: false,
        ..Default::default()
    };
    let output = model
        .run_diurn(&config)
        .expect("public DIURN smoke run failed");
    assert_eq!(output.boxes.len(), 1);
    assert!(!output.time_series[0].steps.is_empty());
    let steps = &output.time_series[0].steps;
    assert_eq!(steps.first().unwrap().elapsed_seconds, 0.0);
    assert!((steps.last().unwrap().elapsed_seconds - 86_400.0).abs() < 1.0);
    assert!(steps
        .windows(2)
        .all(|pair| pair[0].elapsed_seconds < pair[1].elapsed_seconds));
    assert_eq!(steps.first().unwrap().time_hhmm, 1200);
    assert_eq!(steps.last().unwrap().time_hhmm, 1200);

    #[cfg(feature = "fortran-parity")]
    assert_eq!(output.boxes[0].jvalues.o3_o1d, 0.0);

    #[cfg(not(feature = "fortran-parity"))]
    assert!(output.boxes[0].jvalues.o3_o1d > 0.0);
}

#[test]
fn representative_latitude_season_matrix_stays_finite() {
    // Exercise the public DIURN path across tropics, mid-latitudes, and polar
    // cases.  This is an invariant matrix rather than a Fortran golden: the
    // legacy executable has no portable fixture for these API configurations.
    let model = PratmoModel::with_defaults();
    let cases = [
        (0.0, 15),
        (0.0, 172),
        (60.0, 75),
        (-60.0, 289),
        (89.0, 355),
        (-89.0, 172),
    ];
    for (latitude_deg, julian_day) in cases {
        let config = DiurnConfig {
            latitude_deg,
            julian_day,
            integration_days: 1,
            boxes: vec![DiurnBoxSpec {
                altitude_level: 20,
                altitude_km: None,
                aerosol_surface_area_um2_cm3: 0.0,
                sea_salt_surface_area_um2_cm3: 0.0,
                temp_offset_k: 0.0,
            }],
            iodine: false,
            ..Default::default()
        };
        let output = model.run_diurn(&config).unwrap_or_else(|error| {
            panic!("DIURN failed at {latitude_deg}°, day {julian_day}: {error}")
        });
        assert_eq!(output.boxes.len(), 1);
        let snapshot = &output.boxes[0];
        assert!(snapshot.altitude_km.is_finite());
        assert!(snapshot.pressure_mb.is_finite());
        assert!(snapshot.temperature_k.is_finite());
        assert!(snapshot.air_density_cm3.is_finite());
        assert!(snapshot.implicit.o3.is_finite());
        assert!(snapshot.implicit.oh.is_finite());
        assert!(snapshot.long_lived.noy.is_finite());
        assert!(snapshot.long_lived.n2o.is_finite());
        assert!(snapshot.jvalues.o3_o1d.is_finite());
        assert!(output.time_series[0]
            .steps
            .iter()
            .all(|step| step.implicit.o3.is_finite() && step.implicit.oh.is_finite()));
    }
}

#[cfg(feature = "fortran-parity")]
#[test]
fn parity_ctm_smoke_completes_reference_grid() {
    // Keep one end-to-end CTM guard in the feature build.  The detailed
    // gfortran golden comparison remains documented separately because its
    // atmospheric policy is intentionally different from normal Rust mode.
    let dir = input_dir();
    let mut state = ModelState::new();
    state.cinpdir = dir.to_string_lossy().into_owned();
    let mut reader = FortranReader::new(dir);
    reader.read_all(&mut state).expect("read_all failed");
    ctmlfq(&mut state).expect("parity CTM run failed");

    assert_eq!(state.ntimdo, 34);
    assert!((33_000.0..=33_100.0).contains(&state.radcount));
    let nbox = state.nbox.min(state.z.len());
    assert!(nbox > 0);
    assert!(state.z[..nbox].iter().all(|value| value.is_finite()));
    assert!(state.dno[..nbox].iter().all(|value| value.is_finite()));
}

#[test]
fn diurn_parity_zeros_legacy_heterogeneous_rates() {
    let dir = input_dir();
    let mut state = ModelState::new();
    let mut reader = FortranReader::new(dir);
    reader.read_all(&mut state).expect("read_all failed");

    state.nd216 = 0;
    state.ibox = 0;
    state.ialt = state.nboxdo[0].unsigned_abs().saturating_sub(1) as usize;
    setupr(&mut state);

    #[cfg(feature = "fortran-parity")]
    assert!(state.ratek[169..177].iter().all(|&rate| rate == 0.0));

    #[cfg(not(feature = "fortran-parity"))]
    assert!(state.ratek[169..177].iter().any(|&rate| rate > 0.0));
}
