// Custom-atmosphere tests cover the normal physical DIURN API.  They are not
// compiled under `fortran-parity`, which intentionally disables the legacy
// reference executable's uninitialized photolysis multipliers.
#![cfg(not(feature = "fortran-parity"))]

use pratmo_core::api::{
    CustomAtmosphereProfile, DiurnBoxSpec, DiurnConfig, LongLivedMixingRatios,
    No2ConstrainedDiurnConfig, O3InputKind, PratmoModel,
};

#[test]
fn custom_atmosphere_diurn_runs() {
    let model = PratmoModel::with_defaults();
    let cfg = DiurnConfig {
        latitude_deg: 0.0,
        julian_day: 120,
        integration_days: 1,
        bromine: true,
        iodine: false,
        atmosphere: Some(CustomAtmosphereProfile {
            pressure_mb: vec![50.0],
            temperature_k: vec![220.0],
            altitude_km: Some(vec![21.5]),
            o3: vec![5.0e-6],
            o3_kind: O3InputKind::MixingRatio,
        }),
        ..Default::default()
    };

    let out = model.run_diurn(&cfg).expect("custom DIURN run failed");
    assert_eq!(out.boxes.len(), 1);
    assert_eq!(out.time_series.len(), 1);
    assert!((out.boxes[0].pressure_mb - 50.0).abs() < 1e-12);
    assert!((out.boxes[0].temperature_k - 220.0).abs() < 1e-12);
    assert!((out.boxes[0].altitude_km - 21.5).abs() < 1e-12);
    assert!(out.boxes[0].implicit.o3 > 0.0);
}

#[test]
fn no2_constrained_diurn_returns_scale() {
    let model = PratmoModel::with_defaults();
    let cfg = DiurnConfig {
        latitude_deg: 0.0,
        julian_day: 120,
        integration_days: 1,
        bromine: true,
        iodine: false,
        atmosphere: Some(CustomAtmosphereProfile {
            pressure_mb: vec![50.0],
            temperature_k: vec![220.0],
            altitude_km: Some(vec![21.5]),
            o3: vec![5.0e-6],
            o3_kind: O3InputKind::MixingRatio,
        }),
        ..Default::default()
    };
    let constrained = No2ConstrainedDiurnConfig {
        diurn: cfg,
        observed_no2_cm3: vec![1.0e6],
        target_hhmm: 630,
        iterations: 1,
    };

    let out = model
        .run_diurn_no2_constrained(&constrained)
        .expect("constrained DIURN run failed");
    assert_eq!(out.noy_scale.len(), 1);
    assert_eq!(out.modeled_no2_cm3.len(), 1);
    assert!(out.noy_scale[0].is_finite());
}

#[test]
fn custom_atmosphere_respects_box_altitude_level() {
    let model = PratmoModel::with_defaults();
    let cfg = DiurnConfig {
        latitude_deg: 0.0,
        julian_day: 120,
        integration_days: 1,
        boxes: vec![DiurnBoxSpec {
            altitude_level: 3,
            aerosol_surface_area_um2_cm3: 0.0,
            sea_salt_surface_area_um2_cm3: 0.0,
            temp_offset_k: 0.0,
        }],
        bromine: true,
        iodine: false,
        atmosphere: Some(CustomAtmosphereProfile {
            pressure_mb: vec![80.0, 50.0, 30.0],
            temperature_k: vec![225.0, 220.0, 215.0],
            altitude_km: Some(vec![18.0, 21.0, 24.0]),
            o3: vec![3.5e-6, 5.0e-6, 6.0e-6],
            o3_kind: O3InputKind::MixingRatio,
        }),
        ..Default::default()
    };

    let out = model.run_diurn(&cfg).expect("custom DIURN run failed");
    assert_eq!(out.boxes.len(), 1);
    assert!((out.boxes[0].pressure_mb - 30.0).abs() < 1e-12);
    assert!((out.boxes[0].temperature_k - 215.0).abs() < 1e-12);
    assert!((out.boxes[0].altitude_km - 24.0).abs() < 1e-12);
}

#[test]
fn custom_atmosphere_level_is_independent_of_box_slot() {
    let model = PratmoModel::with_defaults();
    let atmosphere = CustomAtmosphereProfile {
        pressure_mb: vec![80.0, 50.0, 30.0],
        temperature_k: vec![225.0, 220.0, 215.0],
        altitude_km: Some(vec![18.0, 21.0, 24.0]),
        o3: vec![3.5e-6, 5.0e-6, 6.0e-6],
        o3_kind: O3InputKind::MixingRatio,
    };
    let initial = LongLivedMixingRatios {
        o3: 6.0e-6,
        n2o: 3.0e-7,
        noy: 1.0e-8,
        ch4: 1.8e-6,
        co: 5.0e-8,
        clx: 2.0e-9,
        cf2cl2: 5.0e-10,
        cfcl3: 2.5e-10,
        ccl4: 1.0e-10,
        ch3cl: 5.0e-10,
        ch3ccl3: 1.0e-10,
        h2: 5.0e-7,
        h2o: 4.0e-6,
        nh3: 1.0e-12,
        c5h8: 0.0,
        brx: 1.0e-11,
        ch3br: 1.0e-11,
        ocs: 5.0e-10,
        iodx: 0.0,
    };

    let single = DiurnConfig {
        latitude_deg: 0.0,
        julian_day: 120,
        integration_days: 1,
        boxes: vec![DiurnBoxSpec {
            altitude_level: 3,
            aerosol_surface_area_um2_cm3: 0.0,
            sea_salt_surface_area_um2_cm3: 0.0,
            temp_offset_k: 0.0,
        }],
        bromine: true,
        iodine: false,
        atmosphere: Some(atmosphere.clone()),
        initial_mixing_ratios: Some(vec![initial.clone()]),
        ..Default::default()
    };
    let multi = DiurnConfig {
        boxes: vec![
            DiurnBoxSpec {
                altitude_level: 1,
                aerosol_surface_area_um2_cm3: 0.0,
                sea_salt_surface_area_um2_cm3: 0.0,
                temp_offset_k: 0.0,
            },
            DiurnBoxSpec {
                altitude_level: 2,
                aerosol_surface_area_um2_cm3: 0.0,
                sea_salt_surface_area_um2_cm3: 0.0,
                temp_offset_k: 0.0,
            },
            DiurnBoxSpec {
                altitude_level: 3,
                aerosol_surface_area_um2_cm3: 0.0,
                sea_salt_surface_area_um2_cm3: 0.0,
                temp_offset_k: 0.0,
            },
        ],
        initial_mixing_ratios: Some(vec![initial.clone(), initial.clone(), initial]),
        ..single.clone()
    };

    let single_out = model
        .run_diurn(&single)
        .expect("single custom DIURN run failed");
    let multi_out = model
        .run_diurn(&multi)
        .expect("multi custom DIURN run failed");
    let single_no2 = single_out.boxes[0].implicit.no2;
    let multi_no2 = multi_out.boxes[2].implicit.no2;
    let rel = (single_no2 - multi_no2).abs() / single_no2.max(1.0);

    assert!(
        rel < 0.10,
        "NO2 changed with box slot: single={single_no2:e}, multi={multi_no2:e}, rel={rel:e}"
    );
}
