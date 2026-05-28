use pratmo_core::api::{
    CustomAtmosphereProfile, DiurnConfig, No2ConstrainedDiurnConfig, O3InputKind, PratmoModel,
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
