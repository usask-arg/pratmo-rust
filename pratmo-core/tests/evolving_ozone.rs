//! Ozone evolution is an opt-in normal-mode extension, outside the fixed-input
//! Fortran differential policy.
#![cfg(not(feature = "fortran-parity"))]

use pratmo_core::api::{
    DiurnBoxSpec, DiurnConfig, DiurnOutput, ImplicitSpecies, JValues, LongLivedMixingRatios,
    PratmoModel,
};

macro_rules! values {
    ($object:expr; $($field:ident),+ $(,)?) => {
        [$( (stringify!($field), $object.$field) ),+]
    };
}

fn species_values(s: &ImplicitSpecies) -> [(&'static str, f64); 40] {
    values!(s;
        no, no2, no3, n2o5, hno3, h, oh, ho2, h2o2, o, o3, bro, br, hbr,
        hno2, hcl, cl, cl2, clo, clono2, hno4, hocl, brono2, hobr, h2co,
        ch3o2, ch3o2h, oclo, cl2o2, brcl, i, io, hoi, iono2, hi, oio,
        i2, i2o2, i2o3, i2o4,
    )
}

fn jvalues(j: &JValues) -> [(&'static str, f64); 52] {
    values!(j;
        no, o2, o3, o3_o1d, h2co_a, h2co_b, h2o2, rooh, no2, no3_x, no3_l,
        n2o5, hno2, hno3, hno4, clono2, cl2, hocl, oclo, cl2o2, clo, bro,
        brono2, hobr, n2o, cfc11, cfc12, cfc113, cfc114, cfc115, ccl4,
        ch3cl, ch3ccl3, ch3br, h1211, h1301, h2402, hcfc22, hcfc123,
        hcfc141b, chbr3, ch3i, cf3i, ocs, io, hoi, iono2, oio, i2, i2o2,
        i2o3, i2o4,
    )
}

fn long_lived_values(s: &LongLivedMixingRatios) -> [(&'static str, f64); 19] {
    values!(s;
        o3, n2o, noy, ch4, co, clx, cf2cl2, cfcl3, ccl4, ch3cl, ch3ccl3,
        h2, h2o, nh3, c5h8, brx, ch3br, ocs, iodx,
    )
}

fn config(levels: &[u8]) -> DiurnConfig {
    DiurnConfig {
        latitude_deg: 0.0,
        julian_day: 120, // 2026-04-30, matching the reported notebook case.
        bromine: true,
        iodine: false,
        boxes: levels
            .iter()
            .map(|&altitude_level| DiurnBoxSpec {
                altitude_level,
                altitude_km: None,
                aerosol_surface_area_um2_cm3: 0.0,
                sea_salt_surface_area_um2_cm3: 0.0,
                temp_offset_k: 0.0,
            })
            .collect(),
        ..Default::default()
    }
}

fn assert_converged(output: &DiurnOutput) {
    assert_eq!(output.diagnostics.newraf_nonconvergence_count, 0);
    assert_eq!(output.diagnostics.rafday_nonconvergence_count, 0);
    assert!(output.diagnostics.rafday_warnings.is_empty());
    for series in &output.time_series {
        let initial = series.steps.first().unwrap();
        let final_step = series.steps.last().unwrap();
        assert_eq!(initial.elapsed_seconds, 0.0);
        assert_eq!(final_step.elapsed_seconds, 86400.0);
        for step in &series.steps {
            for (name, value) in species_values(&step.implicit) {
                assert!(
                    value.is_finite() && value >= 0.0,
                    "box {} {name} at {} seconds: {value}",
                    series.box_index,
                    step.elapsed_seconds,
                );
            }
        }
        for ((name, start), (_, end)) in species_values(&initial.implicit)
            .into_iter()
            .zip(species_values(&final_step.implicit))
        {
            // Check every species, including the fast radicals omitted from
            // RAFDAY's slow-species system. The absolute scale protects only
            // trace concentrations below 0.01 molecules/cm3.
            let relative_change = (end - start).abs() / start.abs().max(1.0e-2);
            assert!(
                relative_change < 3.0e-5,
                "box {} {name}: noon-to-noon change {relative_change:e}",
                series.box_index,
            );
        }
    }
}

fn assert_outputs_identical(left: &DiurnOutput, right: &DiurnOutput) {
    assert_eq!(left.boxes.len(), right.boxes.len());
    assert_eq!(left.time_series.len(), right.time_series.len());
    for (left, right) in left.boxes.iter().zip(&right.boxes) {
        assert_eq!(left.box_index, right.box_index);
        assert_eq!(left.altitude_km, right.altitude_km);
        assert_eq!(left.pressure_mb, right.pressure_mb);
        assert_eq!(left.temperature_k, right.temperature_k);
        assert_eq!(left.air_density_cm3, right.air_density_cm3);
        assert_eq!(
            species_values(&left.implicit),
            species_values(&right.implicit)
        );
        assert_eq!(
            long_lived_values(&left.long_lived),
            long_lived_values(&right.long_lived)
        );
        assert_eq!(jvalues(&left.jvalues), jvalues(&right.jvalues));
    }
    for (left, right) in left.time_series.iter().zip(&right.time_series) {
        assert_eq!(left.box_index, right.box_index);
        assert_eq!(left.altitude_km, right.altitude_km);
        assert_eq!(left.pressure_mb, right.pressure_mb);
        assert_eq!(left.steps.len(), right.steps.len());
        for (left, right) in left.steps.iter().zip(&right.steps) {
            assert_eq!(left.elapsed_seconds, right.elapsed_seconds);
            assert_eq!(left.time_hhmm, right.time_hhmm);
            assert_eq!(
                species_values(&left.implicit),
                species_values(&right.implicit)
            );
        }
    }
    let left = &left.diagnostics;
    let right = &right.diagnostics;
    assert_eq!(left.raxloop, right.raxloop);
    assert_eq!(left.radcount, right.radcount);
    assert_eq!(
        left.newraf_nonconvergence_count,
        right.newraf_nonconvergence_count
    );
    assert_eq!(
        left.rafday_nonconvergence_count,
        right.rafday_nonconvergence_count
    );
    assert_eq!(left.rafday_warnings, right.rafday_warnings);
    assert_eq!(
        left.rafday_max_final_relative_correction,
        right.rafday_max_final_relative_correction
    );
    assert_eq!(
        left.rafday_max_correction_iterations,
        right.rafday_max_correction_iterations
    );
}

#[test]
fn notebook_level_40_closes_the_full_evolving_ozone_orbit() {
    let cfg = DiurnConfig {
        evolve_ozone: true,
        ..config(&[40])
    };
    let output = PratmoModel::with_defaults().run_diurn(&cfg).unwrap();
    assert_converged(&output);
    assert!((70.0..90.0).contains(&output.boxes[0].altitude_km));
    let ozone = output.species_grid(|species| species.o3);
    let min = ozone.iter().copied().fold(f64::INFINITY, f64::min);
    let max = ozone.iter().copied().fold(0.0_f64, f64::max);
    assert!(max > min * 1.1, "ozone should vary over the day");
    let snapshot = &output.boxes[0];
    assert_eq!(
        snapshot.long_lived.o3,
        snapshot.implicit.o3 / snapshot.air_density_cm3
    );
}

#[test]
fn ozone_evolution_is_opt_in_and_keeps_prescribed_photolysis() {
    let model = PratmoModel::with_defaults();
    let cfg = config(&[20]);
    assert!(!cfg.evolve_ozone);
    let default = model.run_diurn(&cfg).unwrap();
    let explicit_fixed = model
        .run_diurn(&DiurnConfig {
            evolve_ozone: false,
            ..cfg.clone()
        })
        .unwrap();
    assert_outputs_identical(&default, &explicit_fixed);
    let ozone = default.species_grid(|species| species.o3);
    assert!(ozone.iter().all(|value| *value == ozone[[0, 0]]));

    let evolving = model
        .run_diurn(&DiurnConfig {
            evolve_ozone: true,
            ..cfg.clone()
        })
        .unwrap();
    assert_converged(&evolving);
    assert_eq!(
        jvalues(&default.boxes[0].jvalues),
        jvalues(&evolving.boxes[0].jvalues)
    );
    let ozone = evolving.species_grid(|species| species.o3);
    assert!(ozone
        .iter()
        .any(|value| (*value - ozone[[0, 0]]).abs() > 1.0e-6 * ozone[[0, 0]]));

    // Reusing the same model must not retain the opt-in chemistry or mutate
    // the prescribed atmosphere used by subsequent default runs.
    let default_again = model.run_diurn(&cfg).unwrap();
    assert_outputs_identical(&default, &default_again);
}

#[test]
fn evolved_boxes_converge_identically_in_serial_and_parallel() {
    let model = PratmoModel::with_defaults();
    let cfg = DiurnConfig {
        evolve_ozone: true,
        ..config(&[12, 20, 30, 40])
    };
    let serial = model.run_diurn(&cfg).unwrap();
    let parallel = model
        .run_diurn(&DiurnConfig {
            parallel_boxes: true,
            ..cfg
        })
        .unwrap();
    assert_converged(&serial);
    assert_converged(&parallel);
    assert_outputs_identical(&serial, &parallel);
}
