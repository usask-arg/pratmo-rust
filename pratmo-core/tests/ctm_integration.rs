/// Full-model integration tests: run PRATMO CTM mode and compare every key
/// species at every active altitude level against the validated gfortran
/// reference output (60°N, March 16, 25 boxes, 40 days).
///
/// Column layout in boxout.dat (fixed-width 13 chars per field):
///   0=z(km), 4=O3, 5=N2O, 6=CH4, 7=H2O, 8=NOy
///   113=HNO3 it0, 215=N2O5 it0, 249=OH it0, 283=HO2 it0

use pratmo_core::{ctm::ctmlfq, reader::{FortranReader, ModelReader}, state::ModelState};
use std::path::{Path, PathBuf};

fn input_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent().unwrap()
        .join("fortran")
}

fn run_ctm() -> Box<ModelState> {
    let dir = input_dir();
    let mut s = ModelState::new();
    s.cinpdir = dir.to_string_lossy().into_owned();
    let mut reader = FortranReader::new(&dir);
    reader.read_all(&mut s).expect("read_all failed");
    ctmlfq(&mut s).expect("ctmlfq failed");
    s
}

/// Parse boxout.dat, return only the 25 active rows (those where NOy ≠ −9).
/// Each row is a fixed-width 13-char per field array of the 825 data columns.
fn active_rows(boxout_path: &Path) -> Vec<Vec<f64>> {
    let text = std::fs::read_to_string(boxout_path).expect("cannot read boxout.dat");
    let lines: Vec<&str> = text.lines().collect();

    // Find the "---" separator
    let sep = lines.iter().position(|l| l.contains("---") && l.chars().filter(|&c| c == '-').count() > 50)
        .expect("separator not found");

    // Each data line has 825 × 13-char fields (no leading offset in data rows)
    lines[sep + 1..].iter().filter_map(|line| {
        if line.len() < 9 * 13 { return None; }
        // NOy is field 8; skip fill rows where it is −9
        let noy = &line[8 * 13..9 * 13];
        if noy.contains("-9.") || noy.trim().is_empty() { return None; }
        let n = line.len() / 13;
        let row: Vec<f64> = (0..n.min(825))
            .map(|i| line[i*13..(i+1)*13].trim().parse().unwrap_or(0.0))
            .collect();
        if row.len() < 9 { return None; }
        Some(row)
    }).collect()
}

fn relerr(got: f64, exp: f64) -> f64 {
    (got - exp).abs() / exp.abs().max(1e-300)
}

fn assert_close(got: f64, exp: f64, tol: f64, label: &str) {
    assert!(
        relerr(got, exp) < tol,
        "{label}: got {got:.5E} expected {exp:.5E} relerr={:.2e}",
        relerr(got, exp)
    );
}

// ── Reference data ────────────────────────────────────────────────────────────
// Values from validated gfortran reference run (60°N, March 16, 40 days).
// Columns: z_km, O3, N2O, CH4, H2O, NOy, HNO3_it0, N2O5_it0, OH_it0, HO2_it0

#[rustfmt::skip]
const REF: &[(f64, f64, f64, f64, f64, f64, f64, f64, f64, f64)] = &[
    (57.87, 1.1805e-6, 1.6185e-9, 2.8775e-7, 6.4245e-6, 3.3432e-9, 1.4204e-13, 8.2148e-21, 1.0703e-9, 4.6573e-10),
    (55.77, 1.4870e-6, 1.6703e-9, 2.8917e-7, 6.4217e-6, 3.6170e-9, 3.3120e-13, 9.7028e-20, 1.0046e-9, 4.3820e-10),
    (53.63, 1.8377e-6, 1.6865e-9, 2.8962e-7, 6.4208e-6, 3.9143e-9, 6.8182e-13, 1.1155e-18, 9.3342e-10, 4.1002e-10),
    (51.46, 2.2178e-6, 1.6910e-9, 2.8974e-7, 6.4205e-6, 4.2120e-9, 1.2328e-12, 1.1903e-17, 8.4447e-10, 3.7635e-10),
    (49.27, 2.6530e-6, 1.7315e-9, 2.9085e-7, 6.4183e-6, 4.5092e-9, 2.2182e-12, 1.1734e-16, 7.3490e-10, 3.3531e-10),
    (47.08, 3.2447e-6, 2.1725e-9, 3.0276e-7, 6.3945e-6, 4.9190e-9, 5.0820e-12, 1.2971e-15, 6.1803e-10, 2.9080e-10),
    (44.90, 4.0275e-6, 2.6132e-9, 3.1442e-7, 6.3712e-6, 5.3470e-9, 1.4978e-11, 4.9627e-13, 4.8901e-10, 2.4274e-10),
    (42.75, 5.0305e-6, 4.2325e-9, 3.5524e-7, 6.2895e-6, 6.0995e-9, 5.0778e-11, 2.8702e-12, 3.5352e-10, 1.9643e-10),
    (40.63, 6.0482e-6, 7.5040e-9, 4.2849e-7, 6.1430e-6, 6.9070e-9, 1.4637e-10, 1.3031e-11, 2.2761e-10, 1.5903e-10),
    (38.55, 6.7967e-6, 1.0776e-8, 4.9050e-7, 6.0190e-6, 8.0557e-9, 3.6063e-10, 4.5326e-11, 1.3552e-10, 1.3278e-10),
    (36.51, 7.1185e-6, 1.8327e-8, 5.9783e-7, 5.8043e-6, 9.2612e-9, 7.5956e-10, 1.1731e-10, 7.8260e-11, 1.1055e-10),
    (34.51, 7.1062e-6, 2.7303e-8, 6.7879e-7, 5.6424e-6, 1.0664e-8, 1.4874e-9,  2.4505e-10, 4.5761e-11, 8.7686e-11),
    (32.55, 6.9042e-6, 3.7290e-8, 7.4055e-7, 5.5189e-6, 1.2098e-8, 2.6801e-9,  4.2850e-10, 2.7274e-11, 6.5483e-11),
    (30.63, 6.5462e-6, 5.8262e-8, 8.5113e-7, 5.2977e-6, 1.3131e-8, 4.1718e-9,  6.1434e-10, 1.6034e-11, 4.6777e-11),
    (28.73, 6.1027e-6, 7.9234e-8, 9.3488e-7, 5.1302e-6, 1.4098e-8, 5.7942e-9,  7.8026e-10, 9.1443e-12, 3.2096e-11),
    (26.85, 5.6297e-6, 1.0612e-7, 1.0268e-6, 4.9464e-6, 1.3858e-8, 6.7576e-9,  8.2894e-10, 4.9639e-12, 2.1764e-11),
    (24.98, 5.1382e-6, 1.4126e-7, 1.1235e-6, 4.7531e-6, 1.3415e-8, 7.2916e-9,  8.4897e-10, 2.6231e-12, 1.3858e-11),
    (23.13, 4.5720e-6, 1.7639e-7, 1.2056e-6, 4.5888e-6, 1.1673e-8, 6.7620e-9,  7.7868e-10, 1.3866e-12, 8.9650e-12),
    (21.28, 3.9165e-6, 2.0372e-7, 1.2595e-6, 4.4810e-6, 9.7163e-9, 5.8582e-9,  6.7788e-10, 7.4823e-13, 5.7320e-12),
    (19.43, 3.0925e-6, 2.2846e-7, 1.3675e-6, 4.2650e-6, 7.6195e-9, 4.7442e-9,  5.4358e-10, 4.1627e-13, 3.5576e-12),
    (17.59, 2.1792e-6, 2.5237e-7, 1.4752e-6, 4.0496e-6, 5.5620e-9, 3.5770e-9,  3.8775e-10, 2.4340e-13, 2.1494e-12),
    (15.73, 1.3365e-6, 2.6723e-7, 1.5422e-6, 3.9157e-6, 3.7890e-9, 2.5139e-9,  2.4160e-10, 1.5261e-13, 1.2793e-12),
    (13.88, 8.4050e-7, 2.8209e-7, 1.6040e-6, 3.7920e-6, 2.5147e-9, 1.7173e-9,  1.4541e-10, 1.0205e-13, 8.6810e-13),
    (12.03, 5.7350e-7, 2.9502e-7, 1.6610e-6, 3.6780e-6, 1.7133e-9, 1.2005e-9,  9.1247e-11, 6.9293e-14, 6.5035e-13),
    (10.18, 2.8475e-7, 3.0625e-7, 1.7105e-6, 3.5789e-6, 1.1057e-9, 7.9517e-10, 4.7214e-11, 5.2380e-14, 4.7270e-13),
];

// ── Tests ─────────────────────────────────────────────────────────────────────

#[test]
fn test_radcount() {
    let s = run_ctm();
    // 40 days × 25 boxes × 33 Newton-Raphson calls per time step
    assert_eq!(s.radcount as i64, 33000, "RADCOUNT mismatch");
}

#[test]
fn test_ntimdo() {
    let s = run_ctm();
    assert_eq!(s.ntimdo, 34, "NTIMDO should be 34 for standard diurnal run");
}

#[test]
fn test_active_box_count() {
    let boxout = input_dir().join("boxout.dat");
    let rows = active_rows(&boxout);
    assert_eq!(rows.len(), 25, "Expected 25 active altitude boxes");
}

#[test]
fn test_long_lived_species() {
    let boxout = input_dir().join("boxout.dat");
    let rows = active_rows(&boxout);
    assert_eq!(rows.len(), REF.len());

    let tol = 5e-4; // 0.05% — tight enough to catch any 4-sig-fig regression
    for (i, (row, &(z_exp, o3, n2o, ch4, h2o, noy, _, _, _, _))) in rows.iter().zip(REF).enumerate() {
        let label = format!("box {i} z≈{z_exp:.1} km");
        assert_close(row[0],  z_exp, 1e-3, &format!("{label} z"));
        assert_close(row[4],  o3,   tol, &format!("{label} O3"));
        assert_close(row[5],  n2o,  tol, &format!("{label} N2O"));
        assert_close(row[6],  ch4,  tol, &format!("{label} CH4"));
        assert_close(row[7],  h2o,  tol, &format!("{label} H2O"));
        assert_close(row[8],  noy,  tol, &format!("{label} NOy"));
    }
}

#[test]
fn test_hno3_profile() {
    let boxout = input_dir().join("boxout.dat");
    let rows = active_rows(&boxout);
    let tol = 5e-4;
    for (i, (row, &(z, _, _, _, _, _, hno3, _, _, _))) in rows.iter().zip(REF).enumerate() {
        assert_close(row[113], hno3, tol, &format!("box {i} z≈{z:.1} km HNO3"));
    }
}

#[test]
fn test_n2o5_profile() {
    let boxout = input_dir().join("boxout.dat");
    let rows = active_rows(&boxout);
    let tol = 5e-4;
    for (i, (row, &(z, _, _, _, _, _, _, n2o5, _, _))) in rows.iter().zip(REF).enumerate() {
        assert_close(row[215], n2o5, tol, &format!("box {i} z≈{z:.1} km N2O5"));
    }
}

#[test]
fn test_oh_profile() {
    let boxout = input_dir().join("boxout.dat");
    let rows = active_rows(&boxout);
    let tol = 5e-4;
    for (i, (row, &(z, _, _, _, _, _, _, _, oh, _))) in rows.iter().zip(REF).enumerate() {
        assert_close(row[249], oh, tol, &format!("box {i} z≈{z:.1} km OH"));
    }
}

#[test]
fn test_ho2_profile() {
    let boxout = input_dir().join("boxout.dat");
    let rows = active_rows(&boxout);
    let tol = 5e-4;
    for (i, (row, &(z, _, _, _, _, _, _, _, _, ho2))) in rows.iter().zip(REF).enumerate() {
        assert_close(row[283], ho2, tol, &format!("box {i} z≈{z:.1} km HO2"));
    }
}

#[test]
fn test_oh_diurnal_at_41km() {
    // Verify the full 34-step OH diurnal cycle at the ~40.6 km box.
    // Reference values from gfortran (col 249..282 in the active row for 40.6 km).
    #[rustfmt::skip]
    const OH_41KM: [f64; 34] = [
        2.2761e-10, 2.2769e-10, 2.2270e-10, 2.0752e-10, 1.8188e-10,
        1.4692e-10, 1.1383e-10, 7.4250e-11, 4.5669e-11, 2.9097e-11,
        1.3451e-11, 5.8408e-12, 2.2483e-12, 1.2517e-12, 6.9265e-13,
        4.3048e-13, 1.4228e-13, 5.1203e-14, 2.0290e-14, 9.2403e-15,
        8.1687e-15, 6.0077e-15, 2.3179e-14, 3.1041e-13, 5.2476e-12,
        1.3350e-11, 4.5124e-11, 8.9079e-11, 1.3334e-10, 1.7379e-10,
        2.0225e-10, 2.2045e-10, 2.2696e-10, 2.2761e-10,
    ];

    let boxout = input_dir().join("boxout.dat");
    let rows = active_rows(&boxout);
    // Find the ~40.6 km box (index 8 in REF)
    let row = rows.iter()
        .find(|r| (r[0] - 40.63).abs() < 0.5)
        .expect("40.6 km box not found");

    let tol = 2e-3; // 0.2% — slightly looser for diurnal cycle comparison
    for (it, &exp) in OH_41KM.iter().enumerate() {
        if exp < 1e-16 { continue; } // skip near-zero night steps
        assert_close(row[249 + it], exp, tol,
            &format!("OH at 40.6 km it={it}"));
    }
}

#[test]
fn test_o3_column_decreases_with_altitude() {
    // Physical sanity check: O3 mixing ratio should peak in the mid-stratosphere
    let boxout = input_dir().join("boxout.dat");
    let rows = active_rows(&boxout);

    // Find the peak O3 box (should be around 34-38 km)
    let (peak_idx, _) = rows.iter().enumerate()
        .max_by(|(_, a), (_, b)| a[4].partial_cmp(&b[4]).unwrap())
        .unwrap();
    let peak_z = rows[peak_idx][0];
    assert!(peak_z > 30.0 && peak_z < 45.0,
        "O3 peak at {peak_z:.1} km, expected 30-45 km");
}

#[test]
fn test_n2o_increases_toward_troposphere() {
    // N2O should increase downward (more at lower altitudes from the tropospheric source)
    let boxout = input_dir().join("boxout.dat");
    let rows = active_rows(&boxout);

    let n2o_top = rows[0][5];   // highest active box
    let n2o_bot = rows[rows.len() - 1][5]; // lowest active box
    assert!(n2o_bot > n2o_top,
        "N2O should increase toward surface: bot={n2o_bot:.3E} top={n2o_top:.3E}");
}
