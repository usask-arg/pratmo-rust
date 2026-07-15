// bry_vs_tracers.gen.vers03.f → N2O-based tracer relationship module
// Subroutines: bry_vs_n2o_wamsley_org_2003, bry_vs_n2o_wamsley_org_tropics03,
//              ch4_vs_n2o_michelsen_ml_2003, cly_vs_n2o_solve, ch3cl_vs_n2o_solve

use crate::state::ModelState;

// ── Linear interpolation helper ───────────────────────────────────────────────

/// 1-D linear interpolation on table (ta, xa) at point tpt.
/// Fortran: SUBROUTINE LINEAR(xpt, tpt, xa, ta, npts, klast)
fn linear(tpt: f64, xa: &[f64], ta: &[f64]) -> f64 {
    let n = ta.len().min(xa.len());
    if n == 0 {
        return -999.0;
    }
    if tpt <= ta[0] {
        return if (tpt - ta[0]).abs() < 1e-15 {
            xa[0]
        } else {
            -999.0
        };
    }
    if tpt >= ta[n - 1] {
        return -999.0;
    }
    for i in 0..n - 1 {
        if tpt >= ta[i] && tpt < ta[i + 1] {
            let frac = (tpt - ta[i]) / (ta[i + 1] - ta[i]);
            return xa[i] + frac * (xa[i + 1] - xa[i]);
        }
    }
    -999.0
}

// ── CFC-11 from N2O (Wamsley eq. 16, mid-latitude) ───────────────────────────

fn cfc11_from_n2o_midlat(xn2otmp: f64) -> f64 {
    if (130.0..=310.0).contains(&xn2otmp) {
        252.92 - 4.5591 * xn2otmp + 0.025644 * xn2otmp * xn2otmp
            - 3.4639e-5 * xn2otmp * xn2otmp * xn2otmp
    } else if xn2otmp > 310.0 {
        272.09
    } else {
        17.52 * (xn2otmp / 130.0)
    }
}

// ── Wamsley Bry from CFC-11 (eq. 19) ─────────────────────────────────────────

fn bry_from_cfc11(cfc11_ppt: f64) -> f64 {
    let x1 = cfc11_ppt;
    let x2 = cfc11_ppt * cfc11_ppt;
    16.02 + 0.0033 * x1 + (-5.305e-4) * x2 + 2.55e-6 * x1 * x2 + (-5.37e-9) * x2 * x2
}

// ── Age-of-air correction tables ─────────────────────────────────────────────

const XN2O_AGE: [f64; 8] = [0., 20., 40., 60., 80., 100., 120., 140.];
const AGE: [f64; 8] = [6.65, 6.10, 5.70, 5.50, 5.30, 5.10, 4.90, 4.65];

fn age_correction(xn2otmp: f64, brypt: f64) -> f64 {
    if xn2otmp < 130.0 {
        let age_pt = linear(xn2otmp, &AGE, &XN2O_AGE);
        let age_130 = linear(130.0, &AGE, &XN2O_AGE);
        let age_diff = age_pt - age_130;
        let bry_scale = (1.0 + 2.96 / 100.0_f64).powf(age_diff);
        brypt / bry_scale
    } else {
        brypt
    }
}

// ── Bry vs N2O — mid-latitude, year 2003 ─────────────────────────────────────

/// Bry (ppt) from N2O (ppb), mid-latitude 2003 parameterisation.
/// Fortran: bry_vs_n2o_wamsley_org_2003
pub fn bry_vs_n2o_wamsley_org_2003(fn2o_ppb: f64) -> f64 {
    let xn2otmp = fn2o_ppb * 310.37 / 319.02;
    let cfc11_ppt = cfc11_from_n2o_midlat(xn2otmp);
    let del_year = 2003.5 - 1994.875;
    let bry_scale = (1.0 + 2.96 / 100.0_f64).powf(del_year);
    let brypt = (bry_from_cfc11(cfc11_ppt) * bry_scale).max(0.0);
    age_correction(xn2otmp, brypt).max(0.6)
}

// ── Bry vs N2O — tropics, year 2003 ──────────────────────────────────────────

/// Bry (ppt) from N2O (ppb), tropical 2003 parameterisation.
/// Fortran: bry_vs_n2o_wamsley_org_tropics03
pub fn bry_vs_n2o_wamsley_org_tropics03(fn2o_ppb: f64) -> f64 {
    let xn2otmp = fn2o_ppb * 310.37 / 319.02;
    let cfc11_ppt = if xn2otmp >= 225.0 {
        39.92 - 0.33539 * xn2otmp + 0.003_373_6 * xn2otmp * xn2otmp
    } else if xn2otmp >= 130.0 {
        cfc11_from_n2o_midlat(xn2otmp)
    } else {
        17.52 * (xn2otmp / 130.0)
    };
    let del_year = 2003.5 - 1994.875;
    let bry_scale = (1.0 + 2.96 / 100.0_f64).powf(del_year);
    let brypt = (bry_from_cfc11(cfc11_ppt) * bry_scale).max(0.0);
    age_correction(xn2otmp, brypt)
}

// ── CH4 vs N2O — Michelsen mid-latitude, year 2003 ───────────────────────────

/// CH4 (ppm) from N2O (ppb), Michelsen et al. 1998 + 2003 scaling.
/// Fortran: ch4_vs_n2o_michelsen_ml_2003
pub fn ch4_vs_n2o_michelsen_ml_2003(fn2o_ppb: f64) -> f64 {
    let x = fn2o_ppb;
    let x2 = x * x;
    let ch4 = if (0.0..45.0).contains(&x) {
        0.233_70 + 0.028_568 * x - 6.435_8e-4 * x2 + 6.018_6e-6 * x2 * x
    } else if x < 100.0 {
        0.504_99 + 6.969_0e-3 * x - 3.011_4e-5 * x2 + 7.232_1e-8 * x2 * x
    } else if x <= 205.0 {
        0.625_67 + 4.063_5e-3 * x - 5.673_8e-6 * x2
    } else if x < 280.0 {
        0.327_20 + 4.356_4e-3 * x
    } else if x <= 320.0 {
        0.348_41 + 4.263_8e-3 * x
    } else {
        0.0
    };
    let ch4_growth = (1753.67 - 1730.975) / 4.0 / 1000.0; // ppm/year
    let ch4_1993 = 1.666_f64;
    let ch4_2003 = ch4_1993 + ch4_growth * (2003.5 - 1993.5);
    ch4 * (ch4_2003 / ch4_1993)
}

// ── Cly vs N2O ────────────────────────────────────────────────────────────────

/// Total inorganic chlorine (ppb) from N2O (ppb). McLinden/Solve fit.
/// Fortran: cly_vs_n2o_solve
pub fn cly_vs_n2o_solve(fn2o_ppb: f64) -> f64 {
    if !(0.0..=320.0).contains(&fn2o_ppb) {
        return 0.0;
    }
    let x = fn2o_ppb;
    let cly = 3.538_76 - 2.677_09e-3 * x - 1.916_93e-5 * x * x - 2.405_84e-8 * x * x * x;
    cly.max(0.005)
}

// ── CH3Cl vs N2O ─────────────────────────────────────────────────────────────

/// CH3Cl (ppb) from N2O (ppb). Schauffler et al. 2003 fit.
/// Fortran: ch3cl_vs_n2o_solve
pub fn ch3cl_vs_n2o_solve(fn2o_ppb: f64) -> f64 {
    if !(0.0..=320.0).contains(&fn2o_ppb) {
        return 0.0;
    }
    let x = fn2o_ppb;
    578.98e-3 - 0.585_07e-3 * x - 7.219e-7 * x * x - 9.200_2e-9 * x * x * x
}

// ── High-level tracer setter ──────────────────────────────────────────────────

/// Set mixing ratios for long-lived tracers in box `ib` from the N2O profile.
pub fn set_tracers_from_n2o(s: &mut ModelState, ib: usize) {
    let fn2o_ppb = s.fn2o[ib] * 1.0e9;
    let ialt = (s.nboxdo[ib].unsigned_abs() as usize).saturating_sub(1);
    let densty = s.dm[ialt];

    s.fch4[ib] = ch4_vs_n2o_michelsen_ml_2003(fn2o_ppb) * 1.0e-6;
    s.fh2o[ib] = (7.0e-6 - 2.0 * s.fch4[ib]).max(0.0);
    s.fco[ib] = 3.0e-8;

    let brypt = bry_vs_n2o_wamsley_org_2003(fn2o_ppb);
    s.fbrx[ib] = (brypt * 1.0e-12 + 2.0e-12).max(5.0e-13);

    let clypt = cly_vs_n2o_solve(fn2o_ppb);
    s.fclx[ib] = clypt * 1.0e-9;

    let ch3clpt = ch3cl_vs_n2o_solve(fn2o_ppb);
    s.fch3cl[ib] = ch3clpt * 1.0e-9;

    s.fnoy[ib] = s.dn2oref[ialt] / densty;
    s.fo3[ib] = s.do3ref[ialt] / densty;
}

#[cfg(test)]
mod tests {
    use super::{ch3cl_vs_n2o_solve, cly_vs_n2o_solve};

    #[test]
    fn tracer_fits_return_zero_outside_their_calibrated_domain() {
        // The original Fortran typo assigned an implicit clypt variable in
        // the CH3Cl branch, leaving the output undefined for out-of-range N2O.
        // Rust deliberately returns the documented zero sentinel for both fits.
        for n2o in [-1.0, 321.0] {
            assert_eq!(cly_vs_n2o_solve(n2o), 0.0);
            assert_eq!(ch3cl_vs_n2o_solve(n2o), 0.0);
        }
    }

    #[test]
    fn tracer_fits_cover_both_calibration_endpoints() {
        assert!(cly_vs_n2o_solve(0.0) >= 0.005);
        assert!(cly_vs_n2o_solve(320.0) >= 0.005);
        assert!(ch3cl_vs_n2o_solve(0.0).is_finite());
        assert!(ch3cl_vs_n2o_solve(320.0).is_finite());
    }
}
