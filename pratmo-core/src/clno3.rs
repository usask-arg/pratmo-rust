// ClNO3code.f → alternative heterogeneous chemistry (older Shi et al. formulation)
// This is the earlier version of hetprob that predates bhetp.f.
// The main model uses bhetp.f (heterogeneous.rs); this module is kept for reference.

use crate::state::ModelState;

/// Compute ClONO2/HOCl heterogeneous reaction probabilities.
/// Simplified formulation (ClNO3code.f) using only T, r, pH2O, pHCl, pClONO2.
/// Returns [gClONO2_HCl, gClONO2_H2O, gHOCl, gN2O5_HCl, gN2O5_H2O,
///           gBrONO2, gHOBr_HCl, gHOBr_HBr].
/// Fortran: subroutine hetprob in ClNO3code.f
pub fn hetprob_clno3(t: f64, _r: f64, ph2o: f64) -> [f64; 8] {
    // Initial estimate from wt%
    let wt =
        (t * (0.6246 * ph2o.ln() - 14.458) + 3565.0) / (44.777 + 1.3204 * ph2o.ln() - 0.199_88 * t);
    let g_h2o = 10.0_f64.powf(1.86 - 0.0747 * wt);
    let g_hcl = 0.1 * g_h2o;
    let g_hocl: f64 = 1.0e-13;
    let g_n2o5_hcl: f64 = 1.0e-4;
    let g_n2o5_h2o: f64 = 0.1;
    let g_hobr_hcl: f64 = 0.2;
    let g_hobr_hbr: f64 = 0.25;
    let g_brono2: f64 = 0.8;

    if t >= 250.0 {
        return [1e-10, 1e-10, 1e-10, 1e-10, 0.1, 0.8, 0.0, 0.0];
    }

    [
        g_hcl, g_h2o, g_hocl, g_n2o5_hcl, g_n2o5_h2o, g_brono2, g_hobr_hcl, g_hobr_hbr,
    ]
}

/// Stub for any ClNO3-specific chemistry applied from bchem.f.
/// Not called in the standard model path (superseded by bhetp.f).
pub fn clno3_chemistry(_s: &mut ModelState) {}
