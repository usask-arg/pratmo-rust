// bhetp.f → heterogeneous chemistry module
// hetprob: Shi et al. / Worsnop stratospheric reaction probabilities

use crate::state::ModelState;

/// Compute reaction probabilities γ for 8 heterogeneous reactions.
///
/// Arguments:
///   pres    — total pressure (hPa or same units as ph2o)
///   temp    — temperature (K)
///   r       — aerosol radius (cm), nominally 1e-5
///   ph2o    — H2O partial pressure
///   phcl    — HCl partial pressure (atm)
///   pclono2 — ClONO2 partial pressure (atm)
///   phbr    — HBr partial pressure (atm)
///   lice    — true for PSC ice conditions (use fixed coefficients)
///
/// Returns [γ1..γ8]:
///   [0] gClONO2_HCl, [1] gClONO2_H2O, [2] gHOCl_HCl,
///   [3] gN2O5_HCl,   [4] gN2O5_H2O,   [5] gBrONO2,
///   [6] gHOBr_HCl,   [7] gHOBr_HBr
///
/// Fortran: SUBROUTINE HETPROB(PRES,T,r,PH2O,PHCl,PClONO2,PHBR,g1..g8,lice)
pub fn hetprob(
    pres: f64,
    temp: f64,
    r: f64,
    ph2o: f64,
    phcl: f64,
    pclono2: f64,
    phbr: f64,
    lice: bool,
) -> [f64; 8] {
    // PSC ice conditions — fixed laboratory values
    if lice {
        return [
            0.3,   // gClONO2_HCl
            0.3,   // gClONO2_H2O
            0.2,   // gHOCl
            0.05,  // gN2O5_HCl
            0.024, // gN2O5_H2O
            0.3,   // gBrONO2
            0.3,   // gHOBr_HCl
            0.0,   // gHOBr_HBr
        ];
    }

    // Liquid sulfate aerosol not applicable at high T or high H2O/pressure ratio
    if ph2o / pres > 20.0e-6 || temp >= 250.0 {
        return [
            0.0, 0.0, 0.0, 0.0,
            0.1,  // gN2O5_H2O
            0.8,  // gBrONO2
            0.0, 0.0,
        ];
    }

    let t = temp;

    // ── Table 1: H2SO4 wt% from T and pH2O ──────────────────────────────────

    let p0h2o = (18.452_406_985
        - 3_505.157_880_7 / t
        - 330_918.550_82 / (t * t)
        + 12_725_068.262 / (t * t * t))
        .exp();

    let aw = ph2o / p0h2o;

    let (a1, b1, c1, d1, a2, b2, c2, d2) = if aw < 0.05 {
        (
            12.372_089_32, -0.161_255_161_14, -30.490_657_554, -2.113_311_424_1,
            13.455_394_705, -0.192_131_225_5, -34.285_174_607, -1.762_007_307_8,
        )
    } else if aw < 0.85 {
        (
            11.820_654_354, -0.207_864_042_44, -4.807_306_373, -5.172_754_034_8,
            12.891_938_068, -0.232_338_477_08, -6.426_123_775_7, -4.900_547_131_9,
        )
    } else {
        (
            -180.065_410_28, -0.386_011_025_92, -93.317_846_778, 273.881_322_45,
            -176.958_140_97, -0.362_570_481_54, -90.469_744_201, 267.455_099_88,
        )
    };

    let y1 = a1 * aw.powf(b1) + c1 * aw + d1;
    let y2 = a2 * aw.powf(b2) + c2 * aw + d2;
    let m  = y1 + (t - 190.0) * (y2 - y1) / 70.0;
    let wt = 9800.0 * m / (98.0 * m + 1000.0);

    // ── Table 2: Parameters for H2SO4 solution ───────────────────────────────

    let z1  = 0.123_64 - 5.6e-7 * t * t;
    let z2  = -0.029_54 + 1.814e-7 * t * t;
    let z3  = 2.343e-3 - 1.487e-6 * t - 1.324e-8 * t * t;
    let rho = 1.0 + z1 * m + z2 * m.powf(1.5) + z3 * m * m;
    let mso4 = rho * wt / 9.8;
    let x   = wt / (wt + (100.0 - wt) * 98.0 / 18.0);
    let aa  = 169.5 + 5.18 * wt - 0.0825 * wt * wt + 3.27e-3 * wt * wt * wt;
    let t_o = 144.11 + 0.166 * wt - 0.015 * wt * wt + 2.18e-4 * wt * wt * wt;
    let h   = aa * t.powf(-1.43) * (448.0 / (t - t_o)).exp();
    let ah  = (60.51 - 0.095 * wt + 0.0077 * wt * wt - 1.61e-5 * wt * wt * wt
        - (1.76 + 2.52e-4 * wt * wt) * t.sqrt()
        + (-805.89 + 253.05 * wt.powf(0.076)) / t.sqrt())
        .exp();

    // ── Table 3: ClONO2 + H2O and ClONO2 + HCl ──────────────────────────────

    let c_clono2  = 1474.0 * t.sqrt();
    let sclono2   = 0.306 + 24.0 / t;
    let hclono2   = 1.6e-6 * (4710.0 / t).exp() * (-sclono2 * mso4).exp();
    let dclono2   = 5.0e-8 * t / h;
    let kh2o      = 1.95e10 * (-2800.0 / t).exp();
    let kh        = 1.22e12 * (-6200.0 / t).exp();
    let khydr     = kh2o * aw + kh * ah * aw;
    let gbh2o     = 4.0 * hclono2 * 0.082 * t * (dclono2 * khydr).sqrt() / c_clono2;
    let hhcl      = (0.094 - 0.61 * x + 1.2 * x * x)
        * (-8.68 + (8515.0 - 10718.0 * x.powf(0.7)) / t).exp();
    let mhcl      = hhcl * phcl;
    let khcl      = 7.9e11 * ah * dclono2 * mhcl;
    let lclono2   = (dclono2 / (khydr + khcl)).sqrt();
    let fclono2   = 1.0 / (r / lclono2).tanh() - lclono2 / r;
    let gclono2rxn = fclono2 * gbh2o * (1.0 + khcl / khydr).sqrt();
    let gbhcl     = gclono2rxn * khcl / (khcl + khydr);
    let gs        = 66.12 * (-1374.0 / t).exp() * hclono2 * mhcl;
    let fhcl      = 1.0 / (1.0 + 0.612 * (gs + gbhcl) * pclono2 / phcl.max(1e-300));
    let gsp       = fhcl * gs;
    let gbhclp    = fhcl * gbhcl;
    let gb        = gbhclp + gclono2rxn * khydr / (khcl + khydr);
    let gclono2   = 1.0 / (1.0 + 1.0 / (gsp + gb));

    let g_clono2_hcl = gclono2 * (gsp + gbhclp) / (gsp + gb);
    let g_clono2_h2o = gclono2 - g_clono2_hcl;

    // ── Table 4: HOCl + HCl ─────────────────────────────────────────────────

    let c_hocl   = 2009.0 * t.sqrt();
    let shocl    = 0.0776 + 59.18 / t;
    let hhocl    = 1.91e-6 * (5862.4 / t).exp() * (-shocl * mso4).exp();
    let dhocl    = 6.4e-8 * t / h;
    let khocl_hcl = 1.25e9 * ah * dhocl * mhcl;
    let ghoclrxn = 4.0 * hhocl * 0.082 * t * (dhocl * khocl_hcl).sqrt() / c_hocl;
    let lhocl    = (dhocl / khocl_hcl).sqrt();
    // Numerically stable approximation (Nick Lloyd 2015-05-26):
    let xx       = r / lhocl;
    let fhocl    = xx * 0.5;
    let g_hocl   = 1.0 / (1.0 + 1.0 / (fhocl * ghoclrxn * fhcl));

    // ── HOBr + HCl and HOBr + HBr ───────────────────────────────────────────

    let c_hobr    = 1477.0 * t.sqrt();
    let khobr_hcl = 10.0_f64.powf(6.08 - 1050.0 / t + 0.0747 * wt) * hhcl * phcl / 5.0;
    let hhbr      = hhcl * 200.0;
    let khobr_hbr = 10.0_f64.powf(6.08 - 1050.0 / t + 0.0747 * wt) * hhbr * phbr / 15.0;
    let lhobr     = (dhocl / (khobr_hcl + khobr_hbr)).sqrt();
    let fhobr     = 1.0 / (r / lhobr).tanh() - lhobr / r;
    let ghobrrxn  = ghoclrxn * (c_hocl / c_hobr)
        * ((khobr_hcl + khobr_hbr) / khocl_hcl).sqrt()
        * 10.0;
    let g_hobr    = 1.0 / (1.0 + 1.0 / (fhobr * ghobrrxn));
    let pp        = khobr_hbr / khobr_hcl.max(1e-300);
    let g_hobr_hcl = g_hobr / (1.0 + pp);
    let g_hobr_hbr = g_hobr * pp / (1.0 + pp);

    // Temperature-independent values (JPL recommendations)
    let g_n2o5_h2o = 0.1;
    let g_brono2   = 0.75;
    let g_n2o5_hcl = 1.0e-4;

    [
        g_clono2_hcl, // [0]
        g_clono2_h2o, // [1]
        g_hocl,       // [2]
        g_n2o5_hcl,   // [3]
        g_n2o5_h2o,   // [4]
        g_brono2,     // [5]
        g_hobr_hcl,   // [6]
        g_hobr_hbr,   // [7]
    ]
}

/// Apply heterogeneous reaction rates for the current box.
/// Stub — actual application is done inside chemistry::apply_het_rates_setupr.
pub fn apply_het_rates(s: &mut ModelState) {
    let _ = s;
}
