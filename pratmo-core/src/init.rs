// bctin.f → box density initialisation module
// CTINIT, GEM_SP

use crate::{
    constants::NDEN,
    solver::rplace,
    state::ModelState,
    tracers::{
        bry_vs_n2o_wamsley_org_2003, bry_vs_n2o_wamsley_org_tropics03, ch3cl_vs_n2o_solve,
        ch4_vs_n2o_michelsen_ml_2003, cly_vs_n2o_solve,
    },
};

/// Initialise species mixing ratios for box `ib` based on N2O, O3, and family
/// mixing ratio empirical relationships.
///
/// Arguments:
///   ib      — box index (0-based)
///   densbx  — number density at the box altitude (cm⁻³)
///   lat     — latitude (degrees, integer sign)
///   mon     — month (1–12)
///
/// Fortran: SUBROUTINE CTINIT(IB, DENSBX, LAT, MON)
pub fn ctinit(s: &mut ModelState, ib: usize, densbx: f64, lat: i32, _mon: i32) {
    // Fixed mixing ratios not set by empirical relationships
    s.fco[ib] = 20.0e-9;
    s.fcf2cl[ib] = 100.0e-12;
    s.fcfcl3[ib] = 100.0e-12;
    s.fccl4[ib] = 100.0e-12;
    s.fch3cl[ib] = 100.0e-12;
    s.fmecl[ib] = 100.0e-12;
    s.fch3br[ib] = 10.0e-12;
    s.fh2[ib] = 0.5e-6;

    // FO3, FNOY, FN2O are already set by the caller (CTMLFQ) before calling CTINIT.
    // fn2o_ppb is used to drive the empirical tracer relationships.
    let ialt = (s.nboxdo[ib].unsigned_abs() as usize).saturating_sub(1);
    let fn2o_ppb = s.fn2o[ib] * 1.0e9;

    // CH4 from N2O (Michelsen 2003)
    s.fch4[ib] = ch4_vs_n2o_michelsen_ml_2003(fn2o_ppb) * 1.0e-6;

    // H2O from 2-CH4 parameterisation
    s.fh2o[ib] = (7.0e-6 - 2.0 * s.fch4[ib]).max(0.0);
    s.fco[ib] = 3.0e-8;

    // Bry — mid-latitude relation (Wamsley 2003)
    let brypt = bry_vs_n2o_wamsley_org_2003(fn2o_ppb);
    s.fbrx[ib] = brypt * 1.0e-12;
    let dfbry = s.fbrx[ib];

    // Override with tropical relation for |lat| ≤ 30
    if lat.abs() <= 30 {
        let brypt_trop = bry_vs_n2o_wamsley_org_tropics03(fn2o_ppb);
        s.fbrx[ib] = brypt_trop * 1.0e-12;
        let _dfbry = s.fbrx[ib] - dfbry; // diagnostic only
    }
    s.fbrx[ib] += 2.0e-12;
    s.fbrx[ib] = s.fbrx[ib].max(5.0e-13);

    // Cly (Solve fit)
    let clypt = cly_vs_n2o_solve(fn2o_ppb);
    s.fclx[ib] = clypt * 1.0e-9;

    // Override Bry from reference profile if IWBRY==1
    if s.iwbry == 1 && s.fbrx_ref[ialt] > 1.0e-30 {
        s.fbrx[ib] = s.fbrx_ref[ialt];
    }

    // CH3Cl (Schauffler 2003)
    let ch3clpt = ch3cl_vs_n2o_solve(fn2o_ppb);
    s.fch3cl[ib] = ch3clpt * 1.0e-9;

    // ── Uniform initialisation of number densities ────────────────────────

    let xxxnoy = s.fnoy[ib] * densbx;
    s.dhno3[ib] = 0.970 * xxxnoy;
    s.dno[ib] = 0.010 * xxxnoy;
    s.dno2[ib] = 0.010 * xxxnoy;
    s.dno3[ib] = 0.001 * xxxnoy;
    s.dn2o5[ib] = 0.001 * xxxnoy;
    s.dhno2[ib] = 0.001 * xxxnoy;
    s.dhno4[ib] = 0.005 * xxxnoy;
    s.dclno3[ib] = 0.000_99 * xxxnoy;
    s.dbrno3[ib] = 0.000_01 * xxxnoy;

    let xxxcly = s.fclx[ib] * densbx - s.dclno3[ib] - s.dbrno3[ib];
    if xxxcly < 0.0 {
        eprintln!(
            "  CTINIT: negative Cly at box/lat/mon={}/{}/{}  xxxcly={}",
            ib, lat, _mon, xxxcly
        );
        // In the Fortran this is STOP — we continue with zero
    }
    let xxxcly = xxxcly.max(0.0);

    s.dhcl[ib] = 0.950 * xxxcly;
    s.dcl[ib] = 0.002 * xxxcly;
    s.dcl2[ib] = 0.001 * xxxcly;
    s.dclo[ib] = 0.025 * xxxcly;
    s.dhocl[ib] = 0.020 * xxxcly;
    s.dcl2o2[ib] = 0.001 * xxxcly;
    s.doclo[ib] = 0.000_9 * xxxcly;
    s.dbrcl[ib] = 0.000_1 * xxxcly;

    let xxxbry = s.fbrx[ib] * densbx - s.dbrcl[ib] - s.dbrno3[ib];
    if xxxbry < 0.0 {
        eprintln!(
            "  CTINIT: negative Bry at box/lat/mon={}/{}/{}  xxxbry={}",
            ib, lat, _mon, xxxbry
        );
    }
    let xxxbry = xxxbry.max(0.0);

    s.dbro[ib] = 0.800 * xxxbry;
    s.dbr[ib] = 0.050 * xxxbry;
    s.dhbr[ib] = 0.100 * xxxbry;
    s.dhobr[ib] = 0.050 * xxxbry;

    // Iy — total inorganic iodine from fiodx (default 1 ppt)
    // Partition following midday lower-stratosphere estimates
    if s.liod {
        let xxxiody = (s.fiodx[ib] * densbx).max(0.0);
        let seed = 1.0e-30 * densbx;
        s.di_[ib] = 0.200 * xxxiody; // I
        s.dio[ib] = 0.500 * xxxiody; // IO
        s.dhoi[ib] = 0.150 * xxxiody; // HOI
        s.diono2[ib] = 0.100 * xxxiody; // IONO2
        s.dhi[ib] = 0.050 * xxxiody; // HI
        s.doio[ib] = seed;
        s.di2[ib] = seed;
        s.di2o2[ib] = seed;
        s.di2o3[ib] = seed;
        s.di2o4[ib] = seed;
    }

    let xxxppt = 1.0e-12 * densbx;
    s.dh[ib] = 0.000_001 * xxxppt;
    s.doh[ib] = 0.010 * xxxppt;
    s.dho2[ib] = 0.010 * xxxppt;
    s.dh2o2[ib] = 0.010 * xxxppt;
    s.dh2co[ib] = 0.010 * xxxppt;
    s.droo[ib] = 0.010 * xxxppt;
    s.drooh[ib] = 0.010 * xxxppt;
    s.do3[ib] = s.fo3[ib] * densbx;
    s.do_[ib] = 0.0001 * s.do3[ib];

    // Initialise all diurnal time steps with first-guess values
    let mut xn = [0.0f64; NDEN];
    rplace(s, &mut xn, ib);
    let ntotx = s.ntotx;
    let ntimdo = s.ntimdo;
    for j in 0..ntotx {
        for i in 0..ntimdo {
            s.xxnoft[[j, i, ib]] = xn[j];
        }
    }
    // Also prime XNOFT column 0
    for j in 0..ntotx {
        s.xnoft[[j, 0]] = xn[j];
    }
}

// ── GEM_SP ───────────────────────────────────────────────────────────────────

/// Saturation vapour pressure and relative humidity.
/// Fortran: SUBROUTINE GEM_SP(TTT, PRS, FOEW8, rh)
pub fn gem_sp(ttt: f64, prs: f64) -> (f64, f64) {
    const TRPL: f64 = 273.16;
    const EPS1: f64 = 0.621_948_002_210_14;
    const EPS2: f64 = 1.0 - EPS1;

    let sign_val = if ttt >= TRPL { 17.2690 } else { 21.8750 };
    let foew8 = 610.780
        * (sign_val * (ttt - TRPL).abs()
            / (ttt - 35.860 + 0.0_f64.max(28.20 * (TRPL - ttt).signum().max(0.0))))
        .min(sign_val * (ttt - TRPL).abs() / (ttt - 35.860))
        .exp();

    let rh = 1.0 / (1.0_f64.max(prs / foew8) - EPS2);
    (foew8, rh)
}
