// bjval.f → photolysis (J-value) module
// SOL, JVALUE, OPTAU, SCATTR, MATIN3, SPHERE + helper functions

use crate::{
    constants::{NL, NTAU},
    state::ModelState,
};

// ── SOL — entry point for J-value update ─────────────────────────────────────

/// Compute J-values at all altitudes (LPRTJV mode) or all boxes (normal mode)
/// for the current cos(SZA) stored in s.gmu.
/// Fortran: SUBROUTINE SOL
pub fn sol(s: &mut ModelState, gmu: f64) {
    s.gmu = gmu;
    s.rflect = s.clouds.clamp(0.0, 1.0);
    s.u0 = gmu;
    s.sza = gmu.acos() * 57.29578;

    // Find matching diurnal time step index for albedo lookup
    let ntimdo = s.ntimdo;
    for i in 0..ntimdo {
        if (s.sza - s.ztime[i]).abs() < 1e-10 {
            s.izalb = i;
        }
    }

    // Zero all J-value profiles. In full-profile mode the leading dimension is
    // altitude; in box mode it is box number.
    let nc = s.nc;
    let njval = s.njval;
    let nlocations = if s.lprtjv { nc } else { s.nbox };
    for iv in 0..njval {
        for j in 0..nlocations {
            s.jval_set(j, iv, 0.0);
        }
    }

    if s.u0 < s.gmu0 {
        return; // night — all J = 0
    }

    // Compute actinic flux FFF(K,J) and base J-values JO2, JNO via JVALUE
    jvalue(s);

    // Compute remaining J-values from temperature-dependent cross-sections
    let nw1 = s.nw1.saturating_sub(1); // 0-based
    let nw2 = s.nw2.saturating_sub(1);

    if s.lprtjv {
        // Full altitude profile mode
        for j in 0..nc {
            let tt = s.t[j];
            // O3 total cross-section and O(1D) quantum yield
            for k in nw1..=nw2 {
                let qo3tot = xseco3(k, tt, s);
                let qo31d = xsec1d(k, tt, s);
                let fff_kj = s.fff[[k, j]];
                let vo3_j = s.vo3[j] + fff_kj * qo3tot;
                let vo3d_j = s.vo3d[j] + fff_kj * qo3tot * qo31d;
                s.vo3[j] = vo3_j;
                s.vo3d[j] = vo3d_j;
            }
            // Additional J-values (IV = 5..NJVAL, 1-based; index 4..njval-1, 0-based)
            for iv in 4..njval {
                let tfact = interp_tfact(tt, s.tqq[[0, iv]], s.tqq[[1, iv]]);
                let mut sum = 0.0;
                for k in nw1..=nw2 {
                    let jq = iv - 4; // index into QQQ third dimension (0-based)
                    let qqqt = s.qqq[[k, 0, jq]] + (s.qqq[[k, 1, jq]] - s.qqq[[k, 0, jq]]) * tfact;
                    sum += qqqt * s.fff[[k, j]];
                }
                let cur = s.jval_get(j, iv);
                s.jval_set(j, iv, cur + sum);
            }
        }
    } else {
        // Per-box mode: J-values indexed by box number (JB → 0-based jb)
        let nbox = s.nbox;
        for jb in 0..nbox {
            let j = (s.nboxdo[jb].unsigned_abs() as usize).saturating_sub(1); // 0-based altitude level
            let tt = s.t[j] + s.boxtt[jb];
            for k in nw1..=nw2 {
                let qo3tot = xseco3(k, tt, s);
                let qo31d = xsec1d(k, tt, s);
                let fff_kj = box_actinic_flux(s, jb, k);
                s.vo3[jb] += fff_kj * qo3tot;
                s.vo3d[jb] += fff_kj * qo3tot * qo31d;
            }
            for iv in 4..njval {
                let tfact = interp_tfact(tt, s.tqq[[0, iv]], s.tqq[[1, iv]]);
                let mut sum = 0.0;
                for k in nw1..=nw2 {
                    let jq = iv - 4;
                    let qqqt = s.qqq[[k, 0, jq]] + (s.qqq[[k, 1, jq]] - s.qqq[[k, 0, jq]]) * tfact;
                    sum += qqqt * box_actinic_flux(s, jb, k);
                }
                let cur = s.jval_get(jb, iv);
                s.jval_set(jb, iv, cur + sum);
            }
        }
    }
}

/// Temperature interpolation factor clamped to [0, 1].
#[inline(always)]
fn interp_tfact(tt: f64, tq1: f64, tq2: f64) -> f64 {
    if tq2 <= tq1 {
        0.0
    } else {
        ((tt - tq1) / (tq2 - tq1)).clamp(0.0, 1.0)
    }
}

/// Resolve the two radiative shells used for one chemistry box. Legacy runs
/// do not populate the explicit mapping, so their chemistry level remains the
/// fallback.
#[inline(always)]
fn box_flux_indices(s: &ModelState, jb: usize) -> (usize, usize, f64) {
    let chemistry_level = (s.nboxdo[jb].unsigned_abs() as usize)
        .saturating_sub(1)
        .min(s.nc.saturating_sub(1));
    let lower = s.box_flux_lower[jb];
    let upper = s.box_flux_upper[jb];
    let weight = s.box_flux_upper_weight[jb];
    let mapping_is_configured = lower < s.nc
        && upper < s.nc
        && lower <= upper
        && weight.is_finite()
        && (0.0..=1.0).contains(&weight)
        && (lower != 0 || upper != 0 || chemistry_level == 0);
    if mapping_is_configured {
        (lower, upper, weight)
    } else {
        (chemistry_level, chemistry_level, 0.0)
    }
}

#[inline(always)]
fn box_actinic_flux(s: &ModelState, jb: usize, wavelength: usize) -> f64 {
    let (lower, upper, upper_weight) = box_flux_indices(s, jb);
    let lower_flux = s.fff[[wavelength, lower]];
    lower_flux + upper_weight * (s.fff[[wavelength, upper]] - lower_flux)
}

// ── JVALUE — wavelength loop + actinic flux ───────────────────────────────────

/// Compute actinic flux FFF and J(O2)/J(NO) by iterating over wavelengths.
/// Fortran: SUBROUTINE JVALUE
fn jvalue(s: &mut ModelState) {
    let nc = s.nc;
    let nw2 = s.nw2.saturating_sub(1); // 0-based

    // Zero actinic flux and JO2/JNO
    for j in 0..nc {
        for k in 0..=nw2 {
            s.fff[[k, j]] = 0.0;
        }
    }
    let nlocations = if s.lprtjv { nc } else { s.nbox };
    for j in 0..nlocations {
        s.vno[j] = 0.0;
        s.vo2[j] = 0.0;
    }

    if s.sza > 98.0 {
        return;
    }

    // Spherical geometry air mass weights
    sphere(s);

    let nw1 = s.nw1.saturating_sub(1); // 0-based
    let nwsrb = s.nwsrb; // 0-based count
    let nsr = s.nsr;
    // ── Main wavelength loop ──────────────────────────────────────────────────
    for k in nw1..=nw2 {
        let wave = s.wl[k];
        // SRB band index (0-based): KSR=k-NWSRB in Fortran (1-based)
        // Fortran: KSR = K-NWSRB; if KSR<1 or KSR>NSR then non-SRB
        let ksr_1based = (k + 1) as i64 - nwsrb as i64;
        let (ksr, n_odf, is_srb) = if ksr_1based >= 1 && ksr_1based <= nsr as i64 {
            let ksr0 = (ksr_1based - 1) as usize; // 0-based
            (ksr0, s.isr[ksr0], true)
        } else {
            (0usize, 1usize, false)
        };

        // O3 cross-sections at each altitude level
        let mut xqo3 = [0.0f64; NL];
        for j in 0..nc {
            xqo3[j] = xseco3(k, s.t[j], s);
        }

        // Rayleigh and Mie aerosol cross-sections
        let xqray = raylay(wave);
        let xqaer = sctmie(wave, s.turbdx);

        // ── Loop over opacity distribution functions (=1 outside SRB) ────────
        for kodf in 0..n_odf {
            // O2 cross-sections
            let mut xqo2 = [0.0f64; NL];
            let ksr_arg = if is_srb { ksr } else { 0 }; // 0-based, but xseco2 needs 1-based originally
            for j in 0..nc {
                xqo2[j] = xseco2(k, s.t[j], ksr_arg, kodf, is_srb, s);
            }

            // Compute optical depths and attenuated beam
            let trans = optau(s, &xqo2, &xqo3, xqray, xqaer);

            let j1 = s.nlbatm.saturating_sub(1); // 0-based

            // Accumulate actinic flux
            let fl_k = s.fl[k];
            let flscal = s.flscal;
            let ssf_k = s.ssf[k];
            let odf_ksr = if is_srb { s.odf[[kodf, ksr]] } else { 1.0 };

            for j in j1..nc {
                let fopt = fl_k * trans[j] * flscal * ssf_k * if is_srb { odf_ksr } else { 1.0 };
                s.fff[[k, j]] += fopt;

                if s.lprtjv {
                    s.vo2[j] += xqo2[j] * fopt;
                    if is_srb && s.fno[ksr] > 0.0 {
                        s.vno[j] += s.qno[[kodf, ksr]] * fopt;
                    }
                }
            }

            // Per-box JO2/JNO (normal mode)
            if !s.lprtjv {
                let nbox = s.nbox;
                for jb in 0..nbox {
                    let j = (s.nboxdo[jb].unsigned_abs() as usize).saturating_sub(1); // 0-based
                    let (lower, upper, upper_weight) = box_flux_indices(s, jb);
                    let box_trans = trans[lower] + upper_weight * (trans[upper] - trans[lower]);
                    let fl_k = s.fl[k];
                    let fopt = fl_k
                        * box_trans
                        * s.flscal
                        * s.ssf[k]
                        * if is_srb { s.odf[[kodf, ksr]] } else { 1.0 };
                    let tt = s.t[j] + s.boxtt[jb];
                    let xqo2b = xseco2(k, tt, ksr_arg, kodf, is_srb, s);
                    s.vo2[jb] += xqo2b * fopt;
                    if is_srb && s.fno[ksr] > 0.0 {
                        s.vno[jb] += s.qno[[kodf, ksr]] * fopt;
                    }
                }
            }
        }
    }
}

// ── OPTAU — optical depth and attenuated beam ─────────────────────────────────

/// Compute slant-path transmission TRANS(J) for each altitude level.
/// Fortran: SUBROUTINE OPTAU(XQO2, XQO3, XQRAY, XQAER, TRANS)
fn optau(
    s: &mut ModelState,
    xqo2: &[f64; NL],
    xqo3: &[f64; NL],
    xqray: f64,
    xqaer: f64,
) -> [f64; NL] {
    let nc = s.nc;
    let j1 = s.nlbatm.saturating_sub(1); // 0-based

    let mut trans = [0.0f64; NL];
    let mut dtau = [0.0f64; NL];
    let mut piray = [0.0f64; NTAU]; // local work arrays, mirror s.piray/piaer
    let mut piaer = [0.0f64; NTAU];

    // Build differential optical depths
    for j in j1..nc {
        let xlo3 = s.do3ref[j] * xqo3[j];
        let xlo2 = s.dm[j] * s.po2 * xqo2[j];
        let xlray = s.dm[j] * xqray;
        let xlaer = if s.radiative_aerosol {
            s.aer[j].max(0.0) * xqaer
        } else {
            0.0
        };
        dtau[j] = xlo3 + xlo2 + xlray + xlaer;
        let d = dtau[j].max(1e-300);
        piray[j] = xlray / d;
        piaer[j] = xlaer / d;
    }

    // Copy boundary values (Fortran: PIRAY(1)=PIRAY(2), using NLBATM indexing)
    piray[j1] = piray[j1 + 1];
    piaer[j1] = piaer[j1 + 1];

    // Build cumulative TAU array (in local storage mirroring s.tau)
    // TAU(1)=0 = top; TAU(2)=DTAU(NC)*ZZHT; TAU(JI)=TAU(JI-1)+trapezoidal
    let ntt = nc + 2 - j1 - 1; // Fortran: NTT = NC+2-J1 (1-based J1)
    let mut tau_local = vec![0.0f64; ntt + 2];
    tau_local[0] = 0.0;
    tau_local[1] = dtau[nc - 1] * s.zzht;
    for ji in 2..ntt {
        let j = nc + 1 - ji - 1; // 0-based
        let jp1 = j + 1;
        tau_local[ji] = tau_local[ji - 1] + (s.z[jp1] - s.z[j]) * (dtau[j] + dtau[jp1]) * 0.5;
    }

    // Store for SCATTR and for per-layer access
    s.ntt = ntt;
    for i in 0..ntt.min(NTAU) {
        s.tau[i] = tau_local[i];
        // The original Fortran passes PIRAY/PIAER in atmospheric-level order
        // to SCATTR (rather than reversing them with TAU's optical-depth
        // order). Preserve that indexing quirk only in parity mode; the
        // normal build uses the physically aligned order.
        if cfg!(feature = "fortran-parity") {
            s.piray[i] = piray[i].max(0.0);
            s.piaer[i] = piaer[i].max(0.0);
        } else {
            s.piray[i] = piray[nc + 1 - i - 1].max(0.0);
            s.piaer[i] = piaer[nc + 1 - i - 1].max(0.0);
        }
    }
    // Keep piray/piaer boundaries
    s.piray[0] = s.piray[1];
    s.piaer[0] = s.piaer[1];

    // Average mode (NDAY=2or3, LDIURN=false): iterate over cos(SZA) points
    let wting_base = 1.0_f64;
    let nloop = if s.ldiurn { 1 } else { s.nmu };

    for iloop in 0..nloop {
        let (wting, u0_loop) = if s.ldiurn {
            (wting_base, s.u0)
        } else {
            // Fortran: WTIME(ILOOP) + WTIME(NTIM+1-ILOOP).
            // Both indices are zero-based here, so the mirrored slot is
            // `NTIM - 1 - ILOOP`.
            let w = s.wtime[iloop] + s.wtime[s.ntim - 1 - iloop];
            let u = s.utime[iloop];
            if u.acos() * 57.29578 > 98.0 {
                continue;
            }
            // Redo SPHERE for new angle
            s.u0 = u;
            s.sza = u.acos() * 57.29578;
            sphere(s);
            (w, u)
        };

        // Compute attenuated beam exp(-TAU_slant/U0) at each level
        s.ttfbot = 0.0;
        for j in j1..nc {
            let ji = nc + 1 - j - 1; // 0-based into tau/fltau
            if s.wtau[[j, j]] > 0.0 {
                let mut xltau = 0.0;
                for i in 0..nc {
                    xltau += dtau[i] * s.wtau[[i, j]];
                }
                let xftau = (-xltau).exp();
                s.fltau[ji] = xftau;
                if j == j1 && u0_loop > 0.0 {
                    s.ttfbot = xftau * u0_loop;
                }
            } else {
                s.fltau[ji] = 0.0;
            }
        }
        s.fltau[0] = 1.0;

        // Multiple scattering
        scattr(s);

        // Accumulate into TRANS
        for j in j1..nc {
            let ji = nc + 1 - j - 1; // 0-based
            trans[j] += wting * s.fltau[ji];
        }

        // Restore U0/SZA if we changed them
        if !s.ldiurn {
            s.u0 = s.gmu;
            s.sza = s.gmu.acos() * 57.29578;
        }
    }

    trans
}

// ── SCATTR — Eddington 2-stream multiple scattering ──────────────────────────

/// Compute diffuse field and add to FLTAU.
/// Fortran: SUBROUTINE SCATTR
fn scattr(s: &mut ModelState) {
    const NWT: usize = 3;
    // Gauss-Legendre quadrature points and weights
    const TTU: [f64; 3] = [0.1127016654, 0.5000000000, 0.8872983346];
    const TTU2: [f64; 3] = [0.0127016654, 0.2500000000, 0.7872983346];
    const TTW: [f64; 3] = [0.2777777778, 0.4444444444, 0.2777777778];
    // Scattering matrix (Rayleigh phase function, column-major in Fortran)
    const TTP: [[f64; 3]; 3] = [
        [0.3099042361, 0.4578040972, 0.2322916667],
        [0.2861275609, 0.4479166666, 0.2659557725],
        [0.2322916667, 0.4255292360, 0.3421790972],
    ];

    let ntt = s.ntt;
    let u0 = s.u0;
    let rflect = s.rflect;

    // Phase function source from direct beam
    let mut ttp0 = [0.0f64; 3];
    for j in 0..NWT {
        ttp0[j] = 0.375 * (3.0 - TTU2[j] - u0 * u0 + 3.0 * TTU2[j] * u0 * u0);
    }

    // Storage for the specific intensities and their propagators
    let mut ttc = vec![[0.0f64; 3]; ntt + 1];
    let mut ttd = vec![[[0.0f64; 3]; 3]; ntt + 1];

    // ── Top boundary (II=1) ───────────────────────────────────────────────────
    let dtau_top = s.tau[1] - s.tau[0];
    let mut b = [[0.0f64; 3]; 3];
    let mut ttq = [0.0f64; 3];

    for j in 0..NWT {
        let atemp = TTU[j] / dtau_top;
        let btemp = 0.5 / atemp;
        ttq[j] = btemp * 0.25 * s.fltau[0] * (s.piray[0] * ttp0[j] + s.piaer[0]);
        for k in 0..NWT {
            b[j][k] = -btemp * (s.piray[0] * TTP[j][k] + s.piaer[0] * TTW[k]);
        }
        b[j][j] += 1.0 + atemp + btemp;
    }
    matin3(&mut b);
    for k in 0..NWT {
        let atemp = TTU[k] / dtau_top;
        ttc[0][k] = 0.0;
        for j in 0..NWT {
            ttd[0][j][k] = b[j][k] * atemp;
            ttc[0][k] += b[k][j] * ttq[j];
        }
    }

    // ── Interior levels (II=2..NTT-2) ─────────────────────────────────────────
    for ii in 1..ntt - 1 {
        let dtemp = 2.0 / (s.tau[ii + 1] - s.tau[ii - 1]);
        let dta = dtemp / (s.tau[ii] - s.tau[ii - 1]);
        let dtc = dtemp / (s.tau[ii + 1] - s.tau[ii]);
        let mut b = [[0.0f64; 3]; 3];
        let mut ttq2 = [0.0f64; 3];
        for j in 0..NWT {
            let atemp = TTU2[j] * dta;
            let ctemp = TTU2[j] * dtc;
            ttq2[j] =
                0.25 * s.fltau[ii] * (s.piray[ii] * ttp0[j] + s.piaer[ii]) + atemp * ttc[ii - 1][j];
            for k in 0..NWT {
                b[j][k] =
                    -atemp * ttd[ii - 1][j][k] - s.piray[ii] * TTP[j][k] - s.piaer[ii] * TTW[k];
            }
            b[j][j] += 1.0 + atemp + ctemp;
        }
        matin3(&mut b);
        for j in 0..NWT {
            let ctemp = TTU2[j] * dtc;
            ttc[ii][j] = 0.0;
            for k in 0..NWT {
                ttc[ii][j] += b[j][k] * ttq2[k];
                ttd[ii][k][j] = b[k][j] * ctemp;
            }
        }
    }

    // ── Bottom boundary ────────────────────────────────────────────────────────
    let ii = ntt - 1;
    let dtau_bot = s.tau[ii] - s.tau[ii - 1];
    let factor = 4.0 * rflect / (1.0 + rflect);
    let ctemp_b = factor * 0.25 * s.ttfbot;
    let mut b = [[0.0f64; 3]; 3];
    let mut ttq3 = [0.0f64; 3];

    for j in 0..NWT {
        let atemp = TTU[j] / dtau_bot;
        let btemp = 0.5 / atemp;
        ttq3[j] = atemp * ttc[ii - 1][j]
            + ctemp_b
            + btemp * 0.25 * s.fltau[ii] * (s.piray[ii] * ttp0[j] + s.piaer[ii]);
        for k in 0..NWT {
            b[j][k] = -atemp * ttd[ii - 1][j][k]
                - btemp * s.piray[ii] * TTP[j][k]
                - btemp * s.piaer[ii] * TTW[k]
                - factor * TTU[k] * TTW[k];
        }
        b[j][j] += atemp + btemp + 1.0;
    }
    matin3(&mut b);
    ttc[ii] = [0.0; 3];
    for j in 0..NWT {
        for k in 0..NWT {
            ttc[ii][j] += b[j][k] * ttq3[k];
        }
    }

    // ── Back-substitution ─────────────────────────────────────────────────────
    for iii in 1..ntt {
        let idx = ntt - 1 - iii;
        let ttc_next = ttc[idx + 1];
        for j in 0..NWT {
            for k in 0..NWT {
                ttc[idx][j] += ttd[idx][j][k] * ttc_next[k];
            }
        }
    }

    // ── Add diffuse field to direct beam ──────────────────────────────────────
    for i in 0..ntt {
        let mut diffuse = 0.0;
        for j in 0..NWT {
            diffuse += 4.0 * TTW[j] * ttc[i][j];
        }
        s.fltau[i] += diffuse;
    }
}

// ── MATIN3 — 3×3 LU matrix inversion ─────────────────────────────────────────

/// In-place 3×3 LU inversion. Row-major: b[row][col].
/// Fortran: SUBROUTINE MATIN3 (uses column-major B(3,3))
fn matin3(b: &mut [[f64; 3]; 3]) {
    // LU decomposition
    b[1][0] /= b[0][0];
    b[1][1] -= b[1][0] * b[0][1];
    b[1][2] -= b[1][0] * b[0][2];
    b[2][0] /= b[0][0];
    b[2][1] = (b[2][1] - b[2][0] * b[0][1]) / b[1][1];
    b[2][2] -= b[2][0] * b[0][2] + b[2][1] * b[1][2];

    // Invert L
    b[2][1] = -b[2][1];
    b[2][0] = -b[2][0] - b[2][1] * b[1][0];
    b[1][0] = -b[1][0];

    // Invert U
    b[2][2] = 1.0 / b[2][2];
    b[1][2] = -b[1][2] * b[2][2] / b[1][1];
    b[1][1] = 1.0 / b[1][1];
    b[0][2] = -(b[0][1] * b[1][2] + b[0][2] * b[2][2]) / b[0][0];
    b[0][1] = -b[0][1] * b[1][1] / b[0][0];
    b[0][0] = 1.0 / b[0][0];

    // Multiply U⁻¹ × L⁻¹
    b[0][0] += b[0][1] * b[1][0] + b[0][2] * b[2][0];
    b[0][1] += b[0][2] * b[2][1];
    b[1][0] = b[1][1] * b[1][0] + b[1][2] * b[2][0];
    b[1][1] += b[1][2] * b[2][1];
    b[2][0] *= b[2][2];
    b[2][1] *= b[2][2];
}

// ── SPHERE — spherical geometry air mass weights ──────────────────────────────

/// Fill s.wtau[layer][level] with slant-path weights.
/// Fortran: SUBROUTINE SPHERE
pub fn sphere(s: &mut ModelState) {
    let nc = s.nc;
    let rad = s.rad;
    let u0 = s.u0;

    // Radii at each level
    let mut rz = [0.0f64; NL];
    let mut rq = [0.0f64; NL]; // (RZ(II-1)/RZ(II))^2
    for ii in 0..nc {
        rz[ii] = rad + s.z[ii];
    }
    for ii in 1..nc {
        rq[ii - 1] = (rz[ii - 1] / rz[ii]).powi(2);
    }

    let nlbatm = s.nlbatm.saturating_sub(1); // 0-based

    // Tangent height
    s.tanht = if u0 < 0.0 {
        rz[nlbatm] / (1.0 - u0 * u0).sqrt()
    } else {
        rz[nlbatm]
    };
    let tanht = s.tanht;

    // Zero WTAU
    let nc = s.nc;
    for j in 0..nc {
        for k in 0..nc {
            s.wtau[[k, j]] = 0.0;
        }
    }

    // Weighted slant path from each level J to the sun
    for j in 0..nc {
        if rz[j] < tanht {
            continue;
        }

        let mut xmu1 = u0.abs();

        // Upward path to top
        for i in j + 1..nc {
            let xmu2 = (1.0 - rq[i - 1] * (1.0 - xmu1 * xmu1)).max(0.0).sqrt();
            let xl = rz[i] * xmu2 - rz[i - 1] * xmu1;
            s.wtau[[i - 1, j]] += xl * 0.5;
            s.wtau[[i, j]] += xl * 0.5;
            xmu1 = xmu2;
        }

        // Scale height at top
        let xmu2_top = xmu1;
        let airmas_top = airmas(xmu2_top, s.zzht / rad);
        s.wtau[[nc - 1, j]] += s.zzht * airmas_top;

        if u0 >= 0.0 {
            continue;
        }

        // Twilight: downward path below tangent height
        let mut xmu1 = u0.abs();
        let mut ji = j as i64;
        while ji >= 1 {
            let jj = ji as usize;
            let diff = rz[jj] * (1.0 - xmu1 * xmu1).sqrt() - rz[jj - 1];
            if diff < 0.0 {
                let xmu2 = (1.0 - (1.0 - xmu1 * xmu1) / rq[jj - 1]).max(0.0).sqrt();
                let xl = (rz[jj] * xmu1 - rz[jj - 1] * xmu2).abs();
                s.wtau[[jj, j]] += xl;
                s.wtau[[jj - 1, j]] += xl;
                xmu1 = xmu2;
            } else {
                let xl = rz[jj] * xmu1 * 2.0;
                let wting = diff / (rz[jj] - rz[jj - 1]);
                s.wtau[[jj, j]] += xl * 0.5 * (1.0 + wting);
                s.wtau[[jj - 1, j]] += xl * 0.5 * (1.0 - wting);
                break;
            }
            ji -= 1;
        }
    }
}

/// Air mass function: Chapman function approximation.
/// Fortran: AIRMAS(U,H) statement function in SPHERE
#[inline(always)]
fn airmas(u: f64, h: f64) -> f64 {
    (1.0 + h)
        / (u * u
            + 2.0
                * h
                * (1.0
                    - 0.6817 * (-57.3 * u.abs() / (1.0 + 5500.0 * h).sqrt()).exp()
                        / (1.0 + 0.625 * h)))
            .sqrt()
}

// ── Cross-section helper functions ───────────────────────────────────────────

/// O3 cross-section via 3-point T interpolation.
/// Fortran: FUNCTION XSECO3(K, TTT)
fn xseco3(k: usize, ttt: f64, s: &ModelState) -> f64 {
    flint(
        ttt,
        s.tqq[[0, 2]],
        s.tqq[[1, 2]],
        s.tqq[[2, 2]],
        s.qo3[[k, 0]],
        s.qo3[[k, 1]],
        s.qo3[[k, 2]],
    )
}

/// O3→O(1D) quantum yield via 3-point T interpolation.
/// Fortran: FUNCTION XSEC1D(K, TTT)
fn xsec1d(k: usize, ttt: f64, s: &ModelState) -> f64 {
    flint(
        ttt,
        s.tqq[[0, 3]],
        s.tqq[[1, 3]],
        s.tqq[[2, 3]],
        s.q1d[[k, 0]],
        s.q1d[[k, 1]],
        s.q1d[[k, 2]],
    )
}

/// O2 cross-section. ksr=0-based SRB band index (or 0 for non-SRB).
/// Fortran: FUNCTION XSECO2(K, TTT, KSR, KODF)
fn xseco2(k: usize, ttt: f64, ksr: usize, kodf: usize, is_srb: bool, s: &ModelState) -> f64 {
    if !is_srb {
        flint(
            ttt,
            s.tqq[[0, 1]],
            s.tqq[[1, 1]],
            s.tqq[[2, 1]],
            s.qo2[[k, 0]],
            s.qo2[[k, 1]],
            s.qo2[[k, 2]],
        )
    } else {
        flint(
            ttt,
            s.tqq[[0, 1]],
            s.tqq[[1, 1]],
            s.tqq[[2, 1]],
            s.o2x[[kodf, ksr, 0]],
            s.o2x[[kodf, ksr, 1]],
            s.o2x[[kodf, ksr, 2]],
        )
    }
}

/// 3-point linear interpolation.
/// Fortran: FUNCTION FLINT(TINT, T1, T2, T3, F1, F2, F3)
pub fn flint(tint: f64, t1: f64, t2: f64, t3: f64, f1: f64, f2: f64, f3: f64) -> f64 {
    if tint <= t2 {
        if tint <= t1 {
            f1
        } else {
            f1 + (f2 - f1) * (tint - t1) / (t2 - t1)
        }
    } else if tint >= t3 {
        f3
    } else {
        f2 + (f3 - f2) * (tint - t2) / (t3 - t2)
    }
}

/// Rayleigh cross-section (cm²).
/// Fortran: FUNCTION RAYLAY(WAVE) — WAVE in nm
pub fn raylay(wave: f64) -> f64 {
    if wave < 170.0 {
        1.0e-24
    } else {
        let wsqi = 1.0e6 / (wave * wave);
        let refrm1 = 1.0e-6 * (64.328 + 29498.1 / (146.0 - wsqi) + 255.4 / (41.0 - wsqi));
        5.40e-21 * (refrm1 * wsqi).powi(2)
    }
}

/// Mie aerosol cross-section (relative, normalized at 300 nm).
/// Fortran: FUNCTION SCTMIE(WAVE)
pub fn sctmie(wave: f64, _turbdx: f64) -> f64 {
    (wave / 300.0_f64).powf(-0.11)
}
