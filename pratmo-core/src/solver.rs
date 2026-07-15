// butil.f → solver module
// NEWRAF, NEWRAX, LINSLV, RESOLV, FIXRAT, FIXMIX, RPLACE, SPLACE

use anyhow::Result;
use std::sync::OnceLock;

use crate::chemistry::{chems, rhslhs_jacobian, rhslhs_rhs};
use crate::constants::NDEN;
use crate::state::ModelState;

// ── RPLACE / SPLACE — scatter/gather species arrays ──────────────────────────

/// Load XN[0..NDEN] from named density arrays at box j (0-based).
/// Fortran: SUBROUTINE RPLACE(XN, J)
pub fn rplace(s: &ModelState, xn: &mut [f64; NDEN], j: usize) {
    let n = s.n;
    xn.iter_mut().for_each(|v| *v = 1.0e-36);

    let put = |xn: &mut [f64; NDEN], ni: usize, val: f64| {
        if (1..=NDEN).contains(&ni) {
            xn[ni - 1] = val;
        }
    };

    put(xn, n[0], s.dno[j]);
    put(xn, n[1], s.dno2[j]);
    put(xn, n[2], s.dno3[j]);
    put(xn, n[3], s.dn2o5[j]);
    put(xn, n[4], s.dhno3[j]);
    put(xn, n[14], s.dhno2[j]);
    put(xn, n[20], s.dhno4[j]);
    put(xn, n[5], s.dh[j]);
    put(xn, n[6], s.doh[j]);
    put(xn, n[7], s.dho2[j]);
    put(xn, n[8], s.dh2o2[j]);
    put(xn, n[9], s.do_[j]);
    put(xn, n[10], s.do3[j]);
    put(xn, n[15], s.dhcl[j]);
    put(xn, n[16], s.dcl[j]);
    put(xn, n[17], s.dcl2[j]);
    put(xn, n[18], s.dclo[j]);
    put(xn, n[19], s.dclno3[j]);
    put(xn, n[21], s.dhocl[j]);
    put(xn, n[27], s.doclo[j]);
    put(xn, n[28], s.dcl2o2[j]);
    put(xn, n[11], s.dbro[j]);
    put(xn, n[12], s.dbr[j]);
    put(xn, n[13], s.dhbr[j]);
    put(xn, n[22], s.dbrno3[j]);
    put(xn, n[23], s.dhobr[j]);
    put(xn, n[24], s.dh2co[j]);
    put(xn, n[25], s.droo[j]);
    put(xn, n[26], s.drooh[j]);
    put(xn, n[29], s.dbrcl[j]);
    if s.liod {
        put(xn, n[30], s.di_[j]);
        put(xn, n[31], s.dio[j]);
        put(xn, n[32], s.dhoi[j]);
        put(xn, n[33], s.diono2[j]);
        put(xn, n[34], s.dhi[j]);
        put(xn, n[35], s.doio[j]);
        put(xn, n[36], s.di2[j]);
        put(xn, n[37], s.di2o2[j]);
        put(xn, n[38], s.di2o3[j]);
        put(xn, n[39], s.di2o4[j]);
    }
}

/// Store XN[0..NDEN] back into named density arrays at box j (0-based).
/// Fortran: SUBROUTINE SPLACE(XN, J)
pub fn splace(s: &mut ModelState, xn: &[f64; NDEN], j: usize) {
    let n = s.n;
    let get = |xn: &[f64; NDEN], ni: usize| -> f64 {
        if (1..=NDEN).contains(&ni) {
            xn[ni - 1]
        } else {
            0.0
        }
    };

    s.dno[j] = get(xn, n[0]);
    s.dno2[j] = get(xn, n[1]);
    s.dno3[j] = get(xn, n[2]);
    s.dn2o5[j] = get(xn, n[3]);
    s.dhno3[j] = get(xn, n[4]);
    s.dhno2[j] = get(xn, n[14]);
    s.dhno4[j] = get(xn, n[20]);
    s.dh[j] = get(xn, n[5]);
    s.doh[j] = get(xn, n[6]);
    s.dho2[j] = get(xn, n[7]);
    s.dh2o2[j] = get(xn, n[8]);
    s.do_[j] = get(xn, n[9]);
    s.do3[j] = get(xn, n[10]);
    s.dhcl[j] = get(xn, n[15]);
    s.dcl[j] = get(xn, n[16]);
    s.dcl2[j] = get(xn, n[17]);
    s.dclo[j] = get(xn, n[18]);
    s.dclno3[j] = get(xn, n[19]);
    s.dhocl[j] = get(xn, n[21]);
    s.doclo[j] = get(xn, n[27]);
    s.dcl2o2[j] = get(xn, n[28]);
    s.dbro[j] = get(xn, n[11]);
    s.dbr[j] = get(xn, n[12]);
    s.dhbr[j] = get(xn, n[13]);
    s.dbrno3[j] = get(xn, n[22]);
    s.dhobr[j] = get(xn, n[23]);
    s.dh2co[j] = get(xn, n[24]);
    s.droo[j] = get(xn, n[25]);
    s.drooh[j] = get(xn, n[26]);
    s.dbrcl[j] = get(xn, n[29]);
    if s.liod {
        s.di_[j] = get(xn, n[30]);
        s.dio[j] = get(xn, n[31]);
        s.dhoi[j] = get(xn, n[32]);
        s.diono2[j] = get(xn, n[33]);
        s.dhi[j] = get(xn, n[34]);
        s.doio[j] = get(xn, n[35]);
        s.di2[j] = get(xn, n[36]);
        s.di2o2[j] = get(xn, n[37]);
        s.di2o3[j] = get(xn, n[38]);
        s.di2o4[j] = get(xn, n[39]);
    }
}

// ── FIXRAT — family conservation cubic solve ──────────────────────────────────

fn fix_iodine_family(x: &mut [f64; NDEN], s: &ModelState, target: f64) {
    if !s.liod {
        return;
    }
    let weighted_species = [
        (30usize, 1.0),
        (31, 1.0),
        (32, 1.0),
        (33, 1.0),
        (34, 1.0),
        (35, 1.0),
        (36, 2.0),
        (37, 2.0),
        (38, 2.0),
        (39, 2.0),
    ];
    let mut total = 0.0;
    for &(species, weight) in &weighted_species {
        let slot = s.n[species];
        if slot > 0 && slot <= s.ntot {
            total += weight * x[slot - 1];
        }
    }
    if target <= 0.0 {
        for &(species, _) in &weighted_species {
            let slot = s.n[species];
            if slot > 0 && slot <= s.ntot {
                x[slot - 1] = 0.0;
            }
        }
    } else if total > 1.0e-30 {
        let scale = target / total;
        for &(species, _) in &weighted_species {
            let slot = s.n[species];
            if slot > 0 && slot <= s.ntot {
                x[slot - 1] *= scale;
            }
        }
    }
}

/// Rescale species so NOy, Cly, Bry, and gas-phase Iy families match targets.
/// Solves a cubic equation (regula falsi + Newton-Raphson) for the NOy scale factor.
/// Fortran: SUBROUTINE FIXRAT(X, I) — I is the box index (1-based in Fortran)
pub fn fixrat(x: &mut [f64; NDEN], s: &ModelState, ib: usize) {
    let n = s.n;
    let dm = s.dm[s.ialt];

    let totnoy = s.fnoy[ib] * dm;
    let totclx = s.fclx[ib] * dm;
    let totbrx = s.fbrx[ib] * dm;
    let totiyx = s.fiodx[ib] * dm;

    // Apply the prescribed gas-phase Iy constraint independently of the legacy
    // NOy/Cly/Bry cubic. The legacy NOy definition does not include trace IONO2.
    fix_iodine_family(x, s, totiyx);

    // BrCl cap
    let n30 = n[29];
    let brcl = x[n30 - 1].min(totclx / 2.0).min(totbrx / 2.0);
    x[n30 - 1] = brcl;

    let totclx = totclx - brcl;
    let totbrx = totbrx - brcl;

    // Family sums (without BrCl)
    let xxn = x[n[0] - 1]
        + x[n[1] - 1]
        + x[n[2] - 1]
        + x[n[3] - 1] * 2.0
        + x[n[4] - 1]
        + x[n[14] - 1]
        + x[n[20] - 1];
    let xxc = x[n[15] - 1]
        + x[n[16] - 1]
        + x[n[18] - 1]
        + x[n[21] - 1]
        + x[n[27] - 1]
        + 2.0 * (x[n[17] - 1] + x[n[28] - 1]);
    let xxcn = x[n[19] - 1]; // ClONO2 — links Cly and NOy
    let xxb = x[n[11] - 1] + x[n[12] - 1] + x[n[13] - 1] + x[n[23] - 1];
    let xxbn = x[n[22] - 1]; // BrONO2 — links Bry and NOy

    if xxn < 1e-30 || xxcn < 1e-30 || xxbn < 1e-30 {
        return;
    }

    // Cubic coefficients for RNOY (NOy scale factor)
    let xxxc = xxc / xxcn;
    let xxxb = xxb / xxbn;
    let tempa = xxxc + xxxb + (totclx + totbrx - totnoy) / xxn;
    let tempb = xxxc * (totbrx - totnoy) / xxn + xxxb * (totclx - totnoy) / xxn + xxxc * xxxb;
    let tempc = -xxxb * xxxc * totnoy / xxn;
    let tempa3 = tempa / 3.0;

    // Closure functions
    let fcube = |y: f64| tempc + y * (tempb + y * (tempa + y));
    let fcube1 = |y: f64| tempb + y * (2.0 * tempa + 3.0 * y);
    let fcube2 = |y: f64| 2.0 * tempa + 6.0 * y;
    let fcuba =
        |y: f64| tempc.abs() + (tempb * y).abs() + (tempa * y * y).abs() + (y * y * y).abs();

    const EPSM: f64 = 1.0e-12;
    const EPS2: f64 = 1.0e-6;

    // Locate extrema
    let discrm = tempa * tempa - 3.0 * tempb;
    let (x1, x2, mut x_bar);
    if discrm <= 0.0 {
        x1 = -tempa3;
        x2 = x1;
        x_bar = x1 + 1.0;
    } else {
        let xx1_raw = -1.5 * tempb / (tempa * tempa);
        let xx1 = if xx1_raw > EPS2 {
            (1.0 + 2.0 * xx1_raw).sqrt() - 1.0
        } else {
            xx1_raw
        };
        let (mut xa, mut xb) = (xx1 * tempa3, (-xx1 - 2.0) * tempa3);
        if tempa3 < 0.0 {
            std::mem::swap(&mut xa, &mut xb);
        }
        x1 = xa;
        x2 = xb;
        x_bar = if fcube(x1) > 0.0 { x2 - 1.0 } else { x1 + 1.0 };
    }

    // Regula falsi to bracket root
    let f1 = fcube(x1);
    let l_up = f1 > 0.0;
    let mut xlo = 0.0_f64;
    let mut xhi = 0.0_f64;
    let mut l_lo = false;
    let mut l_hi = false;

    for _ in 0..20 {
        let fbar = fcube(x_bar);
        if fbar <= 0.0 {
            xlo = x_bar;
            l_lo = true;
            if l_hi {
                break;
            }
            x_bar = if !l_up {
                x1 + (xlo - x1) * 10.0
            } else {
                x2 - (x2 - xlo) * 0.1
            };
        } else {
            xhi = x_bar;
            l_hi = true;
            if l_lo {
                break;
            }
            x_bar = if !l_up {
                x1 + (xhi - x1) * 0.1
            } else {
                x2 - (x2 - xhi) * 10.0
            };
        }
    }

    if !l_lo || !l_hi {
        return;
    } // failed to bracket

    // Newton-Raphson 2nd-order refinement within [xlo, xhi]
    let mut x0 = xhi;
    for _ in 0..40 {
        let f0 = fcube(x0);
        if f0.abs() < EPSM * fcuba(x0) {
            break;
        }
        let f1v = fcube1(x0);
        let f0f1 = f0 / f1v;
        let delx = -f0f1 * (1.0 + 0.5 * f0f1 * fcube2(x0) / f1v);
        if delx.abs() < EPSM * x0.abs() {
            break;
        }
        x0 = xhi.min(xlo.max(x0 + delx));
    }

    if x0 < 1.0e-30 {
        return;
    }

    let rnoy = x0;
    let rbrx = totbrx / (xxb + x0 * xxbn);
    let rclx = totclx / (xxc + x0 * xxcn);

    // Rescale species
    for ni in [n[0], n[1], n[2], n[3], n[4], n[14], n[20]] {
        x[ni - 1] *= rnoy;
    }
    for ni in [n[15], n[16], n[17], n[18], n[21], n[27], n[28]] {
        x[ni - 1] *= rclx;
    }
    x[n[19] - 1] *= rclx * rnoy; // ClONO2
    for ni in [n[11], n[12], n[13], n[23]] {
        x[ni - 1] *= rbrx;
    }
    x[n[22] - 1] *= rbrx * rnoy; // BrONO2
}

/// Wrapper: rplace → fixrat → splace → update FO3.
/// Fortran: SUBROUTINE FIXMIX
pub fn fixmix(s: &mut ModelState) {
    let ib = s.ibox;
    let ialt = s.ialt;
    let mut xn = [0.0f64; NDEN];
    rplace(s, &mut xn, ib);
    fixrat(&mut xn, s, ib);
    splace(s, &xn, ib);
    let dm = s.dm[ialt];
    s.fo3[ib] = s.do3[ib] / dm;
}

// ── NEWRAF — adaptive time-step driver ────────────────────────────────────────

/// Newton-Raphson driver: halves DELTT until NEWRAX converges, then rebuilds.
/// Fortran: SUBROUTINE NEWRAF(DAMP1, X, N, NDXRAF)
pub fn newraf(s: &mut ModelState, damp1: f64, x: &mut [f64; NDEN], n: usize) -> Result<()> {
    let code = newrax(s, damp1, x, n);
    if code == 0 || s.deltt < 1.0e-20 {
        return Ok(());
    }
    retry_probe_newraf(s, "initial", code, 0);

    // Time-step halving loop
    let deltt0 = s.deltt;
    let mut kut = 0;
    for k in 0..30 {
        s.deltt *= 2.0;
        let c = newrax(s, damp1, x, n);
        if c == 0 {
            kut = k + 1;
            retry_probe_newraf(s, "cut_success", c, kut);
            break;
        }
        retry_probe_newraf(s, "cut_retry", c, k + 1);
        if k == 29 {
            s.lprts = true;
            let final_code = newrax(s, damp1, x, n);
            s.deltt = deltt0;
            if final_code != 0 {
                s.newraf_nonconvergence_count += 1;
            }
            return Ok(());
        }
    }

    // Rebuild back to original time step
    for _ in 0..kut {
        s.deltt /= 2.0;
        let c = newrax(s, damp1, x, n);
        if c != 0 {
            retry_probe_newraf(s, "rebuild_fail", c, kut);
            s.lprts = true;
            let final_code = newrax(s, damp1, x, n);
            s.deltt = deltt0;
            if final_code != 0 {
                s.newraf_nonconvergence_count += 1;
            }
            return Ok(());
        }
    }
    s.deltt = deltt0;
    Ok(())
}

// ── NEWRAX — inner Newton-Raphson ────────────────────────────────────────────

/// Inner NR loop. Returns 0=converged, 1=negative density, 2=non-convergence.
/// Fortran: SUBROUTINE NEWRAX(DAMP1, X, N, NDXRAF)
pub fn newrax(s: &mut ModelState, damp1: f64, x: &mut [f64; NDEN], n: usize) -> i32 {
    const NUMITR: usize = 18;
    s.radcount += 1.0;

    let mut errxo = 1.0_f64;
    let mut last_errpl = 0.0_f64;
    let mut last_errpl_j = 0usize;

    // Load XR from X
    for j in 0..s.ntotx {
        s.xr[j] = x[j];
    }

    for iter in 0..NUMITR {
        s.raxloop += 1.0;
        chems(s);

        if errxo < s.raferr {
            break;
        }
        if iter == NUMITR - 1 {
            retry_probe_newrax_nonconverged(s, iter + 1, last_errpl, last_errpl_j, errxo);
            return 2;
        }

        // RHS
        let fxo = rhslhs_rhs(s);

        // Convergence check on P-L residual
        let ntot = s.ntot;
        let mut errpl = 0.0_f64;
        let mut errpl_j = 0usize;
        for j in 0..n.min(ntot) {
            let denom = s.rp[j] + s.rl[j];
            if denom > 0.0 {
                let err = (fxo[j] / denom).abs();
                if err >= errpl {
                    errpl = err;
                    errpl_j = j;
                }
            }
        }
        last_errpl = errpl;
        last_errpl_j = errpl_j;

        // A residual is formed from several additions/subtractions.  In a
        // cancellation-heavy dark-halogen case, the requested Fortran
        // tolerance can be below the relative error that f64 arithmetic can
        // resolve.  Normal Rust mode raises the threshold to a conservative
        // operation-scale roundoff floor; parity mode deliberately keeps the
        // original Fortran criterion so differential output is unchanged.
        #[cfg(not(feature = "fortran-parity"))]
        let achievable_rafpml = (0..n.min(ntot)).fold(s.rafpml, |tol, j| {
            tol.max(achievable_errpl(
                s.xr[j], s.xnold[j], s.deltt, s.rp[j], s.rl[j],
            ))
        });
        #[cfg(feature = "fortran-parity")]
        let achievable_rafpml = s.rafpml;

        if errpl < achievable_rafpml {
            break;
        }

        // Jacobian + solve
        let mut xoo = [0.0f64; NDEN];
        if s.lsvjac {
            resolv_in_place(s, &fxo[..n], &mut xoo[..n]);
        } else {
            rhslhs_jacobian(s);
            linslv(s, &fxo[..n], &mut xoo[..n], n);
        }

        // Apply correction
        errxo = 0.0;
        for j in 0..n.min(ntot) {
            let dx = if iter == 0 { damp1 * xoo[j] } else { xoo[j] };
            if dx.abs() > 1.0e-16 {
                let rel = (dx / s.xr[j]).abs();
                errxo = errxo.max(rel);
            }
            let mut temp = s.xr[j] + dx;
            if iter < 6 {
                temp = temp.max(s.rafmin * s.xr[j]).min(s.rafmax * s.xr[j]);
            }
            s.xr[j] = temp;
        }

        // Check for negative densities
        if let Some(j) = (0..n.min(ntot)).find(|&j| s.xr[j] <= 0.0) {
            retry_probe_newrax_negative(s, iter + 1, j);
            return 1;
        }
    }

    // Success: write XR back to X
    for j in 0..s.ntotx {
        x[j] = s.xr[j];
    }
    0
}

/// Estimate the smallest trustworthy relative P-L residual for one species.
///
/// `RHSLHS` evaluates `RL - RP + (XR-XNOLD)*DELTT`; when those terms nearly
/// cancel, comparing the result with `RAFPML` alone asks for more precision
/// than f64 can provide.  The factor of eight is a conservative bound for the
/// handful of rounded operations in that expression.  This helper is compiled
/// only for normal Rust mode: `fortran-parity` must retain the original
/// executable's strict `RAFPML` test.
#[cfg(not(feature = "fortran-parity"))]
fn achievable_errpl(xr: f64, xnold: f64, deltt: f64, rp: f64, rl: f64) -> f64 {
    let denom = rp + rl;
    if denom <= 0.0 {
        return 0.0;
    }
    let operation_scale = rp.abs() + rl.abs() + ((xr - xnold) * deltt).abs();
    8.0 * f64::EPSILON * operation_scale / denom
}

fn retry_probe_enabled() -> bool {
    static ENABLED: OnceLock<bool> = OnceLock::new();
    *ENABLED.get_or_init(|| {
        std::env::var("PRATMO_RETRY_PROBE")
            .map(|v| {
                matches!(
                    v.as_str(),
                    "1" | "true" | "TRUE" | "yes" | "YES" | "on" | "ON"
                )
            })
            .unwrap_or(false)
    })
}

fn retry_probe_hhmm(s: &ModelState) -> i32 {
    if s.ittt <= 0 {
        return 0;
    }
    let idx = (s.ittt as usize).saturating_sub(1);
    s.nhhmm.get(idx).copied().unwrap_or(0)
}

fn retry_probe_species(s: &ModelState, j: usize) -> &str {
    let name = s.tname[j].trim();
    if !name.is_empty() {
        return name;
    }
    let name = s.tnamet[j].trim();
    if !name.is_empty() {
        return name;
    }
    "UNKNOWN"
}

fn retry_probe_slot_value(s: &ModelState, species: usize) -> f64 {
    let slot = s.n[species];
    if (1..=NDEN).contains(&slot) {
        s.xr[slot - 1]
    } else {
        f64::NAN
    }
}

fn retry_probe_newraf(s: &ModelState, phase: &str, code: i32, cut_depth: usize) {
    if !retry_probe_enabled() {
        return;
    }
    eprintln!(
        "PRATMO_RETRY newraf phase={} code={} cut_depth={} box={} alt={} ittt={} hhmm={} deltt={:.6e} gmu={:.6e} no2={:.6e} bro={:.6e}",
        phase,
        code,
        cut_depth,
        s.ibox + 1,
        s.ialt + 1,
        s.ittt,
        retry_probe_hhmm(s),
        s.deltt,
        s.gmu,
        retry_probe_slot_value(s, 1),
        retry_probe_slot_value(s, 11),
    );
}

fn retry_probe_newrax_negative(s: &ModelState, iter: usize, j: usize) {
    if !retry_probe_enabled() {
        return;
    }
    eprintln!(
        "PRATMO_RETRY newrax code=1 reason=negative box={} alt={} ittt={} hhmm={} deltt={:.6e} gmu={:.6e} iter={} slot={} species={} value={:.6e} no2={:.6e} bro={:.6e}",
        s.ibox + 1,
        s.ialt + 1,
        s.ittt,
        retry_probe_hhmm(s),
        s.deltt,
        s.gmu,
        iter,
        j + 1,
        retry_probe_species(s, j),
        s.xr[j],
        retry_probe_slot_value(s, 1),
        retry_probe_slot_value(s, 11),
    );
}

fn retry_probe_newrax_nonconverged(
    s: &ModelState,
    iter: usize,
    errpl: f64,
    errpl_j: usize,
    errxo: f64,
) {
    if !retry_probe_enabled() {
        return;
    }
    eprintln!(
        "PRATMO_RETRY newrax code=2 reason=nonconverged box={} alt={} ittt={} hhmm={} deltt={:.6e} gmu={:.6e} iter={} errpl={:.6e} errpl_slot={} errpl_species={} errxo={:.6e} no2={:.6e} bro={:.6e}",
        s.ibox + 1,
        s.ialt + 1,
        s.ittt,
        retry_probe_hhmm(s),
        s.deltt,
        s.gmu,
        iter,
        errpl,
        errpl_j + 1,
        retry_probe_species(s, errpl_j),
        errxo,
        retry_probe_slot_value(s, 1),
        retry_probe_slot_value(s, 11),
    );
}

#[cfg(all(test, not(feature = "fortran-parity")))]
mod tests {
    use super::{achievable_errpl, fix_iodine_family};
    use crate::{constants::NDEN, state::ModelState};

    fn identity_species_map() -> [usize; NDEN] {
        let mut n = [0usize; NDEN];
        for (index, slot) in n.iter_mut().enumerate() {
            *slot = index + 1;
        }
        n
    }

    #[test]
    fn precision_floor_catches_cancellation_limited_residual() {
        // A dark BrCl-like species can have a tiny production/loss denominator
        // while the time-step term is much larger.  Asking for 1e-10 relative
        // residual in that case is below the f64 roundoff floor.
        let floor = achievable_errpl(1.0e-4, 0.0, 1.0, 1.0e-10, 1.0e-10);
        assert!(floor > 1.0e-10, "floor={floor:e}");
        assert!(floor.is_finite());
    }

    #[test]
    fn precision_floor_is_zero_without_positive_reaction_scale() {
        assert_eq!(achievable_errpl(1.0, 0.0, 1.0, -1.0, 1.0), 0.0);
    }

    #[test]
    fn iodine_family_scaling_uses_atom_weights_for_every_species() {
        let mut state = ModelState::new();
        state.n = identity_species_map();
        state.ntot = NDEN;
        state.liod = true;

        for (species, iodine_atoms) in [
            (30usize, 1.0),
            (31, 1.0),
            (32, 1.0),
            (33, 1.0),
            (34, 1.0),
            (35, 1.0),
            (36, 2.0),
            (37, 2.0),
            (38, 2.0),
            (39, 2.0),
        ] {
            let mut x = [0.0; NDEN];
            x[state.n[species] - 1] = 1.0;
            fix_iodine_family(&mut x, &state, 10.0);
            assert_eq!(
                x[state.n[species] - 1],
                10.0 / iodine_atoms,
                "wrong iodine-atom weight for species index {species}"
            );
        }
    }

    #[test]
    fn zero_iodine_target_clears_every_iodine_species() {
        let mut state = ModelState::new();
        state.n = identity_species_map();
        state.ntot = NDEN;
        state.liod = true;
        let mut x = [1.0; NDEN];

        fix_iodine_family(&mut x, &state, 0.0);

        for species in 30..=39 {
            assert_eq!(x[state.n[species] - 1], 0.0, "species index {species}");
        }
    }
}

// ── LINSLV — Crout LU decomposition + solve ───────────────────────────────────

/// Partial-pivoting LU decomposition of s.a_mat, then forward/back solve.
/// Modifies s.a_mat in-place and stores pivots in s.ipa.
/// Fortran: SUBROUTINE LINSLV(B, X, N)
pub fn linslv(s: &mut ModelState, b: &[f64], x: &mut [f64], n: usize) {
    if n == 30 {
        linslv_fixed::<30>(s, b, x);
        return;
    }
    if n == 40 {
        linslv_fixed::<40>(s, b, x);
        return;
    }

    // Crout decomposition with partial pivoting on the A matrix
    // Use the Array2 backing storage as a column-major Newton matrix so the
    // Crout column updates walk contiguous memory, matching the Fortran access pattern.
    let mut s_col = [0.0f64; NDEN];
    let a = s.a_mat.as_slice_mut().expect("a_mat is contiguous");

    for kr in 0..n {
        // Copy column KR of A into S
        let kr_col = kr * NDEN;
        s_col[..n].copy_from_slice(&a[kr_col..kr_col + n]);

        // Apply previous eliminations
        if kr > 0 {
            for j in 0..kr {
                let jp = s.ipa[j];
                a[kr_col + j] = s_col[jp];
                s_col[jp] = s_col[j];
                let ajkr = a[kr_col + j];
                let j_col = j * NDEN;
                for i in j + 1..n {
                    s_col[i] -= a[j_col + i] * ajkr;
                }
            }
        }

        // Find pivot
        let mut smax = s_col[kr].abs();
        let mut krmax = kr;
        for i in kr..n {
            if s_col[i].abs() >= smax {
                krmax = i;
                smax = s_col[i].abs();
            }
        }
        s.ipa[kr] = krmax;
        a[kr_col + kr] = s_col[krmax];
        let div = 1.0 / s_col[krmax];
        s_col[krmax] = s_col[kr];

        if kr < n - 1 {
            for i in kr + 1..n {
                a[kr_col + i] = s_col[i] * div;
            }
        }
    }

    resolv_in_place(s, b, x);
}

fn linslv_fixed<const N: usize>(s: &mut ModelState, b: &[f64], x: &mut [f64]) {
    let mut s_col = [0.0f64; NDEN];
    let a = s.a_mat.as_slice_mut().expect("a_mat is contiguous");

    for kr in 0..N {
        let kr_col = kr * NDEN;
        s_col[..N].copy_from_slice(&a[kr_col..kr_col + N]);

        if kr > 0 {
            for j in 0..kr {
                let jp = s.ipa[j];
                a[kr_col + j] = s_col[jp];
                s_col[jp] = s_col[j];
                let ajkr = a[kr_col + j];
                let j_col = j * NDEN;
                for i in j + 1..N {
                    s_col[i] -= a[j_col + i] * ajkr;
                }
            }
        }

        let mut smax = s_col[kr].abs();
        let mut krmax = kr;
        for i in kr..N {
            if s_col[i].abs() >= smax {
                krmax = i;
                smax = s_col[i].abs();
            }
        }
        s.ipa[kr] = krmax;
        a[kr_col + kr] = s_col[krmax];
        let div = 1.0 / s_col[krmax];
        s_col[krmax] = s_col[kr];

        if kr < N - 1 {
            for i in kr + 1..N {
                a[kr_col + i] = s_col[i] * div;
            }
        }
    }

    resolv_fixed::<N>(s, b, x);
}

fn resolv_fixed<const N: usize>(s: &ModelState, b: &[f64], x: &mut [f64]) {
    let mut sv = [0.0f64; NDEN];
    let a = s.a_mat.as_slice().expect("a_mat is contiguous");
    sv[..N].copy_from_slice(&b[..N]);

    for i in 0..N {
        let ip = s.ipa[i];
        x[i] = sv[ip];
        sv[ip] = sv[i];
        if i < N - 1 {
            let i_col = i * NDEN;
            for j in i + 1..N {
                sv[j] -= a[i_col + j] * x[i];
            }
        }
    }

    for i in (0..N).rev() {
        let mut sum = x[i];
        for j in i + 1..N {
            sum -= a[j * NDEN + i] * x[j];
        }
        x[i] = sum / a[i * NDEN + i];
    }
}

// ── RESOLV — back-solve with stored LU ───────────────────────────────────────

/// Solve using the already-decomposed A matrix and pivot array.
/// Fortran: SUBROUTINE RESOLV(B, X, N)
pub fn resolv_in_place(s: &ModelState, b: &[f64], x: &mut [f64]) {
    let n = b.len();
    let mut sv = [0.0f64; NDEN];
    let a = s.a_mat.as_slice().expect("a_mat is contiguous");
    for i in 0..n {
        sv[i] = b[i];
    }

    // Forward substitution with pivoting
    for i in 0..n {
        let ip = s.ipa[i];
        x[i] = sv[ip];
        sv[ip] = sv[i];
        if i < n - 1 {
            let i_col = i * NDEN;
            for j in i + 1..n {
                sv[j] -= a[i_col + j] * x[i];
            }
        }
    }

    // Back substitution
    for i in (0..n).rev() {
        let mut sum = x[i];
        for j in i + 1..n {
            sum -= a[j * NDEN + i] * x[j];
        }
        x[i] = sum / a[i * NDEN + i];
    }
}

// ── RESOLV public wrapper ─────────────────────────────────────────────────────

/// Public back-solve (reuses previously factored A matrix).
/// Fortran: SUBROUTINE RESOLV(B, X, N)
pub fn resolv(s: &ModelState, b: &[f64]) -> Vec<f64> {
    let mut x = vec![0.0f64; b.len()];
    resolv_in_place(s, b, &mut x);
    x
}
