// butil.f → solver module
// NEWRAF, NEWRAX, LINSLV, RESOLV, FIXRAT, FIXMIX, RPLACE, SPLACE

use anyhow::{bail, Result};

use crate::chemistry::{chems, rhslhs};
use crate::constants::NDEN;
use crate::state::ModelState;

// ── RPLACE / SPLACE — scatter/gather species arrays ──────────────────────────

/// Load XN[0..NDEN] from named density arrays at box j (0-based).
/// Fortran: SUBROUTINE RPLACE(XN, J)
pub fn rplace(s: &ModelState, xn: &mut [f64; NDEN], j: usize) {
    let n = s.n;
    xn.iter_mut().for_each(|v| *v = 1.0e-36);

    let put = |xn: &mut [f64; NDEN], ni: usize, val: f64| {
        if ni >= 1 && ni <= NDEN { xn[ni - 1] = val; }
    };

    put(xn, n[0],  s.dno[j]);
    put(xn, n[1],  s.dno2[j]);
    put(xn, n[2],  s.dno3[j]);
    put(xn, n[3],  s.dn2o5[j]);
    put(xn, n[4],  s.dhno3[j]);
    put(xn, n[14], s.dhno2[j]);
    put(xn, n[20], s.dhno4[j]);
    put(xn, n[5],  s.dh[j]);
    put(xn, n[6],  s.doh[j]);
    put(xn, n[7],  s.dho2[j]);
    put(xn, n[8],  s.dh2o2[j]);
    put(xn, n[9],  s.do_[j]);
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
    }
}

/// Store XN[0..NDEN] back into named density arrays at box j (0-based).
/// Fortran: SUBROUTINE SPLACE(XN, J)
pub fn splace(s: &mut ModelState, xn: &[f64; NDEN], j: usize) {
    let n = s.n;
    let get = |xn: &[f64; NDEN], ni: usize| -> f64 {
        if ni >= 1 && ni <= NDEN { xn[ni - 1] } else { 0.0 }
    };

    s.dno[j]    = get(xn, n[0]);
    s.dno2[j]   = get(xn, n[1]);
    s.dno3[j]   = get(xn, n[2]);
    s.dn2o5[j]  = get(xn, n[3]);
    s.dhno3[j]  = get(xn, n[4]);
    s.dhno2[j]  = get(xn, n[14]);
    s.dhno4[j]  = get(xn, n[20]);
    s.dh[j]     = get(xn, n[5]);
    s.doh[j]    = get(xn, n[6]);
    s.dho2[j]   = get(xn, n[7]);
    s.dh2o2[j]  = get(xn, n[8]);
    s.do_[j]    = get(xn, n[9]);
    s.do3[j]    = get(xn, n[10]);
    s.dhcl[j]   = get(xn, n[15]);
    s.dcl[j]    = get(xn, n[16]);
    s.dcl2[j]   = get(xn, n[17]);
    s.dclo[j]   = get(xn, n[18]);
    s.dclno3[j] = get(xn, n[19]);
    s.dhocl[j]  = get(xn, n[21]);
    s.doclo[j]  = get(xn, n[27]);
    s.dcl2o2[j] = get(xn, n[28]);
    s.dbro[j]   = get(xn, n[11]);
    s.dbr[j]    = get(xn, n[12]);
    s.dhbr[j]   = get(xn, n[13]);
    s.dbrno3[j] = get(xn, n[22]);
    s.dhobr[j]  = get(xn, n[23]);
    s.dh2co[j]  = get(xn, n[24]);
    s.droo[j]   = get(xn, n[25]);
    s.drooh[j]  = get(xn, n[26]);
    s.dbrcl[j]  = get(xn, n[29]);
    if s.liod {
        s.di_[j]    = get(xn, n[30]);
        s.dio[j]    = get(xn, n[31]);
        s.dhoi[j]   = get(xn, n[32]);
        s.diono2[j] = get(xn, n[33]);
        s.dhi[j]    = get(xn, n[34]);
    }
}

// ── FIXRAT — family conservation cubic solve ──────────────────────────────────

/// Rescale species X[0..30] so NOy, Cly, Bry families match target totals.
/// Solves a cubic equation (regula falsi + Newton-Raphson) for the NOy scale factor.
/// Fortran: SUBROUTINE FIXRAT(X, I) — I is the box index (1-based in Fortran)
pub fn fixrat(x: &mut [f64; NDEN], s: &ModelState, ib: usize) {
    let n  = s.n;
    let dm = s.dm[s.ialt];

    let totnoy = s.fnoy[ib] * dm;
    let totclx = s.fclx[ib] * dm;
    let totbrx = s.fbrx[ib] * dm;
    let totiyx = s.fiodx[ib] * dm;

    // BrCl cap
    let n30 = n[29];
    let brcl = x[n30 - 1]
        .min(totclx / 2.0)
        .min(totbrx / 2.0);
    x[n30 - 1] = brcl;

    let totclx = totclx - brcl;
    let totbrx = totbrx - brcl;

    // Family sums (without BrCl)
    let xxn  = x[n[0]-1] + x[n[1]-1] + x[n[2]-1] + x[n[3]-1]*2.0
              + x[n[4]-1] + x[n[14]-1] + x[n[20]-1];
    let xxc  = x[n[15]-1] + x[n[16]-1] + x[n[18]-1] + x[n[21]-1]
              + x[n[27]-1] + 2.0*(x[n[17]-1] + x[n[28]-1]);
    let xxcn = x[n[19]-1]; // ClONO2 — links Cly and NOy
    let xxb  = x[n[11]-1] + x[n[12]-1] + x[n[13]-1] + x[n[23]-1];
    let xxbn = x[n[22]-1]; // BrONO2 — links Bry and NOy

    if xxn < 1e-30 || xxcn < 1e-30 || xxbn < 1e-30 {
        return;
    }

    // Cubic coefficients for RNOY (NOy scale factor)
    let xxxc = xxc / xxcn;
    let xxxb = xxb / xxbn;
    let tempa = xxxc + xxxb + (totclx + totbrx - totnoy) / xxn;
    let tempb = xxxc * (totbrx - totnoy) / xxn + xxxb * (totclx - totnoy) / xxn
               + xxxc * xxxb;
    let tempc = -xxxb * xxxc * totnoy / xxn;
    let tempa3 = tempa / 3.0;

    // Closure functions
    let fcube  = |y: f64| tempc + y * (tempb + y * (tempa + y));
    let fcube1 = |y: f64| tempb + y * (2.0 * tempa + 3.0 * y);
    let fcube2 = |y: f64| 2.0 * tempa + 6.0 * y;
    let fcuba  = |y: f64| tempc.abs() + (tempb * y).abs() + (tempa * y * y).abs() + (y * y * y).abs();

    const EPSM: f64 = 1.0e-12;
    const EPS2: f64 = 1.0e-6;

    // Locate extrema
    let discrm = tempa * tempa - 3.0 * tempb;
    let (x1, x2, mut x_bar);
    if discrm <= 0.0 {
        x1   = -tempa3;
        x2   = x1;
        x_bar = x1 + 1.0;
    } else {
        let xx1_raw = -1.5 * tempb / (tempa * tempa);
        let xx1 = if xx1_raw > EPS2 { (1.0 + 2.0 * xx1_raw).sqrt() - 1.0 } else { xx1_raw };
        let (mut xa, mut xb) = (xx1 * tempa3, (-xx1 - 2.0) * tempa3);
        if tempa3 < 0.0 { std::mem::swap(&mut xa, &mut xb); }
        x1 = xa; x2 = xb;
        x_bar = if fcube(x1) > 0.0 { x2 - 1.0 } else { x1 + 1.0 };
    }

    // Regula falsi to bracket root
    let f1   = fcube(x1);
    let l_up = f1 > 0.0;
    let mut xlo = 0.0_f64;
    let mut xhi = 0.0_f64;
    let mut l_lo = false;
    let mut l_hi = false;

    for _ in 0..20 {
        let fbar = fcube(x_bar);
        if fbar <= 0.0 {
            xlo   = x_bar;
            l_lo  = true;
            if l_hi { break; }
            x_bar = if !l_up { x1 + (xlo - x1) * 10.0 }
                    else      { x2 - (x2 - xlo) * 0.1  };
        } else {
            xhi   = x_bar;
            l_hi  = true;
            if l_lo { break; }
            x_bar = if !l_up { x1 + (xhi - x1) * 0.1  }
                    else      { x2 - (x2 - xhi) * 10.0 };
        }
    }

    if !l_lo || !l_hi { return; } // failed to bracket

    // Newton-Raphson 2nd-order refinement within [xlo, xhi]
    let mut x0 = xhi;
    for _ in 0..40 {
        let f0 = fcube(x0);
        if f0.abs() < EPSM * fcuba(x0) { break; }
        let f1v  = fcube1(x0);
        let f0f1 = f0 / f1v;
        let delx = -f0f1 * (1.0 + 0.5 * f0f1 * fcube2(x0) / f1v);
        if delx.abs() < EPSM * x0.abs() { break; }
        x0 = xhi.min(xlo.max(x0 + delx));
    }

    if x0 < 1.0e-30 { return; }

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
    x[n[19] - 1] *= rclx * rnoy;    // ClONO2
    for ni in [n[11], n[12], n[13], n[23]] {
        x[ni - 1] *= rbrx;
    }
    x[n[22] - 1] *= rbrx * rnoy;    // BrONO2

    // Iodine family closure (Iy linked to NOy via IONO2, analogous to BrONO2/ClONO2).
    if s.liod && n[30] > 0 && n[31] > 0 && n[32] > 0 && n[33] > 0 && n[34] > 0 {
        let xxi = x[n[30] - 1] + x[n[31] - 1] + x[n[32] - 1] + x[n[34] - 1];
        let xxin = x[n[33] - 1];
        let denom = xxi + rnoy * xxin;
        if denom > 1.0e-30 {
            let riyx = totiyx / denom;
            for ni in [n[30], n[31], n[32], n[34]] {
                x[ni - 1] *= riyx;
            }
            x[n[33] - 1] *= riyx * rnoy; // IONO2
        }
    }
}

/// Wrapper: rplace → fixrat → splace → update FO3.
/// Fortran: SUBROUTINE FIXMIX
pub fn fixmix(s: &mut ModelState) {
    let ib   = s.ibox;
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

    // Time-step halving loop
    let deltt0 = s.deltt;
    let mut kut = 0;
    for k in 0..30 {
        s.deltt *= 2.0;
        let c = newrax(s, damp1, x, n);
        if c == 0 {
            kut = k + 1;
            break;
        }
        if k == 29 {
            eprintln!(" =====FAILED TO CONVERGE AFTER 2**-30 CUT");
            s.lprts = true;
            // Fortran: diagnostic NEWRAX may still converge; if so, use result
            let c2 = newrax(s, damp1, x, n);
            if c2 == 0 {
                s.deltt = deltt0;
                return Ok(());
            }
            bail!("Newton-Raphson failed to converge after 2^-30 time step cuts");
        }
    }

    // Rebuild back to original time step
    for _ in 0..kut {
        s.deltt /= 2.0;
        let c = newrax(s, damp1, x, n);
        if c != 0 {
            s.lprts = true;
            // Fortran: diagnostic NEWRAX may still converge; if so, use result
            let c2 = newrax(s, damp1, x, n);
            if c2 == 0 {
                s.deltt = deltt0;
                return Ok(());
            }
            bail!("Newton-Raphson diverged during time-step rebuild");
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

    // Load XR from X
    for j in 0..s.ntotx {
        s.xr[j] = x[j];
    }

    for iter in 0..NUMITR {
        s.raxloop += 1.0;
        chems(s);

        if errxo < s.raferr { break; }
        if iter == NUMITR - 1 {
            return 2;
        }

        // RHS
        let fxo = rhslhs(s, 0);

        // Convergence check on P-L residual
        let ntot = s.ntot;
        let errpl = (0..n.min(ntot)).fold(0.0_f64, |acc, j| {
            let denom = s.rp[j] + s.rl[j];
            if denom > 0.0 { acc.max((fxo[j] / denom).abs()) } else { acc }
        });

        if errpl < s.rafpml { break; }

        // Jacobian + solve
        let mut xoo = [0.0f64; NDEN];
        if s.lsvjac {
            resolv_in_place(s, &fxo.as_slice()[..n], &mut xoo[..n]);
        } else {
            rhslhs(s, 1);
            linslv(s, &fxo.as_slice()[..n], &mut xoo[..n], n);
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
                temp = temp
                    .max(s.rafmin * s.xr[j])
                    .min(s.rafmax * s.xr[j]);
            }
            s.xr[j] = temp;
        }

        // Check for negative densities
        if (0..n.min(ntot)).any(|j| s.xr[j] <= 0.0) {
            return 1;
        }
    }

    // Success: write XR back to X
    for j in 0..s.ntotx {
        x[j] = s.xr[j];
    }
    0
}

// ── LINSLV — Crout LU decomposition + solve ───────────────────────────────────

/// Partial-pivoting LU decomposition of s.a_mat, then forward/back solve.
/// Modifies s.a_mat in-place and stores pivots in s.ipa.
/// Fortran: SUBROUTINE LINSLV(B, X, N)
pub fn linslv(s: &mut ModelState, b: &[f64], x: &mut [f64], n: usize) {
    // Crout decomposition with partial pivoting on the A matrix
    // (Fortran stores A column-major; we use row-major Array2 but the algorithm
    //  applies to the same logical matrix)
    let mut s_col = [0.0f64; NDEN];

    for kr in 0..n {
        // Copy column KR of A into S
        for k in 0..n {
            s_col[k] = s.a_mat[[k, kr]];
        }

        // Apply previous eliminations
        if kr > 0 {
            for j in 0..kr {
                let jp = s.ipa[j];
                s.a_mat[[j, kr]] = s_col[jp];
                s_col[jp] = s_col[j];
                let ajkr = s.a_mat[[j, kr]];
                for i in j + 1..n {
                    s_col[i] -= s.a_mat[[i, j]] * ajkr;
                }
            }
        }

        // Find pivot
        let mut smax   = s_col[kr].abs();
        let mut krmax  = kr;
        for i in kr..n {
            if s_col[i].abs() >= smax {
                krmax = i;
                smax  = s_col[i].abs();
            }
        }
        s.ipa[kr] = krmax;
        s.a_mat[[kr, kr]] = s_col[krmax];
        let div = 1.0 / s_col[krmax];
        s_col[krmax] = s_col[kr];

        if kr < n - 1 {
            for i in kr + 1..n {
                s.a_mat[[i, kr]] = s_col[i] * div;
            }
        }
    }

    resolv_in_place(s, b, x);
}

// ── RESOLV — back-solve with stored LU ───────────────────────────────────────

/// Solve using the already-decomposed A matrix and pivot array.
/// Fortran: SUBROUTINE RESOLV(B, X, N)
pub fn resolv_in_place(s: &ModelState, b: &[f64], x: &mut [f64]) {
    let n = b.len();
    let mut sv = [0.0f64; NDEN];
    for i in 0..n { sv[i] = b[i]; }

    // Forward substitution with pivoting
    for i in 0..n {
        let ip = s.ipa[i];
        x[i] = sv[ip];
        sv[ip] = sv[i];
        if i < n - 1 {
            for j in i + 1..n {
                sv[j] -= s.a_mat[[j, i]] * x[i];
            }
        }
    }

    // Back substitution
    for i in (0..n).rev() {
        let mut sum = x[i];
        for j in i + 1..n {
            sum -= s.a_mat[[i, j]] * x[j];
        }
        x[i] = sum / s.a_mat[[i, i]];
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
