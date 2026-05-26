// bdiel.f → diurnal cycle integration module
// DIURN, DAILY, RAFDAY

use anyhow::{bail, Result};
use ndarray::Array2;

use crate::{
    chemistry::{chems, setupr},
    jvalue::sol,
    output,
    solver::{rplace, splace, fixrat, fixmix, linslv},
    state::ModelState,
    constants::{NPMEAN, NSLOWM, NNDXPQ, NDEN},
};

// ── Helpers ──────────────────────────────────────────────────────────────────

/// RCOLUM(j) accessor — 1-based j, maps into CCRTS arrays.
/// RCOLUM(1..30)=XR, (31..280)=R, (281..310)=RP, (311..340)=RL,
///                   (341..370)=RPF, (371..400)=RLF, (401..430)=RQF
fn rcolum_get(s: &ModelState, j: usize) -> f64 {
    match j {
        1..=30  => s.xr[j - 1],
        31..=280  => s.r[j - 31],
        281..=310 => s.rp[j - 281],
        311..=340 => s.rl[j - 311],
        341..=370 => s.rpf[j - 341],
        371..=400 => s.rlf[j - 371],
        401..=430 => s.rqf[j - 401],
        _ => 0.0,
    }
}

// ── DIURN ────────────────────────────────────────────────────────────────────

/// Top-level 24-hour diurnal cycle driver.
/// Iterates over all boxes, calls RAFDAY or DAILY, stores results.
/// Fortran: SUBROUTINE DIURN
pub fn diurn(s: &mut ModelState) -> Result<()> {
    // Compute total box weight for global averaging
    let ibsum: i32 = (0..s.nbox).map(|ib| s.nboxwt[ib]).sum();
    let mut qmean = [0.0f64; NPMEAN];

    s.lresol = true;

    // Write unit-7 header before box loop (before PUNCH(0,0))
    output::diurn_unit7_header(s);

    for ib in 0..s.nbox {
        s.ibox = ib;
        let ialt_abs = s.nboxdo[ib].unsigned_abs() as usize; // IABS(NBOXDO(IB))
        s.izalt = 0;
        if ialt_abs == 0 || ialt_abs > s.nc {
            continue;
        }
        s.ialt = ialt_abs - 1; // 0-based

        s.lsvjac = false;
        s.lsvday = false;

        if s.nboxdo[ib] > 0 {
            rafday(s, ib)?;
        } else {
            let mut xnold = s.xnold;
            rplace(s, &mut xnold, ib);
            s.xnold = xnold;
            daily(s, 1)?;
        }

        // PUNCH(0, 0): write altitude/name header on first box (IB=1 in Fortran)
        if ib == 0 {
            output::punch(s, 0, 0);
        }
        // PUNCH(IB+1, 1): write time series for this box
        output::punch(s, ib + 1, 1);

        // Store XXNOFT from XNOFT
        for kt in 0..s.ntimdo {
            for kn in 0..s.ntotx {
                let val = s.xnoft[[kn, kt]];
                s.xxnoft[[kn, kt, ib]] = val;
            }
        }

        // PPMEAN(K, IB) = PMEAN(460+K),  K=1..30   (P-L for implicit species)
        for k in 0..30usize {
            s.ppmean[[k, ib]] = s.pmean[460 + k];
        }
        // PPMEAN(K+30, IB) = PMEAN(430+K), K=1..20  (P-L for families)
        for k in 0..20usize {
            s.ppmean[[30 + k, ib]] = s.pmean[430 + k];
        }
        // PPMEAN(K+50, IB) from NDXPP-indexed rates
        for k in 0..NNDXPQ {
            let ndx = s.ndxpp[k];
            if ndx > 0 {
                s.ppmean[[50 + k, ib]] = s.pmean[30 + ndx - 1];
            }
        }

        // Accumulate global QMEAN over boxes
        if s.nboxwt[0] != 0 {
            let boxwt = s.nboxwt[ib] as f64 / ibsum as f64;
            for kk in 0..NPMEAN {
                qmean[kk] += s.pmean[kk] * boxwt;
            }
        }
        // LPRTX per-box printout (NBOXPR > 1)
        if s.lprtx && s.nboxpr[ib] > 1 {
            s.ittt = 1;
            println!("\n ----AVERAGE(OVER 24 HRS)-----L=AVG(P-L)-----Box:{:5}", ib + 1);
            let ntotx = s.ntotx;
            let nboxpr_val = s.nboxpr[ib];
            for ii in 0..430 { s.xr[ii.min(29)] = rcolum_get(s, ii + 1); } // simplified
            output::prtall(s, 2, nboxpr_val - 2, ntotx);
            let nfval = s.nfval as usize;
            output::prtall(s, 11, 0, nfval);
        }
    }

    // Global average (LPRT && NBOXWT(1) != 0) — skipped if NBOXWT[0]==0

    // PRTRAT if NBOXWT(1) != 0
    if s.nboxwt[0] != 0 {
        output::prtrat(s, 1);
    }

    s.lsvjac = false;
    s.lsvday = false;
    Ok(())
}

// ── DAILY ────────────────────────────────────────────────────────────────────

/// Single-day 24-hour time-dependent integration.
/// id < 100: save XNOFT; id >= 100: skip XNOFT update (partial-deriv mode).
/// Fortran: SUBROUTINE DAILY(ID)
pub fn daily(s: &mut ModelState, id: i32) -> Result<()> {
    const DAMP1: f64 = 0.5;
    let n = s.ntot;

    // Zero PMEAN
    for j in 0..NPMEAN {
        s.pmean[j] = 0.0;
    }

    // Compute / store J-values if needed
    if s.lresol {
        s.lresol = false;

        if s.nday > 1 {
            // NDAY=2: 4-step (2 day, 2 night); NDAY=3: 1-step 24h avg
            // Daytime average J-values  → STORJV[:,0,:]
            s.ljzer = false;
            sol(s, s.gmu);
            let nbox = s.nbox;
            let njval = s.njval;
            for jb in 0..nbox {
                for jj in 0..njval {
                    let v = s.jval_get(jb, jj);
                    s.storjv[[jj, 0, jb]] = v;
                }
            }
            // Nighttime (zero) J-values → STORJV[:,1,:]
            s.ljzer = true;
            sol(s, s.gmu);
            for jb in 0..nbox {
                for jj in 0..njval {
                    let v = s.jval_get(jb, jj);
                    s.storjv[[jj, 1, jb]] = v;
                }
            }
        } else {
            // Full diurnal: compute J-values at each NMU solar angle
            let nmu = s.nmu;
            for jn in 0..nmu {
                let gmu = s.utime[jn];
                s.gmu = gmu;
                s.ljzer = gmu < s.gmu0;
                sol(s, gmu);
                let nbox = s.nbox;
                let njval = s.njval;
                for jb in 0..nbox {
                    for jj in 0..njval {
                        let v = s.jval_get(jb, jj);
                        s.storjv[[jj, jn, jb]] = v;
                    }
                }
            }
        }
    }

    // Load starting guess from XNOLD
    let mut xn = [0.0f64; NDEN];
    for j in 0..s.ntotx {
        xn[j] = s.xnold[j];
    }

    // If first call since new altitude, initialize XNOFT from XN
    if !s.lsvday {
        let ntimdo = s.ntimdo;
        let ntotx = s.ntotx;
        for it in 0..ntimdo {
            for j in 0..ntotx {
                s.xnoft[[j, it]] = xn[j];
            }
        }
    }

    // ── Time-step loop ────────────────────────────────────────────────────
    let ntimdo = s.ntimdo;
    let njval = s.njval;
    let nbox = s.nbox;

    for it in 0..ntimdo {
        s.ittt = it as i32 + 1; // 1-based like Fortran ITTT
        s.gmu = s.utime[it];
        s.ljzer = s.gmu < s.gmu0;
        let jn = (s.jtim[it] - 1).max(0) as usize; // 0-based

        // Load J-values for this time step into VVVVVV for current box
        for jj in 0..njval {
            let v = s.storjv[[jj, jn, s.ibox]];
            s.jval_set(s.ibox, jj, v);
        }

        if it == 0 {
            // IT==1 in Fortran: skip to label 30 (just set XR = XN, XNOLD = XN, call CHEMS)
            for j in 0..s.ntotx {
                s.xr[j] = xn[j];
                s.xnold[j] = xn[j];
                if id < 100 {
                    s.xnoft[[j, it]] = xn[j];
                }
            }
            chems(s);
            continue;
        }

        // DELTT = 1/dt
        s.deltt = 1.0 / (s.dtime[it] - s.dtime[it - 1]);

        // Load previous-day solution as first guess
        for j in 0..s.ntotx {
            xn[j] = s.xnoft[[j, it]];
        }

        // Newton-Raphson step (NEWRAF handles time-step halving internally)
        let result = crate::solver::newraf(s, DAMP1, &mut xn, n);
        if result.is_err() {
            s.lprts = true;
            // Fortran loops back (GOTO 22) — retry once more
            crate::solver::newraf(s, DAMP1, &mut xn, n)?;
        }

        // Store AEXTRA and update XR, XNOLD
        for j in 0..n {
            s.aextra[j] = xn[j];
        }
        for j in 0..s.ntotx {
            s.xr[j] = xn[j];
            s.xnold[j] = xn[j];
            if id < 100 {
                s.xnoft[[j, it]] = xn[j];
            }
        }

        chems(s);

        // LPRTX printout section — stub (PRTALL)

        // Accumulate PMEAN
        let daysec = s.daysec;
        let weight = (s.dtime[it] - s.dtime[it - 1]) / daysec;
        for j in 0..430usize {
            s.pmean[j] += weight * rcolum_get(s, j + 1);
        }
        for j in 0..30usize {
            s.pmean[430 + j] += weight * (s.rpf[j] - s.rlf[j]);
        }
        for j in 0..NDEN {
            s.pmean[460 + j] += weight * (s.rp[j] - s.rl[j]);
        }
    }

    s.lsvday = true;
    Ok(())
}

// ── RAFDAY ───────────────────────────────────────────────────────────────────

/// Newton-Raphson steady-state driver for NNRT slow species.
/// Runs DAILY to compute 24h means, then applies NR correction.
/// Fortran: SUBROUTINE RAFDAY(IB)
pub fn rafday(s: &mut ModelState, _ib: usize) -> Result<()> {
    let nnr = s.nnr;

    // Local arrays (NSLOWM = 11 max)
    let mut fxdder = Array2::<f64>::zeros((NSLOWM, NSLOWM));
    let mut fxo = [0.0f64; NSLOWM];
    let mut xo  = [0.0f64; NSLOWM];
    let mut xoo = [0.0f64; NDEN];

    let mut lcnvrg = false;
    let lpsave = s.lprtx;
    s.lprtx = nnr < 1 && s.lprtx;

    let maxraf = s.maxraf;
    let maxrlx = s.maxrlx;

    'outer: for itrraf in 0..=maxraf {
        s.lsvjac = false;
        s.lprtx = nnr < 1 && s.lprtx;

        // Relaxation phase
        if maxrlx >= 1 {
            let mut xnold = s.xnold;
            rplace(s, &mut xnold, s.ibox);
            s.xnold = xnold;

            for _ in 0..maxrlx {
                daily(s, 2)?;
            }

            let xnold_snap = s.xnold;
            splace(s, &xnold_snap, s.ibox);
            fixmix(s);
        }

        {
            let mut xnold = s.xnold;
            rplace(s, &mut xnold, s.ibox);
            s.xnold = xnold;
        }

        s.lprtx = lpsave && (itrraf >= maxraf || lcnvrg);
        daily(s, 3)?;

        if nnr < 1 {
            break 'outer;
        }

        if lcnvrg || itrraf >= maxraf {
            break 'outer;
        }

        // Store 24h-mean P-L in FXO (RHS of NR system)
        for j in 0..nnr {
            let jn = s.nnrt[j]; // 1-based species id in NNRT
            let ntjn1 = s.n[jn - 1]; // 1-based NR slot index (Fortran NT)
            if ntjn1 == 0 {
                continue;
            }
            // Fortran PMEAN(460+NTJN) => 0-based index (459 + ntjn1)
            fxo[j] = s.pmean[459 + ntjn1];
        }

        // Build finite-difference Jacobian
        let save_lprtx = s.lprtx;
        s.lprtx = false;

        for j in 0..nnr {
            let jn = s.nnrt[j];
            let ntjn1 = s.n[jn - 1];
            if ntjn1 == 0 {
                continue;
            }
            let ntjn0 = ntjn1 - 1; // 0-based slot for XNOLD

            let mut xnold = s.xnold;
            rplace(s, &mut xnold, s.ibox);
            s.xnold = xnold;

            let epslon = s.dayeps * s.xnold[ntjn0];
            s.xnold[ntjn0] += epslon;

            {
                let xnold_snap = s.xnold;
                fixrat(&mut { let mut tmp = xnold_snap; tmp }, s, s.ibox);
            }
            daily(s, 102)?;

            for jj in 0..nnr {
                let jjn = s.nnrt[jj];
                let ntjjn1 = s.n[jjn - 1];
                if ntjjn1 == 0 {
                    continue;
                }
                fxdder[[jj, j]] = (s.pmean[459 + ntjjn1] - fxo[jj]) / epslon;
            }
        }

        s.lprtx = save_lprtx;

        // Copy FXDDER into A matrix for LINSLV
        for j in 0..nnr {
            for jj in 0..nnr {
                s.a_mat[[jj, j]] = fxdder[[jj, j]];
            }
        }
        let mut xo_vec = vec![0.0f64; nnr];
        linslv(s, &fxo[..nnr], &mut xo_vec, nnr);
        xo[..nnr].copy_from_slice(&xo_vec);

        // Apply correction to slow species with clamping
        {
            let mut xnold = s.xnold;
            rplace(s, &mut xnold, s.ibox);
            s.xnold = xnold;
        }

        let mut lneg = false;
        let mut relerr = 0.0f64;
        for j in 0..nnr {
            let jn = s.nnrt[j];
            let ntjn1 = s.n[jn - 1];
            if ntjn1 == 0 {
                continue;
            }
            let ntjn0 = ntjn1 - 1;
            xoo[j] = s.xnold[ntjn0];
            let mut temp = xoo[j] - xo[j];
            if itrraf < 7 {
                let lo = s.rafmin * xoo[j];
                let hi = s.rafmax * xoo[j];
                temp = temp.max(lo).min(hi);
            }
            if temp <= 0.0 {
                lneg = true;
            }
            relerr = relerr.max(xo[j].abs() / xoo[j].max(1e-100));
            xoo[j] = temp;
            s.xnold[ntjn0] = temp;
        }

        lcnvrg = relerr < s.dayerr;

        if lneg {
            bail!("RAFDAY: negative density after correction");
        }

        {
            let mut xnold_snap = s.xnold;
            fixrat(&mut xnold_snap, s, s.ibox);
            splace(s, &xnold_snap, s.ibox);
        }
    }

    s.lprtx = lpsave;

    if nnr > 0 && !lcnvrg {
        eprintln!(" ***RAFDAY NON-CONV AT BOX/ALT={}/{}", s.ibox + 1, s.ialt + 1);
    }

    Ok(())
}
