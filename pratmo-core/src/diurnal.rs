// bdiel.f → diurnal cycle integration module
// DIURN, DAILY, RAFDAY

use anyhow::{bail, Result};
use ndarray::Array2;
use rayon::prelude::*;
use std::sync::OnceLock;

use crate::{
    chemistry::chems,
    constants::{NDEN, NNDXPQ, NPMEAN, NSLOWM, PPMEAN_FAMILY_OFFSET, PPMEAN_RATE_OFFSET},
    jvalue::sol,
    output,
    solver::{fixmix, fixrat, linslv, rplace, splace},
    state::ModelState,
};

// ── Helpers ──────────────────────────────────────────────────────────────────

/// RCOLUM(j) accessor — 1-based j, maps into CCRTS arrays.
/// RCOLUM(1..30)=XR, (31..280)=R, (281..310)=RP, (311..340)=RL,
///                   (341..370)=RPF, (371..400)=RLF, (401..430)=RQF
fn rcolum_get(s: &ModelState, j: usize) -> f64 {
    match j {
        1..=30 => s.xr[j - 1],
        31..=280 => s.r[j - 31],
        281..=310 => s.rp[j - 281],
        311..=340 => s.rl[j - 311],
        341..=370 => s.rpf[j - 341],
        371..=400 => s.rlf[j - 371],
        401..=430 => s.rqf[j - 401],
        _ => 0.0,
    }
}

fn rafday_reuse_jacobian_enabled() -> bool {
    static ENABLED: OnceLock<bool> = OnceLock::new();
    *ENABLED.get_or_init(|| {
        std::env::var("PRATMO_RAFDAY_REUSE_JACOBIAN")
            .map(|v| {
                matches!(
                    v.as_str(),
                    "1" | "true" | "TRUE" | "yes" | "YES" | "on" | "ON"
                )
            })
            .unwrap_or(false)
    })
}

// ── DIURN ────────────────────────────────────────────────────────────────────

/// Top-level 24-hour diurnal cycle driver.
/// Iterates over all boxes, calls RAFDAY or DAILY, stores results.
/// Fortran: SUBROUTINE DIURN
pub fn diurn(s: &mut ModelState) -> Result<()> {
    // The original DIURN executable leaves COMMON/CHRIS/SSF uninitialized
    // (gfortran's static storage makes it zero), which suppresses all
    // photolysis in the legacy reference run.  Keep the physically useful
    // Rust default SSF=1 outside parity mode, but reproduce that observable
    // behavior when byte/numerical parity is explicitly requested.
    #[cfg(feature = "fortran-parity")]
    s.ssf.fill(0.0);

    // Compute total box weight for global averaging
    let ibsum: i32 = (0..s.nbox).map(|ib| s.nboxwt[ib]).sum();
    let mut qmean = [0.0f64; NPMEAN];

    s.lresol = true;

    // Write unit-7 header before box loop (before PUNCH(0,0))
    output::diurn_unit7_header(s);

    // Fortran calls PUNCH(0,0) before entering the box loop.  At this point
    // IALT still contains the last altitude selected while loading fort02.x
    // (the final BOXDO entry), which is the legacy altitude written in the
    // unit-7 metadata record.  Parity mode keeps this pre-loop ordering;
    // normal Rust mode writes the metadata after selecting the first valid
    // box, which is the less surprising API behavior but has a different
    // legacy altitude record.
    #[cfg(feature = "fortran-parity")]
    output::punch(s, 0, 0);

    for ib in 0..s.nbox {
        s.ibox = ib;
        let ialt_abs = s.nboxdo[ib].unsigned_abs() as usize; // IABS(NBOXDO(IB))
        s.izalt = 0;
        if ialt_abs == 0 || ialt_abs > s.nlev.max(s.nc) {
            continue;
        }
        s.ialt = ialt_abs - 1; // 0-based

        s.lsvjac = false;
        s.lsvday = false;

        if s.nboxdo[ib] > 0 {
            solve_diurnal_box(s, ib)?;
        } else {
            let mut xnold = s.xnold;
            rplace(s, &mut xnold, ib);
            s.xnold = xnold;
            daily(s, 1)?;
        }

        #[cfg(not(feature = "fortran-parity"))]
        if ib == 0 {
            // Normal Rust mode emits PUNCH metadata here, after the first box
            // selected its altitude.  Parity mode emitted the Fortran
            // pre-loop record above instead.
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

        // P-L for all implicit species.
        for k in 0..NDEN {
            s.ppmean[[k, ib]] = s.pmean[460 + k];
        }
        // P-L for explicit families.
        for k in 0..20usize {
            s.ppmean[[PPMEAN_FAMILY_OFFSET + k, ib]] = s.pmean[430 + k];
        }
        // NDXPP-indexed diagnostic rates.
        for k in 0..NNDXPQ {
            let ndx = s.ndxpp[k];
            if ndx > 0 {
                s.ppmean[[PPMEAN_RATE_OFFSET + k, ib]] = s.pmean[30 + ndx - 1];
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
            println!(
                "\n ----AVERAGE(OVER 24 HRS)-----L=AVG(P-L)-----Box:{:5}",
                ib + 1
            );
            let ntotx = s.ntotx;
            let nboxpr_val = s.nboxpr[ib];
            for ii in 0..430 {
                s.xr[ii.min(29)] = rcolum_get(s, ii + 1);
            } // simplified
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

struct BoxDiurnResult {
    ib: usize,
    state: Box<ModelState>,
}

/// Parallel DIURN variant for structured API runs.
///
/// Each active box is solved in an isolated clone of the prepared model state,
/// then per-box outputs are merged back in deterministic box order. File-backed
/// Fortran-style output is intentionally not supported here; callers with output
/// units attached fall back to the sequential DIURN path.
pub fn diurn_parallel_boxes(s: &mut ModelState) -> Result<()> {
    if s.out_unit7.is_some() || s.out_unit8.is_some() || s.out_unit9.is_some() || s.nbox <= 1 {
        return diurn(s);
    }

    let jobs: Vec<usize> = (0..s.nbox)
        .filter(|&ib| {
            let ialt_abs = s.nboxdo[ib].unsigned_abs() as usize;
            ialt_abs != 0 && ialt_abs <= s.nlev.max(s.nc)
        })
        .collect();

    let base = Box::new(s.clone());
    let results: Result<Vec<BoxDiurnResult>> = jobs
        .into_par_iter()
        .map(|ib| {
            let mut worker = base.clone();
            worker.rafday_warnings.clear();
            worker.out_unit7 = None;
            worker.out_unit8 = None;
            worker.out_unit9 = None;
            worker.ibox = ib;
            worker.izalt = 0;
            worker.ialt = worker.nboxdo[ib].unsigned_abs() as usize - 1;
            worker.lresol = true;
            worker.lsvjac = false;
            worker.lsvday = false;

            if worker.nboxdo[ib] > 0 {
                solve_diurnal_box(&mut worker, ib)?;
            } else {
                let mut xnold = worker.xnold;
                rplace(&worker, &mut xnold, ib);
                worker.xnold = xnold;
                daily(&mut worker, 1)?;
            }

            Ok(BoxDiurnResult { ib, state: worker })
        })
        .collect();

    let mut results = results?;
    results.sort_by_key(|r| r.ib);

    s.raxloop = 0.0;
    s.radcount = 0.0;
    s.newraf_nonconvergence_count = 0;
    s.rafday_nonconvergence_count = 0;
    s.rafday_warnings.clear();
    s.rafday_max_final_relative_correction = 0.0;
    s.rafday_max_correction_iterations = 0;
    for result in results {
        let ib = result.ib;
        let worker = result.state;

        for spec in 0..NDEN {
            let val = worker.den_get(ib, spec);
            s.den_set(ib, spec, val);
        }
        for family in 1..=19 {
            let val = worker.fff_get(ib, family);
            s.fff_set(ib, family, val);
        }
        for jj in 0..s.njval {
            let val = worker.jval_get(ib, jj);
            s.jval_set(ib, jj, val);
            for kt in 0..s.storjv.shape()[1] {
                s.storjv[[jj, kt, ib]] = worker.storjv[[jj, kt, ib]];
            }
        }
        for kt in 0..s.ntimdo {
            for kn in 0..s.ntotx {
                s.xxnoft[[kn, kt, ib]] = worker.xnoft[[kn, kt]];
            }
        }
        for k in 0..NDEN {
            s.ppmean[[k, ib]] = worker.pmean[460 + k];
        }
        for k in 0..20usize {
            s.ppmean[[PPMEAN_FAMILY_OFFSET + k, ib]] = worker.pmean[430 + k];
        }
        for k in 0..NNDXPQ {
            let ndx = s.ndxpp[k];
            if ndx > 0 {
                s.ppmean[[PPMEAN_RATE_OFFSET + k, ib]] = worker.pmean[30 + ndx - 1];
            }
        }

        s.raxloop += worker.raxloop;
        s.radcount += worker.radcount;
        s.newraf_nonconvergence_count += worker.newraf_nonconvergence_count;
        s.rafday_nonconvergence_count += worker.rafday_nonconvergence_count;
        s.rafday_warnings.extend(worker.rafday_warnings);
        s.rafday_max_final_relative_correction = s
            .rafday_max_final_relative_correction
            .max(worker.rafday_max_final_relative_correction);
        s.rafday_max_correction_iterations = s
            .rafday_max_correction_iterations
            .max(worker.rafday_max_correction_iterations);
    }

    s.lsvjac = false;
    s.lsvday = false;
    Ok(())
}

// ── DAILY ────────────────────────────────────────────────────────────────────

fn solve_diurnal_box(s: &mut ModelState, ib: usize) -> Result<()> {
    #[cfg(not(feature = "fortran-parity"))]
    if s.evolve_ozone && !s.cpp_compatibility {
        return evolving_ozone_days(s, ib);
    }
    if s.cpp_compatibility {
        cpp_endpoint_days(s, ib)
    } else {
        rafday(s, ib)
    }
}

/// Close the full diurnal orbit after solving the slow-species daily balance.
/// Updating ozone can substantially change the fast radicals' noon state;
/// RAFDAY alone tests only its slow coordinates and leaves those fast initial
/// values behind. Refresh them from the completed orbit and solve again.
#[cfg(not(feature = "fortran-parity"))]
fn evolving_ozone_days(s: &mut ModelState, ib: usize) -> Result<()> {
    let saved_relaxation = s.maxrlx;
    let result = (|| {
        let max_cycles = s.nboxmx[ib].max(1) as usize;
        let mut residual = f64::INFINITY;
        for cycle in 0..max_cycles {
            let failures = s.rafday_nonconvergence_count;
            let step_failures = s.newraf_nonconvergence_count;
            rafday(s, ib)?;
            let endpoint = s.xnold;
            splace(s, &endpoint, ib);
            s.fo3[ib] = s.do3[ib] / s.dm[s.ialt];
            if s.rafday_nonconvergence_count != failures {
                return Ok(());
            }
            if s.newraf_nonconvergence_count != step_failures {
                s.rafday_nonconvergence_count += 1;
                rafday_warn(s, "timestep failures during ozone relaxation");
                return Ok(());
            }
            residual = (0..s.ntot).fold(0.0_f64, |worst, slot| {
                let start = s.xnoft[[slot, 0]];
                let end = s.xnoft[[slot, s.ntimdo - 1]];
                worst.max((end - start).abs() / start.abs().max(1.0e-2))
            });
            if residual < s.dayerr {
                return Ok(());
            }
            if cycle + 1 < max_cycles {
                fixmix(s);
                s.lsvday = false;
                s.maxrlx = 0;
            }
        }
        s.rafday_nonconvergence_count += 1;
        rafday_warn(s, &format!(
            "noon-to-noon relative change {residual:.3e} exceeds {:.3e} after {max_cycles} convergence cycles",
            s.dayerr
        ));
        Ok(())
    })();
    s.maxrlx = saved_relaxation;
    result
}

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

        if s.cpp_compatibility {
            // The C++ reaction graph conserves coupled NOy/Cly/Bry exactly.
            // The legacy array solver represents those coupled families with
            // FIXRAT; apply it at each accepted step to avoid accumulating a
            // spurious family drift over the C++ multi-day relaxation.
            fixrat(&mut xn, s, s.ibox);
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

/// Later C++ box-model convergence driver.
///
/// Unlike Fortran RAFDAY, the C++ implementation simply advances complete
/// diurnal cycles until every integrated radical changes by less than 0.5%
/// between the two noon endpoints. ``nboxmx`` supplies the maximum number of
/// days, matching ``SetNumDaysForConvergence``.
fn cpp_endpoint_days(s: &mut ModelState, ib: usize) -> Result<()> {
    const TOLERANCE: f64 = 5.0e-3;

    let mut initial = s.xnold;
    rplace(s, &mut initial, ib);
    s.xnold = initial;
    s.lsvday = false;

    let max_days = s.nboxmx[ib].max(1) as usize;
    let mut converged = false;
    for day in 0..max_days {
        let step_failures = s.newraf_nonconvergence_count;
        daily(s, 1)?;
        if s.evolve_ozone && s.newraf_nonconvergence_count != step_failures {
            bail!(
                "evolving ozone: timestep integration failed (box {})",
                ib + 1
            );
        }
        let mut max_ratio = 0.0_f64;
        for slot in 0..s.ntot {
            let start = s.xnoft[[slot, 0]];
            let end = s.xnoft[[slot, s.ntimdo - 1]];
            let ratio = if s.evolve_ozone {
                if !start.is_finite() || !end.is_finite() || start < 0.0 || end < 0.0 {
                    bail!("evolving ozone: invalid daily orbit (box {})", ib + 1);
                }
                (end - start).abs() / start.abs().max(1.0e-2)
            } else if start.abs() < 1.0e-2 {
                0.05 * TOLERANCE
            } else {
                ((end - start) / start).abs()
            };
            max_ratio = max_ratio.max(ratio);
        }
        if max_ratio < TOLERANCE {
            converged = true;
            break;
        }
        if s.evolve_ozone && day + 1 < max_days {
            let endpoint = s.xnold;
            splace(s, &endpoint, ib);
            fixmix(s);
            let mut initial = s.xnold;
            rplace(s, &mut initial, ib);
            s.xnold = initial;
            s.lsvday = false;
        }
    }

    let final_density = s.xnold;
    splace(s, &final_density, ib);
    if s.evolve_ozone {
        s.fo3[ib] = s.do3[ib] / s.dm[s.ialt];
    }
    if !converged {
        s.rafday_nonconvergence_count += 1;
        if s.evolve_ozone {
            s.rafday_warnings.push(format!(
                "DIURN did not converge for box {} at {:.3} km: noon-to-noon endpoint tolerance {TOLERANCE:.3e} not reached after {max_days} days. Returning the last valid diurnal cycle; photochemical equilibrium is not established.",
                ib + 1, s.z[s.ialt] * 1.0e-5
            ));
        }
    }
    Ok(())
}

// ── RAFDAY ───────────────────────────────────────────────────────────────────

#[cfg(not(feature = "fortran-parity"))]
fn rafday_family_valid(s: &ModelState, x: &[f64; NDEN]) -> bool {
    if x[..s.ntotx].iter().any(|v| !v.is_finite() || *v < 0.0) {
        return false;
    }
    let v = |species: usize| x[s.n[species] - 1];
    let noy = v(0) + v(1) + v(2) + 2.0 * v(3) + v(4) + v(14) + v(20) + v(19) + v(22);
    let cly = v(15) + v(16) + 2.0 * v(17) + v(18) + v(19) + v(21) + v(27) + 2.0 * v(28) + v(29);
    let bry = v(11) + v(12) + v(13) + v(22) + v(23) + v(29);
    let mut families = vec![
        (noy, s.fnoy[s.ibox]),
        (cly, s.fclx[s.ibox]),
        (bry, s.fbrx[s.ibox]),
    ];
    if s.liod {
        let iy = (30..36).map(v).sum::<f64>() + 2.0 * (36..40).map(v).sum::<f64>();
        families.push((iy, s.fiodx[s.ibox]));
    }
    families.into_iter().all(|(total, mixing)| {
        let target = mixing * s.dm[s.ialt];
        (total - target).abs() <= 1.0e-8 * target.abs().max(1.0e-30)
    })
}

/// Evaluate the same constrained daily chemistry used after an accepted step.
/// SETUPR reads named HCl/HBr/ClONO2 densities for aerosol uptake, so changing
/// only XNOLD omits that dependence from the finite-difference Jacobian.
/// Cloning also keeps rejected trials out of the orbit cache and diagnostics.
#[cfg(not(feature = "fortran-parity"))]
fn rafday_trial(s: &ModelState, candidate: &[f64; NDEN]) -> Result<Box<ModelState>> {
    if !rafday_family_valid(s, candidate) {
        bail!(
            "RAFDAY: invalid density or family constraint in trial (box {})",
            s.ibox + 1
        );
    }
    let mut trial = Box::new(s.clone());
    splace(&mut trial, candidate, s.ibox);
    trial.xnold = *candidate;
    trial.lprtx = false;
    trial.lsvday = false;
    daily(&mut trial, 102)?;
    if trial.newraf_nonconvergence_count != s.newraf_nonconvergence_count
        || trial.xnold[..s.ntotx]
            .iter()
            .any(|v| !v.is_finite() || *v < 0.0)
        || trial.pmean.iter().any(|v| !v.is_finite())
    {
        bail!(
            "RAFDAY: daily integration failed in trial (box {})",
            s.ibox + 1
        );
    }
    Ok(trial)
}

/// Safeguard the whole Newton direction, including FIXRAT's coupled changes.
/// The merit function is sum((daysec * mean(P-L) / initial_density)^2), with
/// scales held fixed during backtracking. Small damped steps cannot establish
/// convergence: both the *undamped* correction and daily residual must pass.
#[cfg(not(feature = "fortran-parity"))]
fn rafday_correction(s: &mut ModelState, fxo: &[f64], correction: &[f64]) -> Result<(f64, bool)> {
    let mut initial = [0.0; NDEN];
    rplace(s, &mut initial, s.ibox);
    let slots: Vec<_> = s.nnrt[..s.nnr]
        .iter()
        .map(|&species| s.n[species - 1] - 1)
        .collect();
    let scales: Vec<_> = slots.iter().map(|&slot| initial[slot] / s.daysec).collect();
    let norm = |residual: &[f64]| {
        residual
            .iter()
            .zip(&scales)
            .map(|(f, scale)| (f / scale).powi(2))
            .sum::<f64>()
    };
    let before = norm(fxo);
    let mut alpha = 1.0_f64;
    let mut relerr = 0.0_f64;
    let mut residual_error = 0.0_f64;
    for (j, &slot) in slots.iter().enumerate() {
        let dx = correction[j];
        if !initial[slot].is_finite()
            || initial[slot] <= 0.0
            || !dx.is_finite()
            || !scales[j].is_finite()
            || scales[j] <= 0.0
            || !fxo[j].is_finite()
        {
            bail!(
                "RAFDAY: invalid Newton state for {} (box {})",
                s.tnamet[slot],
                s.ibox + 1
            );
        }
        relerr = relerr.max((dx / initial[slot]).abs());
        residual_error = residual_error.max((fxo[j] / scales[j]).abs());
        if dx > 0.0 {
            // Stay strictly inside the positive domain, at every iteration.
            alpha = alpha.min(0.9 * initial[slot] / dx);
        }
    }
    if !before.is_finite() || !rafday_family_valid(s, &initial) {
        bail!(
            "RAFDAY: invalid residual or initial family constraint (box {})",
            s.ibox + 1
        );
    }
    // Do not demand strict decrease below numerical resolution once both
    // independent convergence criteria have already been met.
    if relerr < s.dayerr && residual_error < s.dayerr {
        return Ok((relerr, true));
    }
    for _ in 0..24 {
        let mut candidate = initial;
        for (&slot, &dx) in slots.iter().zip(correction) {
            candidate[slot] -= alpha * dx;
        }
        fixrat(&mut candidate, s, s.ibox);
        if let Ok(trial) = rafday_trial(s, &candidate) {
            let residual: Vec<_> = slots.iter().map(|&slot| trial.pmean[460 + slot]).collect();
            let after = norm(&residual);
            if after.is_finite() && after <= (1.0 - 1.0e-4 * alpha) * before {
                splace(s, &candidate, s.ibox);
                return Ok((relerr, false));
            }
        }
        alpha *= 0.5;
    }
    bail!(
        "RAFDAY: no positive, family-conserving residual-reducing Newton step (box {})",
        s.ibox + 1
    )
}

/// Validate the small slow-species solve without changing LINSLV's legacy
/// behavior in the time-step integrator or Fortran parity builds.
#[cfg(not(feature = "fortran-parity"))]
fn rafday_check_solve(
    s: &ModelState,
    jacobian: &Array2<f64>,
    rhs: &[f64],
    step: &[f64],
) -> Result<()> {
    let n = rhs.len();
    for i in 0..n {
        let pivot = s.a_mat.as_slice().expect("a_mat is contiguous")[i * NDEN + i];
        if !pivot.is_finite() || pivot == 0.0 || !step[i].is_finite() || !rhs[i].is_finite() {
            bail!(
                "RAFDAY: singular or nonfinite Newton system (box {})",
                s.ibox + 1
            );
        }
        let mut residual = -rhs[i];
        let mut scale = rhs[i].abs();
        for j in 0..n {
            let term = jacobian[[i, j]] * step[j];
            residual += term;
            scale += term.abs();
        }
        if !residual.is_finite() || !scale.is_finite() || residual.abs() > 1.0e-10 * scale {
            bail!("RAFDAY: inaccurate Newton solve (box {})", s.ibox + 1);
        }
    }
    Ok(())
}

/// Newton-Raphson steady-state driver for NNRT slow species.
/// Runs DAILY to compute 24h means, then applies NR correction.
/// Fortran: SUBROUTINE RAFDAY(IB)
pub fn rafday(s: &mut ModelState, _ib: usize) -> Result<()> {
    #[cfg(feature = "fortran-parity")]
    {
        rafday_inner(s)
    }
    #[cfg(not(feature = "fortran-parity"))]
    {
        let mut checkpoint = None;
        let lprtx = s.lprtx;
        let result = rafday_inner(s, &mut checkpoint);
        s.lprtx = lprtx;
        if let Err(error) = result {
            let Some(mut valid) = checkpoint else {
                // Without a successfully integrated orbit there is no usable
                // output to return. Input/initial-integration errors stay fatal.
                return Err(error);
            };
            // ModelState::clone intentionally omits file handles.
            valid.out_unit7 = s.out_unit7.take();
            valid.out_unit8 = s.out_unit8.take();
            valid.out_unit9 = s.out_unit9.take();
            std::mem::swap(s, valid.as_mut());
            s.lprtx = lprtx;
            s.rafday_nonconvergence_count += 1;
            rafday_warn(s, &error.to_string());
        }
        Ok(())
    }
}

#[cfg(not(feature = "fortran-parity"))]
fn rafday_warn(s: &mut ModelState, reason: &str) {
    s.rafday_warnings.push(format!(
        "RAFDAY did not converge for box {} at {:.3} km: {}. Returning the last valid diurnal cycle; photochemical equilibrium is not established.",
        s.ibox + 1, s.z[s.ialt] * 1.0e-5, reason
    ));
}

fn rafday_inner(
    s: &mut ModelState,
    #[cfg(not(feature = "fortran-parity"))] checkpoint: &mut Option<Box<ModelState>>,
) -> Result<()> {
    let nnr = s.nnr;

    // Local arrays (NSLOWM = 11 max)
    let mut fxdder = Array2::<f64>::zeros((NSLOWM, NSLOWM));
    let mut fxo = [0.0f64; NSLOWM];
    let mut xo = [0.0f64; NSLOWM];

    let mut lcnvrg = false;
    let lpsave = s.lprtx;
    s.lprtx = nnr < 1 && s.lprtx;

    let maxraf = s.maxraf;
    let maxrlx = s.maxrlx;
    let reuse_jacobian = rafday_reuse_jacobian_enabled();
    let mut have_jacobian = false;
    let mut final_relative_correction = 0.0f64;
    let mut correction_iterations = 0usize;

    'outer: for itrraf in 0..=maxraf {
        s.lsvjac = false;
        s.lprtx = nnr < 1 && s.lprtx;

        // Relaxation phase
        // Seed RAFDAY with a relaxed diurnal orbit once. Repeating this block
        // after every Newton step integrates the slow species for MAXRLX more
        // days and erases the HNO3 correction that FIXRAT couples into BrONO2.
        // The standard DAILY(3) call below already refreshes the orbit after
        // each correction. Strict parity mode retains the legacy repetition
        // because the original Fortran contains the same known defect.
        let relax_this_iteration = maxrlx >= 1 && (itrraf == 0 || cfg!(feature = "fortran-parity"));
        if relax_this_iteration {
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
        #[cfg(not(feature = "fortran-parity"))]
        let failures_before = s.newraf_nonconvergence_count;
        #[cfg(not(feature = "fortran-parity"))]
        {
            // Match the deterministic starting guesses used for Jacobian and
            // line-search trials; changing cached guesses can change adaptive
            // time-step paths and make the residual discontinuous.
            s.lsvday = false;
        }
        daily(s, 3)?;

        #[cfg(not(feature = "fortran-parity"))]
        {
            let mut initial = [0.0; NDEN];
            rplace(s, &mut initial, s.ibox);
            if s.newraf_nonconvergence_count != failures_before
                || !rafday_family_valid(s, &initial)
                || s.pmean.iter().any(|v| !v.is_finite())
                || (0..s.ntimdo).any(|it| {
                    (0..s.ntotx).any(|slot| {
                        let value = s.xnoft[[slot, it]];
                        !value.is_finite() || value < 0.0
                    })
                })
            {
                bail!("RAFDAY: invalid daily orbit (box {})", s.ibox + 1);
            }
            let mut valid = Box::new(s.clone());
            valid.rafday_max_final_relative_correction = valid
                .rafday_max_final_relative_correction
                .max(final_relative_correction);
            valid.rafday_max_correction_iterations = valid
                .rafday_max_correction_iterations
                .max(correction_iterations);
            *checkpoint = Some(valid);
        }

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

        // Build finite-difference Jacobian. In experimental fast mode, reuse
        // the first slow-species Jacobian for later RAFDAY iterations.
        if !reuse_jacobian || !have_jacobian || itrraf % 2 == 0 {
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

                if !epslon.is_finite() || epslon <= 0.0 {
                    bail!(
                        "RAFDAY: invalid Jacobian perturbation for {}",
                        s.tname[jn - 1]
                    );
                }
                let mut fixed = s.xnold;
                fixrat(&mut fixed, s, s.ibox);
                #[cfg(not(feature = "fortran-parity"))]
                let means = rafday_trial(s, &fixed)?.pmean;
                #[cfg(feature = "fortran-parity")]
                let means = {
                    s.xnold = fixed;
                    daily(s, 102)?;
                    s.pmean
                };

                for jj in 0..nnr {
                    let jjn = s.nnrt[jj];
                    let ntjjn1 = s.n[jjn - 1];
                    if ntjjn1 == 0 {
                        continue;
                    }
                    fxdder[[jj, j]] = (means[459 + ntjjn1] - fxo[jj]) / epslon;
                }
            }

            s.lprtx = save_lprtx;
            have_jacobian = true;
        }

        // With bromine disabled, CHEMPL gives HBr identically zero P-L.
        // Hold its Newton coordinate fixed instead of solving a singular row.
        #[cfg(not(feature = "fortran-parity"))]
        if !s.lbrom {
            for j in 0..nnr {
                if s.nnrt[j] == 14 {
                    if fxo[j] != 0.0 {
                        bail!("RAFDAY: nonzero HBr residual with bromine disabled");
                    }
                    for k in 0..nnr {
                        fxdder[[j, k]] = if j == k { 1.0 } else { 0.0 };
                    }
                }
            }
        }
        // Scale the optional ozone solve: ozone and trace reservoirs can
        // differ by many orders of magnitude. Solve for fractional density
        // corrections, equilibrate rows, and verify the original equations.
        let mut column_scale = [1.0_f64; NSLOWM];
        let mut row_scale = [1.0_f64; NSLOWM];
        if s.evolve_ozone {
            for j in 0..nnr {
                column_scale[j] = s.den_get(s.ibox, s.nnrt[j] - 1).abs().max(1.0e-36);
            }
            for i in 0..nnr {
                row_scale[i] = (0..nnr)
                    .map(|j| (fxdder[[i, j]] * column_scale[j]).abs())
                    .fold(fxo[i].abs().max(1.0e-36), f64::max);
            }
        }
        // Copy FXDDER into A matrix for LINSLV.
        let a = s.a_mat.as_slice_mut().expect("a_mat is contiguous");
        for j in 0..nnr {
            for jj in 0..nnr {
                a[j * NDEN + jj] = fxdder[[jj, j]] * column_scale[j] / row_scale[jj];
            }
        }
        let mut xo_vec = vec![0.0f64; nnr];
        let rhs: Vec<_> = (0..nnr).map(|j| fxo[j] / row_scale[j]).collect();
        linslv(s, &rhs, &mut xo_vec, nnr);
        for j in 0..nnr {
            xo[j] = xo_vec[j] * column_scale[j];
        }
        #[cfg(not(feature = "fortran-parity"))]
        {
            rafday_check_solve(s, &fxdder, &fxo[..nnr], &xo[..nnr])?;
            (final_relative_correction, lcnvrg) = rafday_correction(s, &fxo[..nnr], &xo[..nnr])?;
            correction_iterations += 1;
        }

        // Strict parity retains the original component bounds and stopping rule.
        #[cfg(feature = "fortran-parity")]
        {
            {
                let mut xnold = s.xnold;
                rplace(s, &mut xnold, s.ibox);
                s.xnold = xnold;
            }

            let mut xoo = [0.0f64; NDEN];
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
            final_relative_correction = relerr;
            correction_iterations += 1;

            if lneg {
                bail!("RAFDAY: negative density after correction");
            }

            {
                let mut xnold_snap = s.xnold;
                fixrat(&mut xnold_snap, s, s.ibox);
                splace(s, &xnold_snap, s.ibox);
            }
        }
    }

    s.lprtx = lpsave;

    s.rafday_max_final_relative_correction = s
        .rafday_max_final_relative_correction
        .max(final_relative_correction);
    s.rafday_max_correction_iterations = s
        .rafday_max_correction_iterations
        .max(correction_iterations);

    if nnr > 0 && !lcnvrg {
        s.rafday_nonconvergence_count += 1;
        #[cfg(not(feature = "fortran-parity"))]
        rafday_warn(s, "Newton iteration limit reached");
    }

    Ok(())
}

#[cfg(all(test, not(feature = "fortran-parity")))]
mod rafday_tests {
    use super::*;
    use crate::reader::{FortranReader, ModelReader};

    fn prepared_box() -> Box<ModelState> {
        let mut s = ModelState::new();
        FortranReader::embedded().read_all(&mut s).unwrap();
        s.nbox = 1;
        s.ibox = 0;
        s.ialt = 11;
        s.nboxdo[0] = 12;
        s.boxrn[0] = 0.0;
        s.lbrom = true;
        s.lresol = true;
        s.lsvday = false;
        s.lprtx = false;
        fixmix(&mut s);
        let mut initial = [0.0; NDEN];
        rplace(&s, &mut initial, 0);
        s.xnold = initial;
        daily(&mut s, 2).unwrap();
        let relaxed = s.xnold;
        splace(&mut s, &relaxed, 0);
        fixmix(&mut s);
        s
    }

    #[test]
    fn oversized_newton_direction_is_damped_without_clipping_or_family_drift() {
        let mut s = prepared_box();
        let n = s.nnr;
        let mut initial = [0.0; NDEN];
        rplace(&s, &mut initial, 0);
        let base = rafday_trial(&s, &initial).unwrap();
        let slots: Vec<_> = s.nnrt[..n].iter().map(|&id| s.n[id - 1] - 1).collect();
        let rhs: Vec<_> = slots.iter().map(|&slot| base.pmean[460 + slot]).collect();
        for j in 0..n {
            let mut perturbed = initial;
            let eps = s.dayeps * initial[slots[j]];
            perturbed[slots[j]] += eps;
            fixrat(&mut perturbed, &s, 0);
            let trial = rafday_trial(&s, &perturbed).unwrap();
            for i in 0..n {
                s.a_mat.as_slice_mut().unwrap()[j * NDEN + i] =
                    (trial.pmean[460 + slots[i]] - rhs[i]) / eps;
            }
        }
        let mut step = vec![0.0; n];
        linslv(&mut s, &rhs, &mut step, n);
        // Deliberately overshoot along an actual chemistry Newton direction.
        for value in &mut step {
            *value *= 100.0;
        }
        assert!(slots
            .iter()
            .zip(&step)
            .any(|(&slot, &dx)| dx > initial[slot]));
        let cache = s.xnoft.clone();
        let failures = s.newraf_nonconvergence_count;
        let (relative_correction, converged) = rafday_correction(&mut s, &rhs, &step).unwrap();
        assert!(!converged, "a small damped step must not claim convergence");
        assert!(relative_correction > 1.0);
        let mut accepted = [0.0; NDEN];
        rplace(&s, &mut accepted, 0);
        assert!(rafday_family_valid(&s, &accepted));
        assert!(accepted[..s.ntotx].iter().all(|&v| v > 0.0));
        // H2O2 and CH3OOH are not rescaled by FIXRAT: both must have the
        // same alpha, proving that the coupled direction was not clipped.
        let alpha = |species: usize| {
            let j = s.nnrt[..n].iter().position(|&id| id == species).unwrap();
            (initial[slots[j]] - accepted[slots[j]]) / step[j]
        };
        assert!(alpha(9) > 0.0 && alpha(9) < 1.0);
        assert!((alpha(9) - alpha(27)).abs() < 1.0e-12);
        let trial = rafday_trial(&s, &accepted).unwrap();
        let norm = |means: &[f64]| {
            slots
                .iter()
                .map(|&slot| (means[460 + slot] / initial[slot]).powi(2))
                .sum::<f64>()
        };
        assert!(norm(&trial.pmean) < norm(&base.pmean));
        assert_eq!(s.xnoft, cache, "trial orbit leaked into accepted cache");
        assert_eq!(s.newraf_nonconvergence_count, failures);
    }

    #[test]
    fn rejected_direction_leaves_the_box_and_orbit_unchanged() {
        let mut s = prepared_box();
        let mut initial = [0.0; NDEN];
        rplace(&s, &mut initial, 0);
        let base = rafday_trial(&s, &initial).unwrap();
        let rhs: Vec<_> = s.nnrt[..s.nnr]
            .iter()
            .map(|&id| base.pmean[459 + s.n[id - 1]])
            .collect();
        let cache = s.xnoft.clone();
        let failures = s.newraf_nonconvergence_count;
        // No change cannot decrease a nonzero residual.
        let step = vec![0.0; s.nnr];
        assert!(rafday_correction(&mut s, &rhs, &step).is_err());
        let mut after = [0.0; NDEN];
        rplace(&s, &mut after, 0);
        assert_eq!(initial, after);
        assert_eq!(s.xnoft, cache);
        assert_eq!(s.newraf_nonconvergence_count, failures);
    }

    #[test]
    fn singular_nonfinite_and_inaccurate_slow_solves_are_rejected() {
        let mut s = ModelState::new();
        let jac = Array2::<f64>::eye(2);
        let rhs = [1.0, 2.0];
        // Zero pivot, even with a finite-looking supplied solution.
        assert!(rafday_check_solve(&s, &jac, &rhs, &rhs).is_err());
        s.a_mat.as_slice_mut().unwrap()[0] = 1.0;
        s.a_mat.as_slice_mut().unwrap()[NDEN + 1] = 1.0;
        assert!(rafday_check_solve(&s, &jac, &rhs, &rhs).is_ok());
        assert!(rafday_check_solve(&s, &jac, &rhs, &[f64::NAN, 2.0]).is_err());
        assert!(rafday_check_solve(&s, &jac, &[f64::INFINITY, 2.0], &rhs).is_err());
        assert!(rafday_check_solve(&s, &jac, &rhs, &[1.0, 3.0]).is_err());
    }

    #[test]
    fn recovery_requires_a_valid_orbit() {
        let mut s = prepared_box();
        s.maxrlx = 0;
        s.fnoy[0] = f64::NAN;
        assert!(rafday(&mut s, 0).is_err());
        assert!(s.rafday_warnings.is_empty());
    }
}
