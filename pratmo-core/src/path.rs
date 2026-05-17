// bpath.f → trajectory path integration module
// TPATH, DSTEP, NEWATM

use anyhow::Result;

use crate::{
    chemistry::chems,
    jvalue::sol,
    output,
    reader::{hystat, setday},
    solver::{rplace, splace, fixmix},
    state::ModelState,
    constants::NJH2O,
};

// ── TPATH ────────────────────────────────────────────────────────────────────

/// Top-level path integration driver.
/// Sets NBOXDO positive, then calls NEWATM + DSTEP in a loop.
/// Fortran: SUBROUTINE TPATH
/// Note: NEWATM reading from fort01 requires the caller to supply the remaining
/// path records via the `records` iterator. Pass an empty slice for no-path-reset.
pub fn tpath(s: &mut ModelState) -> Result<()> {
    // Initial PRTPTH(0,1,IB) for each box + set NBOXDO positive
    let nbox = s.nbox;
    for ib in 0..nbox {
        output::prtpth(s, 0, 1, ib + 1); // 1-based box index
        s.nboxdo[ib] = s.nboxdo[ib].abs();
    }

    for ipath in 1..=999usize {
        // Read next path record via NEWATM
        s.lresol = true;
        if !newatm(s)? {
            break;
        }

        let nbox = s.nbox;
        for ib in 0..nbox {
            dstep(s, ipath as i32, ib)?;

            if s.lprtx {
                println!("\n----Seg/Day/Box: {:3}{:3}{:3}-----24-hr L=AVG(P-L)",
                    ipath, s.ndaysd, ib + 1);
                // Copy PMEAN into RCOLUM (simplified: just call prtall directly)
                let ntotx = s.ntotx;
                let nprtrr = s.nprtrr;
                let ntot = s.ntot;
                output::prtall(s, 2, 0, ntotx);
                output::prtall(s, 0, nprtrr, ntot);
                let nfval = s.nfval as usize;
                output::prtall(s, 11, 0, nfval);
            }
        }
    }
    Ok(())
}

// ── DSTEP ────────────────────────────────────────────────────────────────────

/// Integrate box `ib` for NDAYSD days.
/// Assumes XXNOFT and J-values are already initialised.
/// Fortran: SUBROUTINE DSTEP(IPATH, IB)
pub fn dstep(s: &mut ModelState, ipath: i32, ib: usize) -> Result<()> {
    let _ipath = ipath;
    const DAMP1: f64 = 0.5;

    s.ibox = ib;
    let ialt = s.nboxdo[ib].unsigned_abs() as usize; // IALT = NBOXDO(IB) (positive)
    if ialt == 0 || ialt > s.nc {
        return Ok(());
    }
    s.ialt = ialt - 1; // 0-based
    let densty = s.dm[s.ialt];
    s.izalt = 0;

    // Load initial XNOFT from saved XXNOFT for this box
    let ntotx = s.ntotx;
    let ntimdo = s.ntimdo;
    for j in 0..ntotx {
        for i in 0..ntimdo {
            s.xnoft[[j, i]] = s.xxnoft[[j, i, ib]];
        }
    }

    // Compute J-values if needed (LRESOL)
    if s.lresol {
        s.lresol = false;
        let njval = s.njval;
        let nbox = s.nbox;

        if s.nday > 1 {
            // Average J-values: daytime → STORJV[:,0,:], nighttime → STORJV[:,1,:]
            s.gmu = 1.0;
            s.ljzer = false;
            sol(s, 1.0);
            for jb in 0..nbox {
                for jj in 0..njval {
                    let v = s.jval_get(jb, jj);
                    s.storjv[[jj, 0, jb]] = v;
                }
            }
            s.gmu = -1.0;
            s.ljzer = true;
            sol(s, -1.0);
            for jb in 0..nbox {
                for jj in 0..njval {
                    let v = s.jval_get(jb, jj);
                    s.storjv[[jj, 1, jb]] = v;
                }
            }
        } else {
            // Full diurnal J-values at NMU solar angles
            let nmu = s.nmu;
            for jn in 0..nmu {
                let gmu = s.utime[jn];
                s.gmu = gmu;
                s.ljzer = gmu < s.gmu0;
                sol(s, gmu);
                for jb in 0..nbox {
                    for jj in 0..njval {
                        let v = s.jval_get(jb, jj);
                        s.storjv[[jj, jn, jb]] = v;
                    }
                }
            }
        }
    }

    // ── Day loop ──────────────────────────────────────────────────────────
    let ndaysd = s.ndaysd;
    let njval = s.njval;

    for iday_idx in 0..ndaysd as usize {
        let _iday = iday_idx;
        // Clear daily averages
        let nndxpq = s.ndxpp.len(); // NNDXPQ
        for k in 0..50 + nndxpq {
            s.ppmean[[k, ib]] = 0.0;
        }
        for i in 0..490usize {
            s.pmean[i] = 0.0;
        }

        // Load starting values from DDDDDD arrays for this box
        let mut xnoft_col0 = [0.0f64; 30];
        rplace(s, &mut xnoft_col0, ib);
        for j in 0..ntotx {
            s.xnoft[[j, 0]] = xnoft_col0[j];
        }

        // Propagate explicitly-integrated species constant across time steps
        for j in s.ntot..ntotx {
            let v0 = s.xnoft[[j, 0]];
            for i in 1..ntimdo {
                s.xnoft[[j, i]] = v0;
            }
        }

        s.xpo3 = 0.0;
        s.xlo3 = 0.0;

        // ── Time-step loop ────────────────────────────────────────────────
        for it in 1..ntimdo {
            s.ittt = it as i32 + 1;
            s.gmu = s.utime[it];
            s.ljzer = s.gmu < s.gmu0;
            let jn = (s.jtim[it] - 1).max(0) as usize; // 0-based

            // Load J-values for current box at this time
            for jj in 0..njval {
                let v = s.storjv[[jj, jn, ib]];
                s.jval_set(ib, jj, v);
            }

            let tdelt = s.dtime[it] - s.dtime[it - 1];
            s.deltt = 1.0 / tdelt;

            // Load guess from stored diurnal cycle
            let mut xn = [0.0f64; 30];
            for j in 0..ntotx {
                xn[j] = s.xnoft[[j, it]];
                s.xnold[j] = s.xnoft[[j, it - 1]];
            }

            // H2O photolysis rate at current altitude and time
            s.xjdo = hunt_xjh2o(s, s.z[s.ialt] * 1e-5).max(0.0)
                * s.utime[it].max(0.0);



            // Newton-Raphson integration
            let result = crate::solver::newraf(s, DAMP1, &mut xn, s.ntot);
            if result.is_err() {
                s.lprts = true;
                crate::solver::newraf(s, DAMP1, &mut xn, s.ntot)?;
            }

            // Update XR and call chemistry for rates
            for j in 0..ntotx {
                s.xr[j] = xn[j];
            }
            chems(s);

            // Accumulate PMEAN
            let weight = tdelt / s.daysec;
            for j in 0..430usize {
                let rcol = rcolum_get(s, j + 1);
                s.pmean[j] += weight * rcol;
            }
            for j in 0..30usize {
                s.pmean[430 + j] += weight * (s.rpf[j] - s.rlf[j]);
            }
            // O3 production/loss tracking
            s.xpo3 += weight * s.r[135]; // R(136) 0-based → r[135]
            s.xlo3 += weight * s.rlf[0];
            for j in 0..30usize {
                s.pmean[460 + j] += weight * (s.rp[j] - s.rl[j]);
            }

            // Store solution
            for j in 0..ntotx {
                s.xnoft[[j, it]] = xn[j];
            }
        } // end time-step loop

        // ── End of day: save noon values and update long-lived species ─────

        // Store noon values back into DDDDDD arrays
        let mut xn_final = [0.0f64; 30];
        for j in 0..ntotx {
            xn_final[j] = s.xnoft[[j, ntimdo - 1]];
        }
        splace(s, &xn_final, ib);

        // Store PPMEAN for this box
        for k in 0..30usize {
            s.ppmean[[k, ib]] = s.pmean[460 + k];
        }
        for k in 0..20usize {
            s.ppmean[[30 + k, ib]] = s.pmean[430 + k];
        }
        for k in 0..nndxpq {
            let ndx = s.ndxpp[k];
            if ndx > 0 {
                s.ppmean[[50 + k, ib]] = s.pmean[30 + ndx - 1];
            }
        }

        // Explicit integration of species flagged by NDIFF
        // ntsav[j] is 1-based slot index; convert to 0-based for xnoft
        for j in 0..ntotx {
            let ntj = s.ntsav[j]; // NT(J): 1-based slot
            if ntj <= s.ntot || s.ndiff[j] == 0 {
                continue;
            }
            let xadd = s.ppmean[[ntj, ib]] * s.daysec;
            let cur = s.den_get(ib, j);
            let new_val = cur + xadd;
            s.den_set(ib, j, new_val);
            if ntj >= 1 { s.xnoft[[ntj - 1, 0]] = new_val; }
        }

        // Explicit integration of long-lived mixing ratios (NFLO > 0)
        for j in 0..s.nfval as usize {
            if s.nflo[j] <= 0 {
                continue;
            }
            let xadd = s.ppmean[[30 + j, ib]] * s.daysec / densty;
            // FFFFFF(IB,J) += XADD → update long-lived mixing ratio j for box ib
            let old_val = s.fff_get(ib, j);
            s.fff_set(ib, j, old_val + xadd);
        }

        // Renormalize O3: DO3(IB) = FO3(IB)*DENSTY or FO3(IB) = DO3(IB)/DENSTY
        // s.n[10] is 1-based slot for O3 (N11 in Fortran), convert to 0-based for xnoft
        let n11 = s.n[10]; // 1-based slot value
        if n11 > s.ntot {
            s.do3[ib] = s.fo3[ib] * densty;
            if n11 >= 1 { s.xnoft[[n11 - 1, 0]] = s.do3[ib]; }
        } else {
            s.fo3[ib] = s.do3[ib] / densty;
        }

        // Renormalize family mixing ratios from current box densities
        // s.n[] contains 1-based XN slot indices; convert to 0-based with -1
        {
            let mut xn_tmp = [0.0f64; 30];
            rplace(s, &mut xn_tmp, ib);

            // Helper: get xn_tmp by 1-based slot index
            let xn = |ni: usize| -> f64 { if ni >= 1 && ni <= 30 { xn_tmp[ni - 1] } else { 0.0 } };
            let n = s.n;
            if s.nflo[2] <= 0 {
                // FNOY = (XN(N1)+XN(N2)+XN(N3)+XN(N4)+XN(N4)+XN(N5)+XN(N15)+XN(N21)+XN(N20)+XN(N23))
                let fnoy = (xn(n[0]) + xn(n[1]) + xn(n[2])
                    + xn(n[3]) + xn(n[3])
                    + xn(n[4]) + xn(n[14]) + xn(n[20])
                    + xn(n[19]) + xn(n[22]))
                    / densty;
                s.fnoy[ib] = fnoy;
            }
            if s.nflo[5] <= 0 {
                let fclx = (xn(n[15]) + xn(n[16]) + xn(n[18])
                    + xn(n[21])
                    + xn(n[27]) + 2.0 * (xn(n[17]) + xn(n[28]))
                    + xn(n[19]) + xn(n[29]))
                    / densty;
                s.fclx[ib] = fclx;
            }
            if s.nflo[15] <= 0 {
                let fbrx = (xn(n[11]) + xn(n[12]) + xn(n[13])
                    + xn(n[23])
                    + xn(n[22]) + xn(n[29]))
                    / densty;
                s.fbrx[ib] = fbrx;
            }
        }

        // Store final diurnal cycle in XXNOFT
        for j in 0..ntotx {
            for i in 0..ntimdo {
                s.xxnoft[[j, i, ib]] = s.xnoft[[j, i]];
            }
        }

        // LPRTX output at end of each day
        if s.lprtx {
            let iday_val = _iday as i32 + 1;
            output::prtpth(s, _ipath, iday_val, ib + 1);
            output::punch(s, ib + 1, iday_val);
        }
    } // end day loop

    Ok(())
}

// ── NEWATM ───────────────────────────────────────────────────────────────────

/// Reset atmosphere, ozone profile, and timing for the next path segment.
/// Reads one record from s.fort01_remaining (populated by READIN).
/// Returns false when NDAYSD==0 (end of path records or EOF).
/// Fortran: SUBROUTINE NEWATM(LAT, MON)
pub fn newatm(s: &mut ModelState) -> Result<bool> {
    // Pop next record from fort01_remaining (reversed so pop gives original order)
    let record = match s.fort01_remaining.pop() {
        Some(line) => line,
        None => {
            s.ndaysd = 0;
            return Ok(false);
        }
    };

    newatm_from_record_inner(s, &record)
}

fn newatm_from_record_inner(s: &mut ModelState, record: &str) -> Result<bool> {
    // Parse: 2X A8 2X A8 I5 3E10.3 I5 2F5.3
    if record.len() < 2 {
        s.ndaysd = 0;
        return Ok(false);
    }
    let titatm = record.get(2..10).unwrap_or("").trim().to_string();
    let tito3  = record.get(12..20).unwrap_or("").trim().to_string();

    let ndaysd: i32 = record.get(20..25).unwrap_or("    0")
        .trim().parse().unwrap_or(0);
    s.ndaysd = ndaysd;
    if ndaysd == 0 {
        return Ok(false);
    }

    let xlatd: f64 = parse_e10(record, 25);
    let xdecd: f64 = parse_e10(record, 35);
    let xo3col: f64 = parse_e10(record, 45);
    let kalt: i32 = record.get(55..60).unwrap_or("    0")
        .trim().parse().unwrap_or(0);
    let xarea: f64 = record.get(60..65).unwrap_or("    0")
        .trim().parse().unwrap_or(0.0);
    let xrain: f64 = record.get(65..70).unwrap_or("    0")
        .trim().parse().unwrap_or(0.0);

    const CRAD: f64 = 57.295_779_513;
    s.xlatd = xlatd;
    s.xdecd = xdecd;
    s.xlat = xlatd / CRAD;
    s.xdec = xdecd / CRAD;

    let mut newat: i32 = 0;

    // Reset atmosphere (fort13) if label is non-blank
    if !titatm.is_empty() {
        newat += 1;
        let fort13 = s.fort13_lines.clone();
        load_atmosphere(s, &titatm, &fort13)?;
    }

    // Reset O3 profile (fort14) if label is non-blank
    if !tito3.is_empty() {
        newat += 1;
        let fort14 = s.fort14_lines.clone();
        load_ozone_profile(s, &tito3, &fort14)?;
    }

    // Rescale densities if atmosphere changed or altitude shift
    if newat > 0 || kalt != 0 {
        // Save old densities before altitude shift
        let mut dmold = [0.0f64; crate::constants::NB];
        for ib in 0..s.nbox {
            let ialt_old = s.nboxdo[ib].unsigned_abs() as usize;
            dmold[ib] = if ialt_old > 0 && ialt_old <= s.nc { s.dm[ialt_old - 1] } else { 1.0 };
        }
        for ib in 0..s.nbox {
            s.nboxdo[ib] += kalt;
            let ialt = (s.nboxdo[ib].unsigned_abs() as usize).min(s.nc);
            let dmfact = if dmold[ib] > 0.0 {
                s.dm[ialt.saturating_sub(1)] / dmold[ib]
            } else { 1.0 };
            for j in 0..s.ndval as usize {
                let old = s.den_get(ib, j);
                s.den_set(ib, j, old * dmfact);
            }
        }
        fixmix(s);
    }

    // Rescale O3 column if XO3COL > 0
    if xo3col > 0.0 {
        let nc = s.nc;
        let zzht = s.zzht;
        let mut xo3c = s.do3ref[nc - 1] * (zzht + 0.5 * (s.z[nc - 1] - s.z[nc - 2]));
        for ii in 1..nc - 1 {
            xo3c += s.do3ref[ii] * 0.5 * (s.z[ii + 1] - s.z[ii - 1]);
        }
        xo3c += s.do3ref[0] * 0.5 * (s.z[1] - s.z[0]);
        let fo3c = xo3col / xo3c;
        for ii in 0..nc {
            s.do3ref[ii] *= fo3c;
        }
    }

    // Rescale aerosol area and rainout efficiency
    if xarea > 0.0 {
        for ib in 0..s.nbox {
            s.boxaa[ib] *= xarea;
        }
    }
    if xrain > 0.0 {
        for ib in 0..s.nbox {
            s.boxrn[ib] *= xrain;
        }
    }

    setday(s);
    s.lresol = true;
    s.izalt = 0;

    Ok(true)
}

// ── Helpers ───────────────────────────────────────────────────────────────────

/// Parse a 10-character E10.3 field at byte offset `off` from a record line.
fn parse_e10(line: &str, off: usize) -> f64 {
    line.get(off..off + 10)
        .unwrap_or("          ")
        .trim()
        .replace(|c: char| c == 'd' || c == 'D', "e")
        .parse::<f64>()
        .unwrap_or(0.0)
}

/// Load a labelled T-profile from fort13 lines into state.
fn load_atmosphere(s: &mut ModelState, label: &str, lines: &[String]) -> Result<()> {
    let nc = s.nc;
    let mut i = 0;
    while i + 1 < lines.len() {
        let hdr = &lines[i];
        let tita = hdr.get(2..10).unwrap_or("").trim();
        if tita == label {
            // Read T(1..46) from next line(s) in 8E10.3 format
            let tline = &lines[i + 1];
            let mut t = Vec::new();
            let mut off = 0;
            while off + 10 <= tline.len() && t.len() < 46 {
                let v: f64 = tline[off..off + 10].trim()
                    .parse().unwrap_or(200.0);
                t.push(v);
                off += 10;
            }
            for j in 0..nc.min(t.len()) {
                s.t[j] = t[j];
            }
            // Z grid: 2 km spacing
            for j in 0..nc {
                s.z[j] = 2.0e5 * j as f64;
            }
            let logp = if label.starts_with('P') { 1 } else { 0 };
            hystat(s, logp);
            return Ok(());
        }
        i += 2;
    }
    anyhow::bail!("Atmosphere label '{}' not found in fort13", label);
}

/// Load a labelled O3 profile from fort14 lines into state.
fn load_ozone_profile(s: &mut ModelState, label: &str, lines: &[String]) -> Result<()> {
    let nc = s.nc;
    let mut i = 0;
    // Skip header line
    i += 1;
    while i + 1 < lines.len() {
        let hdr = &lines[i];
        let titz = hdr.get(2..10).unwrap_or("").trim();
        if titz == label {
            let dline = &lines[i + 1];
            let mut do3 = Vec::new();
            let mut off = 0;
            while off + 10 <= dline.len() && do3.len() < nc {
                let v: f64 = dline[off..off + 10].trim()
                    .parse().unwrap_or(0.0);
                do3.push(v);
                off += 10;
            }
            let ndddz = do3.len();
            for j in 0..ndddz.min(nc) {
                s.do3ref[j] = do3[j];
            }
            if ndddz > 1 {
                let scalez = do3[ndddz - 1] / do3[ndddz - 2];
                for j in ndddz..nc {
                    s.do3ref[j] = s.do3ref[j - 1] * scalez;
                }
            }
            return Ok(());
        }
        i += 2;
    }
    anyhow::bail!("O3 label '{}' not found in fort14", label);
}

/// Linear interpolation into the H2O photolysis rate table.
fn hunt_xjh2o(s: &ModelState, zkm: f64) -> f64 {
    if NJH2O == 0 {
        return 0.0;
    }
    let mut j0 = 0usize;
    for j in 0..NJH2O - 1 {
        if s.zjh2o[j] <= zkm && zkm < s.zjh2o[j + 1] {
            j0 = j;
            break;
        }
    }
    s.xjh2o[j0]
}

/// RCOLUM(j) accessor — 1-based j, mirrors EQUIVALENCE in CCRTS.
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
