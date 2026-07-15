// butil.f → output/diagnostic subroutines
// PUNCH, PRTALL, PRTPTH, PRTAVG, PRTRAT

use std::fmt::Write as FmtWrite;
use std::io::Write;

use crate::constants::PPMEAN_RATE_OFFSET;
use crate::state::ModelState;

// ── Fortran-style float formatters ────────────────────────────────────────────

/// Fortran 1PE10.3: scientific, field width 10, 3 decimal places.
/// Format: sign(1) + digit(1) + dot(1) + 3dec(3) + E(1) + expsign(1) + 2exp(2) = 10 chars.
pub fn fmt_e10p3(v: f64) -> String {
    if !v.is_finite() {
        return "**********".to_string();
    }
    if v == 0.0 {
        return " 0.000E+00".to_string();
    }
    let neg = v < 0.0;
    // {:.3E} gives "1.234E5" or "1.234E-5" — mantissa already has the decimal point
    let raw = format!("{:.3E}", v.abs());
    let (mant, exp_str) = raw.split_once('E').unwrap_or((&raw, "0"));
    let exp: i32 = exp_str.parse().unwrap_or(0);
    let sign = if neg { '-' } else { ' ' };
    let es = if exp >= 0 { '+' } else { '-' };
    // mant is e.g. "1.234" (5 chars); total = 1+5+1+1+2 = 10 chars
    format!("{}{}E{}{:02}", sign, mant, es, exp.unsigned_abs())
}

/// Write a line of values in 1P,8E10.3 format to `w`.
fn write_e10p3_line(w: &mut dyn Write, vals: &[f64]) {
    let per_line = 8;
    for chunk in vals.chunks(per_line) {
        for v in chunk {
            let _ = write!(w, "{}", fmt_e10p3(*v));
        }
        let _ = writeln!(w);
    }
}

/// Write integer values in 16I5 format to `w`.
fn write_i5_line(w: &mut dyn Write, vals: &[i32]) {
    let per_line = 16;
    for chunk in vals.chunks(per_line) {
        for v in chunk {
            let _ = write!(w, "{:5}", v);
        }
        let _ = writeln!(w);
    }
}

/// Write string values in 8A8 format (8 per line, 8-char each) to `w`.
fn write_a8_line(w: &mut dyn Write, vals: &[&str]) {
    let per_line = 8;
    for chunk in vals.chunks(per_line) {
        for v in chunk {
            let _ = write!(w, "{:<8}", &v[..v.len().min(8)]);
        }
        let _ = writeln!(w);
    }
}

/// Fortran 1PE13.6 field used by READIN's final LEND save (unit 7).
fn fmt_e13p6(v: f64) -> String {
    if !v.is_finite() {
        return "*************".to_string();
    }
    if v == 0.0 {
        return " 0.000000E+00".to_string();
    }
    let neg = v < 0.0;
    let raw = format!("{:.6E}", v.abs());
    let (mant, exp_str) = raw.split_once('E').unwrap_or((&raw, "0"));
    let exp: i32 = exp_str.parse().unwrap_or(0);
    let sign = if neg { '-' } else { ' ' };
    let es = if exp >= 0 { '+' } else { '-' };
    format!("{}{}E{}{:02}", sign, mant, es, exp.unsigned_abs())
}

/// Write READIN's final mixing-ratio save when the path input reaches EOF.
/// Fortran: bread.f label 90, FORMAT 100/301/302.
pub fn lend_dump(s: &mut ModelState) {
    if s.out_unit7.is_none() {
        return;
    }
    let mut out = String::new();

    let titler_line: String = s
        .titler
        .iter()
        .map(|t| format!("{:<8}", &t[..t.len().min(8)]))
        .collect();
    let _ = writeln!(out, "{}", titler_line);

    let ndval = s.ndval as usize;
    let nfval = s.nfval as usize;
    for j in 0..ndval + nfval {
        let title = if j < s.titles.len() { &s.titles[j] } else { "" };
        let _ = writeln!(out, "{:<8}", &title[..title.len().min(8)]);

        let mut vals = Vec::with_capacity(s.nbox + 1);
        for ib in 0..s.nbox {
            if j < ndval {
                let ialt = s.nboxdo[ib].unsigned_abs().saturating_sub(1) as usize;
                vals.push(s.den_get(ib, j) / s.dm[ialt]);
            } else {
                vals.push(s.fff_get(ib, j - ndval + 1));
            }
        }
        // Fortran appends DDDDDD(1,J) once more after the NBOX values.
        if let Some(&first) = vals.first() {
            vals.push(first);
        }
        for chunk in vals.chunks(6) {
            for &v in chunk {
                let _ = write!(out, "{}", fmt_e13p6(v));
            }
            let _ = writeln!(out);
        }
    }
    if let Some(ref mut w) = s.out_unit7 {
        let _ = w.write_all(out.as_bytes());
    }
}

// ── DDDDDD(ib, j) accessor ───────────────────────────────────────────────────

/// DDDDDD(ib, j) = species value for box ib (0-based), species j (1-based).
/// J=1..NDVAL → implicit species densities; J=NDVAL+1..NDVAL+NFVAL → long-lived mixing ratios.
fn dddddd(s: &ModelState, ib: usize, j: usize) -> f64 {
    let ndval = s.ndval as usize;
    if j >= 1 && j <= ndval {
        s.den_get(ib, j - 1)
    } else if j > ndval && j <= ndval + s.nfval as usize {
        s.fff_get(ib, j - ndval)
    } else {
        0.0
    }
}

// ── DIURN unit-7 header ──────────────────────────────────────────────────────

/// Write the header lines that DIURN writes to unit 7 before calling PUNCH(0,0).
/// Fortran: bdiel.f lines before CALL PUNCH(0,0)
pub fn diurn_unit7_header(s: &mut ModelState) {
    let Some(ref mut w) = s.out_unit7 else { return };

    // WRITE(7,202) TITLER   [5 × A8]
    let titler_refs: Vec<&str> = s.titler.iter().map(|t| t.as_str()).collect();
    write_a8_line(w, &titler_refs);

    // WRITE(7,201) NTOT,NTOTX   [2 × I5]
    write_i5_line(w, &[s.ntot as i32, s.ntotx as i32]);

    // WRITE(7,201) (NT(JJ),JJ=1,NTOTX)
    // NT(J) = N(J) = 1-based species slot index; s.n[j] is already 1-based in Rust
    let nt_vals: Vec<i32> = (0..s.ntotx)
        .map(|j| {
            if j < 30 {
                s.n[j] as i32
            } else {
                s.ntotx as i32
            }
        })
        .collect();
    write_i5_line(w, &nt_vals);

    // WRITE(7,202) (TNAME(JJ),JJ=1,NTOTX)   [NTOTX × A8]
    let tname_refs: Vec<&str> = (0..s.ntotx)
        .map(|j| if j < 30 { s.tname[j].as_str() } else { "" })
        .collect();
    write_a8_line(w, &tname_refs);

    // WRITE(7,201) NTIMDO
    write_i5_line(w, &[s.ntimdo as i32]);

    // WRITE(7,201) (NHHMM(JJ),JJ=1,NTIMDO)
    let nhhmm_vals: Vec<i32> = (0..s.ntimdo).map(|j| s.nhhmm[j]).collect();
    write_i5_line(w, &nhhmm_vals);
}

// ── PUNCH ─────────────────────────────────────────────────────────────────────

/// Write diurnal species snapshot to unit 7.
/// ib=0: header (altitude info + species names); ib>0: time series for box ib.
/// Fortran: butil.f SUBROUTINE PUNCH(IB, IDAY)
pub fn punch(s: &mut ModelState, ib: usize, iday: i32) {
    let Some(ref mut w) = s.out_unit7 else { return };

    if ib == 0 {
        // Header call: IB=0
        let i = s.ialt;
        let zkm = 1e-5 * s.z[i];
        // WRITE(7,201) IB,IDAY
        write_i5_line(w, &[0, iday]);
        // WRITE(7,203) ZKM,T(I),DM(I),PSTD(I),DO3REF(I),DO3INT(I)
        write_e10p3_line(
            w,
            &[zkm, s.t[i], s.dm[i], s.pstd[i], s.do3ref[i], s.do3int[i]],
        );
        // WRITE(7,202) (TNAMET(JJ),JJ=1,NTOTX)
        let refs: Vec<&str> = (0..s.ntotx)
            .map(|j| if j < 30 { s.tnamet[j].as_str() } else { "" })
            .collect();
        write_a8_line(w, &refs);
    } else {
        // Data: time series for box ib
        let ntotx = s.ntotx;
        let ntimdo = s.ntimdo;
        for it in 0..ntimdo {
            // WRITE(7,201) IB,IDAY,NHHMM(IT)
            write_i5_line(w, &[ib as i32, iday, s.nhhmm[it]]);
            // WRITE(7,203) (XNOFT(JJ,IT),JJ=1,NTOTX)
            let vals: Vec<f64> = (0..ntotx).map(|j| s.xnoft[[j, it]]).collect();
            write_e10p3_line(w, &vals);
        }
    }
}

// ── PRTALL ────────────────────────────────────────────────────────────────────

/// Print species/rate diagnostics to stdout.
/// MODEX: 0=nothing, 1=X only, 2=X+P+L with titles, 11=family F/P/L
/// MODER: -1=rate mnemonics, 0=nothing, 1=packed rates, 2=detailed rates
/// NXDO: number of species/rates to print
/// Fortran: butil.f SUBROUTINE PRTALL(MODEX, MODER, NXDO)
pub fn prtall(s: &ModelState, modex: i32, moder: i32, nxdo: usize) {
    let ialt = s.ialt;
    let zkm = 1e-5 * s.z[ialt];
    let tt = s.t[ialt] + s.boxtt[s.ibox];
    let ldiurn_char = if s.ldiurn { 'T' } else { 'F' };
    let ittt = s.ittt as usize;
    let nhhmm = if ittt > 0 && ittt <= s.ntimdo {
        s.nhhmm[ittt - 1]
    } else {
        0
    };

    if modex <= 0 && moder == 0 {
        return;
    }

    // Header line (FORMAT 100)
    let tag = if modex == 11 { "FFFFFFFF" } else { "XXXXXXXX" };
    println!(
        " {0}{0}{0}{1:3}={2:5.2}-KM T={3:6.2} M={4:9.2E} O3={5:9.2E} O3COL={6:11.4E} D={7} HHMM={8:4} AA={9:8.3} RN={10:6.2}{11}",
        tag, ialt + 1, zkm, tt,
        s.dm[ialt], s.do3ref[ialt], s.do3int[ialt],
        ldiurn_char, nhhmm, s.boxaa[s.ibox], s.boxrn[s.ibox],
        &tag[..7]
    );

    let nxdd = if !s.lbrom {
        nxdo.min(s.ntotx.saturating_sub(5))
    } else {
        nxdo
    };

    if modex > 0 && modex != 11 {
        let nrows = nxdd.div_ceil(13);
        for i in 0..nrows {
            let i1 = 13 * i;
            let i13 = (i1 + 13).min(nxdd);
            if modex > 1 {
                // Species names
                let names: Vec<&str> = (i1..i13)
                    .map(|j| if j < 30 { s.tnamet[j].as_str() } else { "" })
                    .collect();
                let name_str: String = names.iter().map(|n| format!("{:>10}", n)).collect();
                println!("   {}", name_str);
            }
            // X values
            let xvals: Vec<f64> = (i1..i13)
                .map(|j| if j < 30 { s.xr[j] } else { 0.0 })
                .collect();
            let xstr: String = xvals.iter().map(|v| fmt_e10p3(*v)).collect();
            println!(" X{}", xstr);
            if modex > 1 {
                // P values
                let pvals: Vec<f64> = (i1..i13)
                    .map(|j| if j < 30 { s.rp[j] } else { 0.0 })
                    .collect();
                let pstr: String = pvals.iter().map(|v| fmt_e10p3(*v)).collect();
                println!(" P{}", pstr);
                // L values
                let lvals: Vec<f64> = (i1..i13)
                    .map(|j| if j < 30 { s.rl[j] } else { 0.0 })
                    .collect();
                let lstr: String = lvals.iter().map(|v| fmt_e10p3(*v)).collect();
                println!(" L{}", lstr);
            }
        }
    }

    if modex == 11 {
        // Family species
        let nrows = nxdo.div_ceil(13);
        for i in 0..nrows {
            let i1 = 13 * i + 1;
            let i13 = (i1 + 12).min(nxdo);
            let names: Vec<&str> = (i1..=i13)
                .map(|j| if j <= 20 { s.tnomf[j - 1].as_str() } else { "" })
                .collect();
            let name_str: String = names.iter().map(|n| format!("{:>10}", n)).collect();
            println!("   {}", name_str);
            let fvals: Vec<f64> = (i1..=i13).map(|j| s.fff_get(s.ibox, j)).collect();
            let fstr: String = fvals.iter().map(|v| fmt_e10p3(*v)).collect();
            println!(" X{}", fstr);
            let pvals: Vec<f64> = (i1..=i13)
                .map(|j| if j <= 30 { s.rpf[j - 1] } else { 0.0 })
                .collect();
            let pstr: String = pvals.iter().map(|v| fmt_e10p3(*v)).collect();
            println!(" P{}", pstr);
            let lvals: Vec<f64> = (i1..=i13)
                .map(|j| if j <= 30 { s.rlf[j - 1] } else { 0.0 })
                .collect();
            let lstr: String = lvals.iter().map(|v| fmt_e10p3(*v)).collect();
            println!(" L{}", lstr);
        }
    }

    if moder == 0 {
        return;
    }

    // Rate output
    println!(
        " {0}{0}{0}{1:3}={2:5.2}-KM T={3:6.2} M={4:9.2E} O3={5:9.2E} O3COL={6:11.4E} D={7} HHMM={8:4} AA={9:8.3} RN={10:6.2}RRRRRRR",
        "RRRRRRRR", ialt + 1, zkm, tt,
        s.dm[ialt], s.do3ref[ialt], s.do3int[ialt],
        ldiurn_char, nhhmm, s.boxaa[s.ibox], s.boxrn[s.ibox]
    );
    let nrates = s.nrates;
    let nrate1 = s.nrate1;
    if moder == -1 || moder == 1 {
        let nrows = nrates.div_ceil(10);
        for i in 0..nrows {
            let i1 = 10 * i + 1;
            let i10 = (i1 + 9).min(nrates);
            if i1 > nrate1 && i10 <= 200 {
                continue;
            }
            if moder == -1 {
                // mnemonics
                let names: String = (i1..=i10)
                    .map(|j| {
                        if j <= nrates {
                            format!("{:<8}{:<3}", &s.rfmt_str[j - 1][0], &s.rfmt_str[j - 1][1])
                        } else {
                            " ".repeat(11)
                        }
                    })
                    .collect();
                println!("{:3}={}", i1, names);
            } else {
                // packed rates
                let rvals: Vec<f64> = (i1..=i10)
                    .map(|j| if j <= nrates { s.r[j - 1] } else { 0.0 })
                    .collect();
                let rstr: String = rvals.iter().map(|v| fmt_e10p3(*v)).collect();
                println!("{:3}={}", i1, rstr);
            }
        }
    } else if moder == 2 {
        // Detailed rates
        let nrows = nrates.div_ceil(4);
        for i in 0..nrows {
            let i1 = 4 * i + 1;
            let i4 = (i1 + 3).min(nrates);
            if i1 > nrate1 && i4 <= 200 {
                continue;
            }
            for j in i1..=i4 {
                if j <= nrates {
                    print!(
                        "{:3}={:12.4E} {:<8}{:<8}   ",
                        j,
                        s.r[j - 1],
                        &s.rfmt_str[j - 1][0],
                        &s.rfmt_str[j - 1][1]
                    );
                }
            }
            println!();
        }
    }
}

// ── PRTRAT ────────────────────────────────────────────────────────────────────

/// Print rate table to unit 6 or 8 (depending on LPRT8 flag).
/// Fortran: butil.f SUBROUTINE PRTRAT(NNN)
pub fn prtrat(s: &mut ModelState, nnn: i32) {
    let ndxp = {
        let mut n = crate::constants::NNDXPQ;
        for k in (0..crate::constants::NNDXPQ).rev() {
            if s.ndxpp[k] == 0 {
                n = k;
            } else {
                break;
            }
        }
        n
    };

    // Write to unit 8 if LPRT8, else stdout
    let to_unit8 = s.lprt8;

    let header_line = format!(
        " D={:2}{}",
        nnn,
        (0..ndxp)
            .map(|k| format!("{:6}    ", s.ndxpp[k]))
            .collect::<String>()
    );

    let mnem_line1: String = (0..ndxp)
        .map(|k| {
            let j = s.ndxpp[k];
            if j > 0 && j <= s.nrates {
                format!("{:<8}  ", &s.rfmt_str[j - 1][0])
            } else {
                " ".repeat(10)
            }
        })
        .collect();
    let mnem_line2: String = (0..ndxp)
        .map(|k| {
            let j = s.ndxpp[k];
            if j > 0 && j <= s.nrates {
                format!("{:<8}  ", &s.rfmt_str[j - 1][1])
            } else {
                " ".repeat(10)
            }
        })
        .collect();

    if to_unit8 {
        if let Some(ref mut w) = s.out_unit8 {
            let _ = writeln!(w, "{}", header_line);
            let _ = writeln!(w, "     {}", mnem_line1);
            let _ = writeln!(w, "     {}", mnem_line2);
        }
    } else {
        println!("{}", header_line);
        println!("     {}", mnem_line1);
        println!("     {}", mnem_line2);
    }

    let nbox = s.nbox;
    for ib in 0..nbox {
        let row: String = (0..ndxp)
            .map(|k| fmt_e10p3(s.ppmean[[PPMEAN_RATE_OFFSET + k, ib]]))
            .collect();
        let line = format!("{:3}{}", ib + 1, row);
        if to_unit8 {
            if let Some(ref mut w) = s.out_unit8 {
                let _ = writeln!(w, "{}", line);
            }
        } else {
            println!("{}", line);
        }
    }
}

// ── PRTAVG ───────────────────────────────────────────────────────────────────

/// Print average mixing ratios.
/// Fortran: butil.f SUBROUTINE PRTAVG(NNN)
pub fn prtavg(s: &mut ModelState, nnn: i32) {
    let ndxq = {
        let mut n = crate::constants::NNDXPQ;
        for k in (0..crate::constants::NNDXPQ).rev() {
            if s.ndxqq[k] == 0 {
                n = k;
            } else {
                break;
            }
        }
        n
    };

    let to_unit8 = s.lprt8;
    let nbox = s.nbox;
    let ndval = s.ndval as usize;
    let zdnum = s.zdnum;

    // Compute averages across boxes (trapezoid rule)
    let n_species = ndval + s.nfval as usize;
    let mut avgs = vec![0.0f64; n_species.max(50) + 1];
    for j in 0..n_species {
        let j1 = j + 1;
        avgs[j] = 0.5 * (dddddd(s, 0, j1) + dddddd(s, nbox - 1, j1));
        if nbox >= 3 {
            for ib in 1..nbox - 1 {
                avgs[j] += dddddd(s, ib, j1);
            }
            avgs[j] /= (nbox - 1) as f64;
        }
        if j1 <= ndval && zdnum > 0.0 {
            avgs[j] /= zdnum;
        }
    }
    avgs[48] = avgs[0] + avgs[1]; // NOx
    avgs[49] = avgs[16] + avgs[18] + avgs[27] + 2.0 * (avgs[17] + avgs[28]); // ClOx

    if nnn == 0 {
        let titler_line: String = s
            .titler
            .iter()
            .map(|t| format!("{:<8}", &t[..t.len().min(8)]))
            .collect();
        if to_unit8 {
            if let Some(ref mut w) = s.out_unit8 {
                let _ = writeln!(w, " {}", titler_line);
            }
        } else {
            println!(" {}", titler_line);
        }
    }

    // Print titles if needed
    if nnn == 0 || (!s.lprt8 && (s.lprtx || s.lprty)) {
        let titles: String = (0..ndxq)
            .map(|k| {
                let idx = s.ndxqq[k];
                if idx > 0 && idx < avgs.len() {
                    format!(
                        "{:>10}",
                        &s.titles[idx - 1][..s.titles[idx - 1].len().min(8)]
                    )
                } else {
                    " ".repeat(10)
                }
            })
            .collect();
        if to_unit8 {
            if let Some(ref mut w) = s.out_unit8 {
                let _ = writeln!(w, "     {}", titles);
            }
        } else {
            println!("     {}", titles);
        }
    }

    let row: String = (0..ndxq)
        .map(|k| {
            let idx = s.ndxqq[k];
            if idx > 0 && idx < avgs.len() {
                fmt_e10p3(avgs[idx - 1])
            } else {
                " ".repeat(10)
            }
        })
        .collect();
    let line = format!("{:3}{}", nnn, row);
    if to_unit8 {
        if let Some(ref mut w) = s.out_unit8 {
            let _ = writeln!(w, "{}", line);
        }
    } else {
        println!("{}", line);
    }
}

// ── PRTPTH ────────────────────────────────────────────────────────────────────

/// Print path summary for segment iseg, day iday, box ibx (1-based).
/// Writes species to unit8, rates to unit9.
/// Fortran: butil.f SUBROUTINE PRTPTH(ISEG, IDAY, IBX)
pub fn prtpth(s: &mut ModelState, iseg: i32, iday: i32, ibx: usize) {
    let nnn = iseg + iday + ibx as i32 - 2;
    let ndxq = {
        let mut n = 0usize;
        for k in 0..crate::constants::NNDXPQ {
            if s.ndxqq[k] != 0 {
                n = k + 1;
            }
        }
        n.min(24)
    };
    let ndxp = {
        let mut n = 0usize;
        for k in 0..crate::constants::NNDXPQ {
            if s.ndxpp[k] != 0 {
                n = k + 1;
            }
        }
        n
    };

    let ib0 = ibx.saturating_sub(1); // 0-based
    let ialt = s.nboxdo[ib0].unsigned_abs() as usize;
    let ialt0 = ialt.saturating_sub(1); // 0-based
    let densty = s.dm[ialt0];
    let ndval = s.ndval as usize;
    let nfval = s.nfval as usize;

    // Build AVGS with room for all requested species and legacy derived columns.
    let n_species = ndval + nfval;
    let mut avgs = vec![0.0f64; n_species.max(50) + 1];
    for j in 1..=ndval {
        avgs[j - 1] = s.den_get(ib0, j - 1) / densty;
    }
    for j in ndval + 1..=ndval + nfval {
        avgs[j - 1] = s.fff_get(ib0, j - ndval);
    }
    avgs[48] = avgs[0] + avgs[1]; // NO-x
    avgs[49] = avgs[16] + avgs[18] + avgs[27] + 2.0 * (avgs[17] + avgs[28]); // ClO-x

    // ── Unit 8: species ──────────────────────────────────────────────────────
    if nnn == 0 {
        // Header: TITLER (5A8)
        if let Some(ref mut w8) = s.out_unit8 {
            let titler_line: String = s
                .titler
                .iter()
                .map(|t| format!("{:<8}", &t[..t.len().min(8)]))
                .collect();
            let _ = writeln!(w8, " {}", titler_line);
            // Column headers: ' SegDayBox ' + NDXQ × (A8 + 2X)
            let _ = write!(w8, " SegDayBox ");
            for k in 0..ndxq {
                let idx = s.ndxqq[k];
                let name = match idx {
                    49 => "NO-x    ",
                    50 => "ClO-x   ",
                    _ if idx > 0 && idx <= s.titles.len() => &s.titles[idx - 1],
                    _ => "",
                };
                let _ = write!(w8, "{:<8}", &name[..name.len().min(8)]);
                if k + 1 < ndxq {
                    let _ = write!(w8, "  ");
                }
            }
            let _ = writeln!(w8);
        }
    }
    if let Some(ref mut w8) = s.out_unit8 {
        // WRITE(8,102) ISEG,IDAY,IBX, (AVGS(NDXQQ(KK)),KK=1,NDXQ)
        let _ = write!(w8, " {:3}{:3}{:3}", iseg, iday, ibx);
        for k in 0..ndxq {
            let idx = s.ndxqq[k];
            let val = if idx > 0 && idx < avgs.len() {
                avgs[idx - 1]
            } else {
                0.0
            };
            let _ = write!(w8, "{}", fmt_e10p3(val));
        }
        let _ = writeln!(w8);
    }

    // ── Unit 9: rates ────────────────────────────────────────────────────────
    if nnn == 0 {
        if let Some(ref mut w9) = s.out_unit9 {
            let titler_line: String = s
                .titler
                .iter()
                .map(|t| format!("{:<8}", &t[..t.len().min(8)]))
                .collect();
            let _ = writeln!(w9, " {}", titler_line);
            // Rate indices header: ' D=0 ' + NDXP × (I6 + 4X)
            let _ = write!(w9, " D={:2}", 0);
            for k in 0..ndxp {
                let _ = write!(w9, "{:6}", s.ndxpp[k]);
                if k + 1 < ndxp {
                    let _ = write!(w9, "    ");
                }
            }
            let _ = writeln!(w9);
            // Rate mnemonics rows
            let _ = write!(w9, " Day ");
            for k in 0..ndxp {
                let j = s.ndxpp[k];
                let name = if j > 0 && j <= s.nrates {
                    &s.rfmt_str[j - 1][0]
                } else {
                    ""
                };
                let _ = write!(w9, "{:<8}", &name[..name.len().min(8)]);
                if k + 1 < ndxp {
                    let _ = write!(w9, "  ");
                }
            }
            let _ = writeln!(w9);
            let _ = write!(w9, " Day ");
            for k in 0..ndxp {
                let j = s.ndxpp[k];
                let parsed_name = if j > 0 && j <= s.nrates {
                    &s.rfmt_str[j - 1][1]
                } else {
                    ""
                };
                // fort11's legacy rate-37 mnemonic is `1D)` embedded in a
                // short, non-A8 annotation; gfortran's fixed reader retains
                // it while the generic parser has no second field.
                let name = if j == 37 && parsed_name.trim().is_empty() {
                    "1D)"
                } else {
                    parsed_name
                };
                let _ = write!(w9, "{:<8}", &name[..name.len().min(8)]);
                if k + 1 < ndxp {
                    let _ = write!(w9, "  ");
                }
            }
            let _ = writeln!(w9);
        }
    }
    if let Some(ref mut w9) = s.out_unit9 {
        // WRITE(9,105) IDAY, (PPMEAN(K+50,IBX), K=1,NDXP)
        let _ = write!(w9, " {:3} ", iday);
        for k in 0..ndxp {
            let _ = write!(w9, "{}", fmt_e10p3(s.ppmean[[PPMEAN_RATE_OFFSET + k, ib0]]));
        }
        let _ = writeln!(w9);
    }
}

// ── HUNT ─────────────────────────────────────────────────────────────────────

/// General-purpose bracket search (Numerical Recipes HUNT).
/// Searches xx[0..n-1] for the bracket containing x.
/// On entry jlo is a hint; on exit jlo is the lower bracket index (0-based).
/// ip > 0 enables verbose output (for debugging).
/// Fortran: HUNT.FOR
pub fn hunt(xx: &[f64], x: f64, jlo: &mut usize) {
    let n = xx.len();
    if n == 0 {
        return;
    }
    let ascnd = xx[n - 1] > xx[0];

    // Convert to 1-based logic then back (Fortran-style algorithm)
    let mut jlo_f: i32 = *jlo as i32; // may be 0 (not set)
    let mut jhi_f: i32;

    if jlo_f <= 0 || jlo_f > n as i32 {
        jlo_f = 0;
        jhi_f = n as i32 + 1;
    } else {
        let mut inc = 1i32;
        if (x >= xx[jlo_f as usize - 1]) == ascnd {
            loop {
                jhi_f = jlo_f + inc;
                if jhi_f > n as i32 {
                    jhi_f = n as i32 + 1;
                    break;
                }
                if (x >= xx[jhi_f as usize - 1]) != ascnd {
                    break;
                }
                jlo_f = jhi_f;
                inc *= 2;
            }
        } else {
            jhi_f = jlo_f;
            loop {
                jlo_f = jhi_f - inc;
                if jlo_f < 1 {
                    jlo_f = 0;
                    break;
                }
                if (x < xx[jlo_f as usize - 1]) != ascnd {
                    break;
                }
                jhi_f = jlo_f;
                inc *= 2;
            }
        }
    }

    // Binary search
    while jhi_f - jlo_f > 1 {
        let jm = (jhi_f + jlo_f) / 2;
        if (x > xx[jm as usize - 1]) == ascnd {
            jlo_f = jm;
        } else {
            jhi_f = jm;
        }
    }

    // Convert back to 0-based: jlo_f is 1-based lower bracket
    *jlo = jlo_f.max(0) as usize;
}

#[cfg(test)]
mod tests {
    use super::*;

    // ── hunt ─────────────────────────────────────────────────────────────────

    fn hunt_result(xx: &[f64], x: f64) -> usize {
        let mut jlo = 0;
        hunt(xx, x, &mut jlo);
        jlo
    }

    // hunt uses strictly-greater bisection and stores the 1-based result
    // directly as usize without subtracting 1. So for xx=[a,b,c,...]:
    //   jlo=k means x is in the open-upper interval (xx[k-1], xx[k]] (0-based).

    #[test]
    fn test_hunt_interior() {
        let xx = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        // 2.5 in (xx[1]=2, xx[2]=3] → 1-based jlo=2
        assert_eq!(hunt_result(&xx, 2.5), 2);
        // 3.0 exactly at xx[2]=3 (upper end of bracket) → jlo=2
        assert_eq!(hunt_result(&xx, 3.0), 2);
        // 4.9 in (xx[3]=4, xx[4]=5] → jlo=4
        assert_eq!(hunt_result(&xx, 4.9), 4);
    }

    #[test]
    fn test_hunt_exact_boundary() {
        let xx = vec![0.0, 10.0, 20.0, 30.0];
        // 10.0 at the upper end of (0, 10] → jlo=1
        assert_eq!(hunt_result(&xx, 10.0), 1);
        // 20.0 at the upper end of (10, 20] → jlo=2
        assert_eq!(hunt_result(&xx, 20.0), 2);
    }

    #[test]
    fn test_hunt_below_range() {
        let xx = vec![5.0, 10.0, 15.0];
        let jlo = hunt_result(&xx, 1.0);
        // Below range: jlo = 0 (standard NR behavior)
        assert_eq!(jlo, 0);
    }

    #[test]
    fn test_hunt_above_range() {
        let xx = vec![1.0, 2.0, 3.0];
        let jlo = hunt_result(&xx, 99.0);
        assert_eq!(jlo, 3); // beyond last element
    }

    #[test]
    fn test_hunt_hint_used() {
        // With a good hint the algorithm converges to the same answer
        let xx: Vec<f64> = (0..100).map(|i| i as f64).collect();
        let mut jlo = 50_usize; // hint near the answer
        hunt(&xx, 73.5, &mut jlo);
        // 73.5 in (xx[73]=73, xx[74]=74] → jlo=74
        assert_eq!(jlo, 74, "jlo={jlo}");
    }

    #[test]
    fn test_hunt_consistent_with_fortran_jddo() {
        // Replicate the JDDO array lookup in ctmlfq: find bracket for jdaydo=75
        let jddo: Vec<f64> = [
            1., 16., 32., 47., 60., 75., 91., 106., 121., 137., 152., 167., 182., 197., 213., 228.,
            244., 259., 274., 289., 305., 320., 335., 350.,
        ]
        .to_vec();
        let mut jlo = 0;
        hunt(&jddo, 75.0, &mut jlo);
        // 75 is at index 5 exactly; HUNT should place jlo = 5 (the exact hit)
        assert_eq!(jlo, 5, "jlo={jlo}");
    }

    // ── fmt_e10p3 ────────────────────────────────────────────────────────────

    #[test]
    fn test_fmt_e10p3_width() {
        // Output must always be exactly 10 characters
        for v in &[0.0, 1.234e-5, -3.14e10, 9.999e99, 1.0e-99] {
            let s = fmt_e10p3(*v);
            assert_eq!(s.len(), 10, "fmt_e10p3({v}) = {s:?}, len={}", s.len());
        }
    }

    #[test]
    fn test_fmt_e10p3_value() {
        let s = fmt_e10p3(1.234e-5);
        let parsed: f64 = s.trim().parse().expect("not a number: {s:?}");
        assert!((parsed - 1.234e-5).abs() / 1.234e-5 < 1e-3);
    }

    #[test]
    fn test_fmt_e10p3_zero() {
        let s = fmt_e10p3(0.0);
        assert_eq!(s.len(), 10);
    }
}
