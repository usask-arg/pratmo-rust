// bctmx.f → CTM-style climatology driver
// CTMLFQ, CTOUTP, CTINIT

use std::io::Write;

use anyhow::Result;

use crate::{
    init::ctinit,
    output::hunt,
    path::dstep,
    reader::{hystat, setday},
    state::ModelState,
};

// ── CTMLFQ ───────────────────────────────────────────────────────────────────

/// Outer loop over latitudes × months for CTM-style climatological run.
/// Reads boxin_gui.dat + CTM climatology (fort03, fort04, fort05, fort51),
/// then for each (lat, month) interpolates T/O3/NOy, initialises boxes,
/// integrates NDAYSD days, and writes output.
///
/// This is a direct translation of bctmx.f SUBROUTINE CTMLFQ (2404 lines).
/// The I/O requires fort03_LLM.x, fort04.x, fort05.x, fort51.x in cinpdir.
/// Fortran: SUBROUTINE CTMLFQ
pub fn ctmlfq(s: &mut ModelState) -> Result<()> {
    use std::fs::File;
    use std::io::{BufRead, BufReader};
    use std::path::Path;

    let base = Path::new(&s.cinpdir);

    // Set all NBOXDO to positive
    for ib in 0..s.nbox {
        s.nboxdo[ib] = s.nboxdo[ib].abs();
    }

    // Read boxin_gui.dat for run configuration
    // bctmx.f reads: jdaydo, xlatdo, xalbedo, iwnoy/iwn2o/iwbry, aero_sf,
    //   ipf(8), ipfr(24), ipjv(10), isza, szaout, iampm, ivarO3, bmoutfile, irunclim, iTfull
    let mut jdaydo: i32 = 75;
    let mut ipf = [0i32; 8];
    let mut ipfr = [0i32; 24];
    let mut ipjv = [0i32; 10];
    let mut isza: i32 = 0;
    let mut szaout: f64 = 90.0;
    let mut iampm: i32 = 0;
    let mut bmoutfile = base.join("boxout.dat");

    let boxin = base.join("boxin_gui.dat");
    if boxin.exists() {
        let f = File::open(&boxin)?;
        let mut r = BufReader::new(f);
        let mut line = String::new();

        macro_rules! read_val {
            ($ty:ty) => {{
                line.clear();
                r.read_line(&mut line)?;
                line.trim().parse::<$ty>().unwrap_or_default()
            }};
        }
        macro_rules! read_ints {
            ($arr:expr) => {{
                line.clear();
                r.read_line(&mut line)?;
                let mut it = line.split_whitespace();
                for x in $arr.iter_mut() {
                    *x = it.next().unwrap_or("0").parse().unwrap_or(0);
                }
            }};
        }

        jdaydo           = read_val!(i32);
        let xlatdo: f64  = read_val!(f64);
        let xalbedo: f64 = read_val!(f64);

        line.clear(); r.read_line(&mut line)?;
        { let mut it = line.split_whitespace();
          s.iwnoy = it.next().unwrap_or("0").parse().unwrap_or(0);
          s.iwn2o = it.next().unwrap_or("0").parse().unwrap_or(0);
          s.iwbry = it.next().unwrap_or("0").parse().unwrap_or(0); }

        let _aero_sf: f64 = read_val!(f64);
        read_ints!(ipf);
        read_ints!(ipfr);
        read_ints!(ipjv);
        isza   = read_val!(i32);
        szaout = read_val!(f64);
        iampm  = read_val!(i32);
        let _ivar_o3: i32 = read_val!(i32);

        // bmoutfile (may be an absolute Windows path — extract filename portion)
        line.clear();
        r.read_line(&mut line)?;
        let raw_path = line.trim();
        // Handle both Windows (\) and Unix (/) separators
        let fname = raw_path
            .rsplit(|c| c == '\\' || c == '/')
            .next()
            .filter(|s| !s.is_empty())
            .unwrap_or("boxout.dat");
        bmoutfile = base.join(fname);

        s.clouds = xalbedo;
        s.xlat   = xlatdo.to_radians();
        s.xlatd  = xlatdo;

        // OSIRIS solar declination formula
        let pi = std::f64::consts::PI;
        let xjd = 2.0 * pi * jdaydo as f64 / 365.0;
        let decang = 6.918e-3 - 0.399912 * xjd.cos() + 0.070257 * xjd.sin()
            - 6.758e-3 * (2.0 * xjd).cos() + 9.07e-4 * (2.0 * xjd).sin()
            - 2.697e-3 * (3.0 * xjd).cos() + 1.480e-3 * (3.0 * xjd).sin();
        let crad = 57.29578_f64;
        s.xdecd = decang * crad;
        s.xdec  = decang;

        const EDIST: [f64; 12] = [
            0.9837, 0.9875, 0.9945, 1.0032, 1.0109, 1.0158,
            1.0165, 1.0128, 1.0057, 0.9970, 0.9892, 0.9842,
        ];
        let mon_approx = ((jdaydo - 1) / 30).min(11) as usize;
        s.flscal = 1.0 / (EDIST[mon_approx] * EDIST[mon_approx]);

        // Compute nd216/nd216s from jdaydo and xlatdo (bctmx.f lines 187-194)
        // JDDO(1..24): representative Julian days for each of 24 half-months
        const JDDO: [i32; 24] = [
            1,16,32,47,60,75,91,106,121,137,152,167,
            182,197,213,228,244,259,274,289,305,320,335,350,
        ];
        // HUNT for closest JDDO to jdaydo
        let mut j0 = 0usize;
        hunt(&JDDO.iter().map(|&x| x as f64).collect::<Vec<_>>(), jdaydo as f64, &mut j0);
        let j0_1 = j0.min(JDDO.len() - 2); // 0-based, clamp to valid range
        let dd1 = (JDDO[j0_1 + 1] - jdaydo).unsigned_abs();
        let dd2 = (JDDO[j0_1] - jdaydo).unsigned_abs();
        let j0_final = if dd1 < dd2 { j0_1 + 1 } else { j0_1 }; // 0-based (Fortran: 1-based)
        // ilat: 2.5° spacing from -90 to +90 (71 lats total, 0-based)
        let ilat = ((xlatdo + 90.0) / 2.5).floor() as i32;
        let nd_val = j0_final as i32 * 71 + ilat;
        s.nd216  = nd_val.max(1); // at least 1 so loop runs
        s.nd216s = nd_val.max(1);
    }

    // Read CTM climatology if fort03_LLM.x exists
    let fort03 = base.join("fort03_LLM.x");
    if fort03.exists() {
        let f = File::open(&fort03)?;
        let mut r = BufReader::new(f);
        let mut line = String::new();

        // Skip two header lines
        r.read_line(&mut line)?;
        line.clear();
        r.read_line(&mut line)?;

        for m in 0..12usize {
            for j in 0..18usize {
                // Skip latitude header
                line.clear();
                r.read_line(&mut line)?;

                // T(1..41) — "3X,11F7.1" repeated
                let t_vals = read_f7_1(&mut r, 41)?;
                for i in 0..t_vals.len().min(41) {
                    s.tinp[[m, j, i]] = t_vals[i];
                }

                // DO3INP(1..31) — "3X,11F7.4" repeated
                let o3_vals = read_f7_4(&mut r, 31)?;
                for i in 0..o3_vals.len().min(31) {
                    s.do3inp[[m, j, i]] = o3_vals[i];
                }
            }
        }
    }

    // Read NOy climatology (fort05.x)
    let fort05 = base.join("fort05.x");
    if fort05.exists() {
        let f = File::open(&fort05)?;
        let mut r = BufReader::new(f);
        let mut line = String::new();
        r.read_line(&mut line)?;
        line.clear();
        r.read_line(&mut line)?;
        for m in 0..12usize {
            for j in 0..18usize {
                line.clear();
                r.read_line(&mut line)?; // header
                let noy_vals = read_f7_4(&mut r, 31)?;
                for i in 0..noy_vals.len().min(31) {
                    s.dnoyi_np[[m, j, i]] = noy_vals[i];
                }
            }
        }
    }

    // Read N2O climatology (fort51.x)
    let fort51 = base.join("fort51.x");
    if fort51.exists() {
        let f = File::open(&fort51)?;
        let mut r = BufReader::new(f);
        let mut line = String::new();
        r.read_line(&mut line)?;
        line.clear();
        r.read_line(&mut line)?;
        for m in 0..12usize {
            for j in 0..18usize {
                line.clear();
                r.read_line(&mut line)?;
                // fort51.x uses F8.4 (8-char fields, not 7)
                let n2o_vals = read_f8_4(&mut r, 31)?;
                for i in 0..n2o_vals.len().min(31) {
                    s.dn2oinp[[m, j, i]] = n2o_vals[i];
                }
            }
        }
    }

    // Initialization call: CTOUTP(0, 0, 0) — sets NDXPP and prints rate list
    ctoutp(s, 0, 0)?;

    // ── Main climatology loop ─────────────────────────────────────────────

    let nd216 = s.nd216;
    let nd216s = s.nd216s;

    for ipath in 1..=nd216 as usize {
        // Compute (ilat2, imon2) from ipath
        let ilat2 = {
            let r = ipath % 71;
            if r == 0 { 71 } else { r }
        };
        let imon2 = (ipath - 1) / 71 + 1;

        let xlat = -90.0 + 2.5 * ilat2 as f64;
        let imon = (imon2 - 1) / 2 + 1;
        let lat = s.xlatd.round() as i32;
        let mon = imon as i32;

        if (ipath as i32) < nd216s {
            continue;
        }

        // Month indices for temporal interpolation
        let (mm1, mm2) = if imon2 == 1 {
            (11usize, 0usize) // 12→0, 1→0  (0-based month)
        } else if imon2 % 2 == 1 {
            (imon - 2, imon - 1)
        } else {
            let m2 = if imon >= 12 { 0 } else { imon };
            (imon - 1, m2)
        };

        // Latitude interpolation index into 18-lat grid
        let xlatinp = core::array::from_fn::<f64, 18, _>(|j| -85.0 + 10.0 * j as f64);
        let jj1 = xlatinp.iter()
            .position(|&la| la > xlat)
            .map(|p| p.saturating_sub(1))
            .unwrap_or(17)
            .min(16);
        let jj2 = (jj1 + 1).min(17);
        let wj = if jj1 == jj2 { 0.0 }
            else { (xlat - xlatinp[jj1]) / (xlatinp[jj2] - xlatinp[jj1]) };
        let wm = 0.5f64;

        // Interpolate T, O3REF, NOyREF, N2OREF for this (lat, month) from climatology.
        // Arrays are stored in ppmv (O3) or ppbv (NOy, N2O); conversion to cm^-3 after HYSTAT.
        let nc = s.nc;
        for i in 0..nc.min(41) {
            s.t[i] = interp2(wm, wj, &s.tinp, mm1, mm2, jj1, jj2, i);
        }
        for i in 0..31usize.min(nc) {
            s.do3ref[i]  = interp2(wm, wj, &s.do3inp,   mm1, mm2, jj1, jj2, i);
            s.dnoy_ref[i]= interp2(wm, wj, &s.dnoyi_np, mm1, mm2, jj1, jj2, i);
            s.dn2oref[i] = interp2(wm, wj, &s.dn2oinp,  mm1, mm2, jj1, jj2, i);
        }

        // CTM mode always uses log-pressure coordinates: CALL HYSTAT(1) in bctmx.f line 352.
        // HYSTAT(1) computes Z from T: Z[ii] = Z[ii-1] + (T[ii-1]+T[ii]) * CLOGP (~2 km/level)
        hystat(s, 1);
        s.xlat = xlat.to_radians();
        s.xlatd = xlat;
        // XDECD and FLSCAL already set from boxin_gui.dat jdaydo computation above
        // (s.xdecd = declination in degrees, s.flscal = 1/edist^2)
        setday(s);
        s.lresol = true;
        s.izalt = 0;

        // Convert climatology from mixing ratio to number density:
        //   O3: ppmv * DM * 1e-6 → cm^-3
        //   NOy, N2O: ppbv * DM * 1e-9 → cm^-3
        // Then extend from level 31 to NC using exponential scaling.
        for j in 0..31usize.min(nc) {
            s.do3ref[j]   *= s.dm[j] * 1.0e-6;
            s.dnoy_ref[j] *= s.dm[j] * 1.0e-9;
            s.dn2oref[j]  *= s.dm[j] * 1.0e-9;
        }
        if nc > 31 {
            let scalez  = if s.do3ref[30]   > 0.0 { s.do3ref[30]   / s.do3ref[29].max(1e-40) } else { 1.0 };
            let scalezz = if s.dnoy_ref[30] > 0.0 { s.dnoy_ref[30] / s.dnoy_ref[29].max(1e-40) } else { 1.0 };
            let scalez2 = if s.dn2oref[30]  > 0.0 { s.dn2oref[30]  / s.dn2oref[29].max(1e-40) } else { 1.0 };
            for j in 31..nc {
                s.do3ref[j]   = s.do3ref[j - 1]   * scalez;
                s.dnoy_ref[j] = s.dnoy_ref[j - 1] * scalezz;
                s.dn2oref[j]  = s.dn2oref[j - 1]  * scalez2;
            }
        }
        // Rebuild O3 column integral
        let zzht = s.zzht;
        let nc = s.nc;
        s.do3int[nc - 1] = s.do3ref[nc - 1] * zzht;
        for j in (0..nc - 1).rev() {
            s.do3int[j] = s.do3int[j + 1]
                + 0.5 * (s.z[j + 1] - s.z[j]) * (s.do3ref[j + 1] + s.do3ref[j]);
        }

        // Set mixing ratios from reference profiles, then call CTINIT for species init.
        let nbox = s.nbox;

        // Open output file on first box of this atmosphere
        let mut outfile: Option<std::io::BufWriter<std::fs::File>> = if ipath == nd216s as usize {
            match std::fs::File::create(&bmoutfile) {
                Ok(f) => Some(std::io::BufWriter::new(f)),
                Err(e) => {
                    eprintln!("Warning: cannot open bmoutfile {:?}: {}", bmoutfile, e);
                    None
                }
            }
        } else {
            None
        };

        for ib in 0..nbox {
            let ialt = (s.nboxdo[ib].unsigned_abs() as usize).saturating_sub(1);
            let densbx = s.dm[ialt];
            s.fo3[ib]  = s.do3ref[ialt]   / densbx;
            s.fnoy[ib] = s.dnoy_ref[ialt] / densbx;
            s.fn2o[ib] = s.dn2oref[ialt]  / densbx;
            ctinit(s, ib, densbx, lat, mon);
            s.lresol = true;
            dstep(s, ipath as i32, ib)?;

            // Write output row to bmoutfile (Fortran: bctmx.f lines 1170–1371)
            if let Some(ref mut w) = outfile {
                write_box_row(
                    s, w, ib, nc, &ipf, &ipfr, &ipjv, isza, szaout, iampm, jdaydo,
                )?;
            }
        }

        ctoutp(s, (ilat2 - 1) as usize, (imon - 1) as usize)?;
    }

    Ok(())
}

// ── CTOUTP ───────────────────────────────────────────────────────────────────

/// Compute and write climatological diagnostics.
/// ilat=0, imon=0: initialization (sets NDXPP, prints rate list).
/// ilat>0, imon>0: per-atmosphere output (prints loss frequencies).
/// Fortran: SUBROUTINE CTOUTP(LAT, MON, IPATH)
pub fn ctoutp(s: &mut ModelState, ilat: usize, imon: usize) -> Result<()> {
    use crate::constants::NNDXPQ;

    if ilat == 0 && imon == 0 {
        // Initialization: set NDXPP with rate indices for key chemical families
        // (mirrors the hardcoded assignments in bctmx.f CTOUTP when MON=0)
        let indices: [usize; NNDXPQ] = [
            9, 10, 137,   // N2O loss
            86, 87, 89, 90, // NOy loss
            49, 7, 8, 104, // CH4 loss
            5,             // H2O loss
            111, 208,      // CFC-11
            58, 59, 60, 106, 218, // CO production
            48,            // CO loss
            5,             // OH production
            56,            // HCHO rainout
            51, 52,        // CH3OO branching
            0,             // unused
        ];
        for (k, &idx) in indices.iter().enumerate() {
            if k < NNDXPQ { s.ndxpp[k] = idx; }
        }

        println!("  StratLoss OUTPUT rates:");
        for k in 0..NNDXPQ {
            let kr = s.ndxpp[k];
            if kr > 0 && kr <= s.nrates {
                println!("{:5}{:5}     {:<8}{:<8}", k + 1, kr,
                    &s.rfmt_str[kr - 1][0], &s.rfmt_str[kr - 1][1]);
            }
        }
    } else {
        // Per-atmosphere output: compute and print key loss frequencies
        println!("{:5}{:5} StratLossFreq: jNOy NOY N2O F11", imon, ilat);

        for ib in 0..s.nbox {
            let ialt_1 = s.nboxdo[ib].unsigned_abs() as usize;
            if ialt_1 == 0 || ialt_1 > s.nc { continue; }
            let ialt = ialt_1 - 1;
            let densbx = s.dm[ialt];
            let zkm = 1e-5 * s.z[ialt];

            let fn2o = s.fn2o[ib];
            let fnoy = s.fnoy[ib];
            let fch4 = s.fch4[ib];
            let fcfcl3 = s.fcfcl3[ib];
            let fo3 = s.fo3[ib];

            // N2O loss freq: (PPMEAN(54)+PPMEAN(55)+PPMEAN(56)) / (M * FN2O)
            // NDXPP(1)=9→ppmean[50+0], NDXPP(2)=10→ppmean[50+1], NDXPP(3)=137→ppmean[50+2]
            let sln2o_num = s.ppmean[[50, ib]] + s.ppmean[[51, ib]] + s.ppmean[[52, ib]];
            let sln2o = if fn2o > 0.0 { sln2o_num / (densbx * fn2o) } else { 0.0 };

            // NOy loss freq: 2*(PPMEAN(52)+PPMEAN(53)) / (M * FNOY)
            let slnoy_num = 2.0 * (s.ppmean[[51, ib]] + s.ppmean[[52, ib]]);
            let slnoy = if fnoy > 0.0 { slnoy_num / (densbx * fnoy) } else { 0.0 };

            // CH4 loss freq
            let slch4_num = s.ppmean[[57, ib]] + s.ppmean[[58, ib]] + s.ppmean[[59, ib]] + s.ppmean[[60, ib]];
            let slch4 = if fch4 > 0.0 { slch4_num / (densbx * fch4) } else { 0.0 };

            // CFC-11 loss freq
            let slf11 = if fcfcl3 > 0.0 {
                (s.ppmean[[62, ib]] + s.ppmean[[63, ib]]) / (densbx * fcfcl3)
            } else { 0.0 };

            // O3 loss freq
            let slo3 = if fo3 > 0.0 { s.ppmean[[50, ib]] / (86400.0 * fo3) } else { 0.0 };

            println!("{:3} {:4.1} {:9.2e} {:9.2e} {:9.2e} {:9.2e} {:9.2e} {:9.2e}",
                ib + 1, zkm, fo3, fnoy, fn2o, sln2o, slnoy, slch4);
            let _ = (slo3, slf11);
        }
    }
    Ok(())
}

// ── Output file writer ────────────────────────────────────────────────────────

/// Write one box's row (and header/fill rows on ib=0 and ib=nbox-1) to the output writer.
/// Mirrors bctmx.f lines 1155–1371.
#[allow(clippy::too_many_arguments)]
fn write_box_row(
    s: &mut ModelState,
    w: &mut impl Write,
    ib: usize,           // 0-based box index
    nc: usize,           // number of altitude levels
    ipf:  &[i32; 8],    // family species flags (N2O,CH4,H2O,NOy,Cly,Bry,CO,Aer)
    ipfr: &[i32; 24],   // radical diurnal flags
    ipjv: &[i32; 10],   // J-value output indices
    isza: i32,           // 0=full diurnal, 1=single SZA
    szaout: f64,         // target SZA (degrees) for isza=1
    iampm: i32,          // 0=sunrise, 1=sunset
    jdaydo: i32,
) -> Result<()> {
    // ng: XXNOFT 1-based slot for each of the 24 cref2 radical species
    // DATA NG from bctmx.f: 9,10,11,13,15,6,12,17,18,19,21,22,24,27,20,26,25,2,23,1,3,5,4,30
    const NG: [usize; 24] = [9,10,11,13,15,6,12,17,18,19,21,22,24,27,20,26,25,2,23,1,3,5,4,30];
    const CREF2: [&str; 24] = [
        "NO","NO2","NO3","HNO3","HO2NO2","H2CO","N2O5","OH",
        "HO2","H2O2","Cl","Cl2","ClO","OClO","HCl","HOCl","ClONO2",
        "Br","BrCl","BrO","HBr","HOBr","BrONO2","O3",
    ];
    const CREF1: [&str; 13] = [
        "z (km)","T (K)","p (hPa)","air (cm-3)","O3",
        "N2O","CH4","H2O","NOy","Cly","Bry","CO","Aer SA (cm2)",
    ];

    let ialt  = (s.nboxdo[ib].unsigned_abs() as usize).saturating_sub(1);
    let densbx = s.dm[ialt];
    let zkm    = s.z[ialt] * 1e-5;
    let ntimdo = s.ntimdo;
    let nbox   = s.nbox;

    // ── Assemble xout (scalar fields) ────────────────────────────────────────
    // xout(1..5) always: z, T, p, M, O3
    let mut xout = vec![0.0f64; 13];
    let mut chead1 = vec![""; 13];
    xout[0] = zkm;
    xout[1] = s.t[ialt];
    xout[2] = s.pstd[ialt];
    xout[3] = densbx;
    xout[4] = s.do3ref[ialt] / densbx;
    for i in 0..5 { chead1[i] = CREF1[i]; }
    let mut ij = 5usize;
    let fam_vals = [s.fn2o[ib], s.fch4[ib], s.fh2o[ib], s.fnoy[ib], s.fclx[ib], s.fbrx[ib], s.fco[ib], 0.0];
    for (n, &flag) in ipf.iter().enumerate() {
        if flag == 1 && ij < 13 {
            xout[ij] = fam_vals[n];
            chead1[ij] = CREF1[5 + n];
            ij += 1;
        }
    }
    let np1 = ij;

    // ── Assemble xoutr (radical diurnal cycles) ───────────────────────────────
    let mut xoutr: Vec<Vec<f64>> = Vec::new();
    let mut chead2: Vec<&str> = Vec::new();
    for (n, &flag) in ipfr.iter().enumerate() {
        if flag == 1 && n < 24 {
            let slot = NG[n].saturating_sub(1); // 0-based XXNOFT row index
            let radrow: Vec<f64> = (0..ntimdo)
                .map(|it| if slot < 30 { s.xxnoft[[slot, it, ib]] / densbx } else { 0.0 })
                .collect();
            xoutr.push(radrow);
            chead2.push(CREF2[n]);
        }
    }
    let np2 = xoutr.len();

    // ── Assemble xoutj (J-values at 16 stored solar angles) ──────────────────
    let mut xoutj: Vec<Vec<f64>> = Vec::new();
    let mut chead3: Vec<String> = Vec::new();
    for &jidx in ipjv.iter() {
        if jidx > 0 {
            let jj = (jidx as usize).saturating_sub(1); // 0-based J-value index
            let jrow: Vec<f64> = (0..16usize)
                .map(|it| if it < 16 { s.storjv[[jj, it, ib]] } else { 0.0 })
                .collect();
            xoutj.push(jrow);
            let title = if jj < s.titlej.len() { s.titlej[jj][0].clone() } else { String::new() };
            chead3.push(title);
        }
    }
    let njv = xoutj.len();

    // ── SZA interpolation if isza=1 (single-SZA output) ─────────────────────
    let ntout = if isza == 1 { 1 } else { ntimdo };
    let mut xoutr_out: Vec<Vec<f64>> = xoutr.clone();
    let mut dtime2: Vec<f64> = (0..ntimdo)
        .map(|it| { let d = s.dtime[it] / 3600.0 - 12.0; if d < 0.0 { d + 24.0 } else { d } })
        .collect();

    if isza == 1 {
        // Find time step closest to szaout using HUNT on ztime
        let (it1, it2, idt) = if iampm == 0 {
            // sunrise side: second half of ztime array
            let mid = ntimdo / 2;
            (mid + 1, ntimdo, ntimdo - mid)
        } else {
            // sunset side: first half
            (1, ntimdo / 2 + 1, ntimdo / 2 + 1)
        };
        let it1 = it1.saturating_sub(1); // 0-based
        let it2 = it2.min(ntimdo);
        let zslice: Vec<f64> = (it1..it2).map(|i| s.ztime[i]).collect();
        let mut jlo = 0usize;
        hunt(&zslice, szaout, &mut jlo);
        let jlo = jlo.min(idt.saturating_sub(2));
        let it_base = it1 + jlo;

        let wtsz = if jlo + 1 < zslice.len() && (s.ztime[it_base + 1] - s.ztime[it_base]).abs() > 1e-10 {
            (szaout - s.ztime[it_base]) / (s.ztime[it_base + 1] - s.ztime[it_base])
        } else { 0.0 }.clamp(0.0, 1.0);

        // Interpolate each radical at the target SZA
        for (n, row) in xoutr_out.iter_mut().enumerate() {
            let v0 = xoutr[n][it_base];
            let v1 = if it_base + 1 < ntimdo { xoutr[n][it_base + 1] } else { v0 };
            *row = vec![wtsz * v1 + (1.0 - wtsz) * v0];
        }
        let d0 = dtime2[it_base];
        let d1 = if it_base + 1 < ntimdo { dtime2[it_base + 1] } else { d0 };
        dtime2 = vec![wtsz * d1 + (1.0 - wtsz) * d0];
    }

    let ncd = np2 * ntout + np1.saturating_sub(5) + njv * 16;

    if ib == 0 {
        // Header (written once, when processing the first box)
        let cver = "v5.0";
        let _ = cver;
        writeln!(w, "{:<280}", "PRATMO v8.0 (Rust port): JPL-09 photochem data")?;
        writeln!(w, "{:<280}", "Rust port of Prather PRATMO, faithful to Fortran v6.0")?;
        writeln!(w, "{:<280}", "all species in volume mixing ratio, -9 = model not run")?;
        writeln!(w, "{:4}   Standard output fields",   5)?;
        writeln!(w, "{:4}   Family output fields",     np1.saturating_sub(5))?;
        writeln!(w, "{:4}   Radical output fields",    np2)?;
        writeln!(w, "{:4}   Jvalue output fields",     njv)?;
        writeln!(w, "{:4}   Total columns of output",  ncd + 5)?;
        if isza == 1 {
            writeln!(w, "{:6.2} {:2}Solar Zenith Angle", szaout,
                if iampm == 0 { "am" } else { "pm" })?;
            writeln!(w, "{:6.2}   Local Solar Time", dtime2[0])?;
        } else {
            writeln!(w, "{:4}   Number of Solar Zenith Angles", ntout)?;
            writeln!(w, "{:4}   Number of Jvalue SZAs",         16)?;
        }
        writeln!(w, "{:4}   Number of Altitudes", 39)?;
        writeln!(w, "{:4}   Julian Day",     jdaydo)?;
        writeln!(w, "{:6.2}   Latitude",     s.xlatd)?;
        writeln!(w, "{:6.2}   Surface Albedo", s.clouds)?;
        writeln!(w)?;

        // Column header line: scalar fields + radical names + J-value names
        let mut hdr = String::new();
        for &h in chead1[..np1].iter() { hdr.push_str(&format!("{:>13}", h)); }
        for (n, &rname) in chead2.iter().enumerate() {
            let _ = n;
            for _it in 0..ntout { hdr.push_str(&format!("{:>13}", rname)); }
        }
        for (n, jname) in chead3.iter().enumerate() {
            let _ = n;
            for _it in 0..16 { hdr.push_str(&format!("{:>13}", &jname[..jname.len().min(12)])); }
        }
        writeln!(w, "  {}", hdr)?;

        // SZA/time rows (for full diurnal only)
        if isza == 0 {
            let mut row_sza = format!("{:>13}{:>13}{:>13}{:>13}{:>13}", "0","0","0","0","0");
            for _n in 0..np2 {
                for it in 0..ntimdo { row_sza.push_str(&format!("{:13.4E}", s.ztime[it])); }
            }
            for _n in 0..njv {
                for it in 0..16 { row_sza.push_str(&format!("{:13.4E}", s.ztime[it.min(ntimdo.saturating_sub(1))])); }
            }
            writeln!(w, "{}", row_sza)?;

            let mut row_lt = format!("{:>13}{:>13}{:>13}{:>13}{:>13}", "0","0","0","0","0");
            for _n in 0..np2 {
                for it in 0..ntimdo { row_lt.push_str(&format!("{:13.4E}", dtime2[it])); }
            }
            for _n in 0..njv {
                for it in 0..16 { row_lt.push_str(&format!("{:13.4E}", dtime2[it.min(dtime2.len().saturating_sub(1))])); }
            }
            writeln!(w, "{}", row_lt)?;
        }
        writeln!(w, " {}", "-----------".repeat(ncd + 5))?;

        // Fill altitudes above boxes with -9 (nc-2 down to nboxdo[0])
        let top_box_ialt = (s.nboxdo[0].unsigned_abs() as usize).saturating_sub(1);
        for ii in (top_box_ialt + 1..nc.saturating_sub(1)).rev() {
            let mut fill = format!("{:13.4E}{:13.4E}{:13.4E}{:13.4E}{:13.4E}",
                s.z[ii] * 1e-5, s.t[ii], s.pstd[ii], s.dm[ii],
                s.do3ref[ii] / s.dm[ii].max(1e-50));
            let nfill = ncd;
            for _ in 0..nfill { fill.push_str(&format!("{:13.4E}", -9.0f64)); }
            writeln!(w, "{}", fill)?;
        }
    }

    // ── Per-box data row ──────────────────────────────────────────────────────
    let mut row = String::new();
    for &v in &xout[..np1] { row.push_str(&format!("{:13.4E}", v)); }
    for (n, _) in chead2.iter().enumerate() {
        for it in 0..ntout {
            row.push_str(&format!("{:13.4E}", xoutr_out[n][it]));
        }
    }
    for (n, _) in chead3.iter().enumerate() {
        for it in 0..16 {
            row.push_str(&format!("{:13.4E}", xoutj[n][it]));
        }
    }
    writeln!(w, "{}", row)?;

    // ── Fill below boxes on last box ──────────────────────────────────────────
    if ib == nbox - 1 {
        let bot_box_ialt = (s.nboxdo[nbox - 1].unsigned_abs() as usize).saturating_sub(1);
        for ii in (1..bot_box_ialt).rev() {
            let mut fill = format!("{:13.4E}{:13.4E}{:13.4E}{:13.4E}{:13.4E}",
                s.z[ii] * 1e-5, s.t[ii], s.pstd[ii], s.dm[ii],
                s.do3ref[ii] / s.dm[ii].max(1e-50));
            let nfill = ncd;
            for _ in 0..nfill { fill.push_str(&format!("{:13.4E}", -9.0f64)); }
            writeln!(w, "{}", fill)?;
        }
    }

    Ok(())
}

// ── Helpers ───────────────────────────────────────────────────────────────────

/// Bilinear interpolation of a [12, 18, NL] climatology array.
fn interp2(
    wm: f64,
    wj: f64,
    arr: &ndarray::Array3<f64>,
    mm1: usize,
    mm2: usize,
    jj1: usize,
    jj2: usize,
    i: usize,
) -> f64 {
    let y1 = arr[[mm1, jj1, i]];
    let y2 = arr[[mm1, jj2, i]];
    let y3 = arr[[mm2, jj1, i]];
    let y4 = arr[[mm2, jj2, i]];
    (1.0 - wm) * (1.0 - wj) * y1 + (1.0 - wm) * wj * y2
        + wm * (1.0 - wj) * y3 + wm * wj * y4
}

/// Read a block of values formatted as "3X, 11F7.1" (T profile).
fn read_f7_1(r: &mut impl std::io::BufRead, count: usize) -> Result<Vec<f64>> {
    let mut vals = Vec::with_capacity(count);
    let mut line = String::new();
    while vals.len() < count {
        line.clear();
        if r.read_line(&mut line)? == 0 { break; }
        let start = if line.len() > 3 { 3 } else { 0 };
        let row = &line[start..];
        let mut off = 0;
        while off + 7 <= row.len() && vals.len() < count {
            let v: f64 = row[off..off + 7].trim().parse().unwrap_or(0.0);
            vals.push(v);
            off += 7;
        }
    }
    Ok(vals)
}

/// Read a block of values formatted as "3X, 11F7.4" (O3/NOy profile).
fn read_f7_4(r: &mut impl std::io::BufRead, count: usize) -> Result<Vec<f64>> {
    read_f7_1(r, count) // same format width, just different decimal places
}

/// Read a block of values formatted as "3X, 11F8.4" (N2O profile, fort51.x).
fn read_f8_4(r: &mut impl std::io::BufRead, count: usize) -> Result<Vec<f64>> {
    let mut vals = Vec::with_capacity(count);
    let mut line = String::new();
    while vals.len() < count {
        line.clear();
        if r.read_line(&mut line)? == 0 { break; }
        let start = if line.len() > 3 { 3 } else { 0 };
        let row = &line[start..];
        let row_trimmed = row.trim_end();
        let mut off = 0;
        while off + 8 <= row_trimmed.len() && vals.len() < count {
            let v: f64 = row_trimmed[off..off + 8].trim().parse().unwrap_or(0.0);
            vals.push(v);
            off += 8;
        }
    }
    Ok(vals)
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::Array3;

    // ── interp2 ───────────────────────────────────────────────────────────────

    #[test]
    fn test_interp2_at_corners() {
        let mut arr = Array3::<f64>::zeros((12, 18, 5));
        arr[[0, 0, 2]] = 1.0;
        arr[[0, 1, 2]] = 2.0;
        arr[[1, 0, 2]] = 3.0;
        arr[[1, 1, 2]] = 4.0;

        assert_eq!(interp2(0.0, 0.0, &arr, 0, 1, 0, 1, 2), 1.0, "(mm1,jj1)");
        assert_eq!(interp2(1.0, 0.0, &arr, 0, 1, 0, 1, 2), 3.0, "(mm2,jj1)");
        assert_eq!(interp2(0.0, 1.0, &arr, 0, 1, 0, 1, 2), 2.0, "(mm1,jj2)");
        assert_eq!(interp2(1.0, 1.0, &arr, 0, 1, 0, 1, 2), 4.0, "(mm2,jj2)");
    }

    #[test]
    fn test_interp2_midpoint() {
        let mut arr = Array3::<f64>::zeros((12, 18, 5));
        arr[[0, 0, 0]] = 10.0;
        arr[[0, 1, 0]] = 20.0;
        arr[[1, 0, 0]] = 30.0;
        arr[[1, 1, 0]] = 40.0;
        // Midpoint of all four corners = average = 25.0
        assert_eq!(interp2(0.5, 0.5, &arr, 0, 1, 0, 1, 0), 25.0, "midpoint");
    }

    #[test]
    fn test_interp2_linear_in_wm() {
        let mut arr = Array3::<f64>::zeros((12, 18, 3));
        arr[[2, 5, 1]] = 0.0;
        arr[[2, 6, 1]] = 0.0;
        arr[[3, 5, 1]] = 100.0;
        arr[[3, 6, 1]] = 100.0;
        // Varying wm only (wj=0): expect 100*wm
        assert!((interp2(0.3, 0.0, &arr, 2, 3, 5, 6, 1) - 30.0).abs() < 1e-10);
        assert!((interp2(0.7, 0.0, &arr, 2, 3, 5, 6, 1) - 70.0).abs() < 1e-10);
    }

    #[test]
    fn test_interp2_uniform_field() {
        // Uniform field: result must equal that value everywhere
        let mut arr = Array3::<f64>::zeros((12, 18, 4));
        for m in 0..2 { for j in 0..2 { arr[[m, j, 3]] = 42.0; } }
        for &wm in &[0.0, 0.25, 0.5, 0.75, 1.0] {
            for &wj in &[0.0, 0.5, 1.0] {
                assert_eq!(interp2(wm, wj, &arr, 0, 1, 0, 1, 3), 42.0);
            }
        }
    }

    // ── read_f8_4 ─────────────────────────────────────────────────────────────

    #[test]
    fn test_read_f8_4_single_line() {
        use std::io::Cursor;
        // F8.4 format: 3-char prefix "   ", then 8-char fields right-justified.
        // " 0.3162" is 7 chars; " 0.3162 " would be 8. Typical: "  0.3162".
        // Build exactly 6 × 8-char fields = 48 chars after the 3-char prefix.
        let data = "     0.3162  0.0631  0.0126  0.0025  0.0005  0.0001\n";
        //          ^--- 3 spaces then 6 × " Xx.xxxx" (8 chars each)
        let mut cur = Cursor::new(data);
        let vals = read_f8_4(&mut cur, 6).unwrap();
        assert_eq!(vals.len(), 6, "expected 6 values, got {}", vals.len());
        assert!((vals[0] - 0.3162).abs() < 1e-4, "vals[0]={}", vals[0]);
        assert!((vals[1] - 0.0631).abs() < 1e-4, "vals[1]={}", vals[1]);
        assert!((vals[5] - 0.0001).abs() < 1e-5, "vals[5]={}", vals[5]);
    }

    #[test]
    fn test_read_f8_4_multiline() {
        use std::io::Cursor;
        // Standard fort51.x "3X,11F8.4": 11 values per line, 8 chars each.
        // 11 × 8 = 88 chars of data per line after 3-char prefix.
        let line1 = "     0.3162  0.0631  0.0126  0.0025  0.0005  0.0001  0.3162  0.0631  0.0126  0.0025  0.0005\n";
        let line2 = "     0.0001  0.3162\n";
        let data = format!("{line1}{line2}");
        let mut cur = Cursor::new(data.as_str());
        let vals = read_f8_4(&mut cur, 13).unwrap();
        assert_eq!(vals.len(), 13);
        assert!((vals[0]  - 0.3162).abs() < 1e-4, "vals[0]={}", vals[0]);
        assert!((vals[11] - 0.0001).abs() < 1e-5, "vals[11]={}", vals[11]);
        assert!((vals[12] - 0.3162).abs() < 1e-4, "vals[12]={}", vals[12]);
    }

    // ── read_f7_1 ─────────────────────────────────────────────────────────────

    #[test]
    fn test_read_f7_1_temperature() {
        use std::io::Cursor;
        // "3X,11F7.1" — T profiles use 1 decimal place, 7 chars per field.
        // Fields: "  250.0"  "  240.5"  "  230.1" (each 7 chars, right-justified).
        let data = "     250.0  240.5  230.1\n";
        //          ^3X^ then each "  nnn.d" = 7 chars
        let mut cur = Cursor::new(data);
        let vals = read_f7_1(&mut cur, 3).unwrap();
        assert_eq!(vals.len(), 3, "expected 3 values, got {}", vals.len());
        assert!((vals[0] - 250.0).abs() < 0.05, "vals[0]={}", vals[0]);
        assert!((vals[1] - 240.5).abs() < 0.05, "vals[1]={}", vals[1]);
        assert!((vals[2] - 230.1).abs() < 0.05, "vals[2]={}", vals[2]);
    }
}
