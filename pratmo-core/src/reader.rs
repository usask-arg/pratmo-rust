// bread.f → I/O module
// ModelReader trait + FortranReader implementation

use std::collections::HashMap;
use std::io::{BufRead, BufReader, Cursor};
use std::path::Path;

use anyhow::{bail, Context, Result};

use crate::constants::NDEN;
use crate::solver::fixmix;
use crate::state::ModelState;

// ── Trait ────────────────────────────────────────────────────────────────────

/// Abstraction over input file sources so alternative formats (TOML, JSON,
/// netCDF) can implement this trait without touching model logic.
pub trait ModelReader {
    fn read_spectral_data(&mut self, s: &mut ModelState) -> Result<()>;
    fn read_rate_constants(&mut self, s: &mut ModelState) -> Result<()>;
    fn read_atmosphere(&mut self, s: &mut ModelState) -> Result<()>;
    fn read_ozone_profile(&mut self, s: &mut ModelState) -> Result<()>;
    fn read_control_params(&mut self, s: &mut ModelState) -> Result<()>;
    fn read_jh2o(&mut self, s: &mut ModelState) -> Result<()>;

    /// Read initial species mixing ratios (fort02.x) for non-CTM (DIURN) mode.
    fn read_initial_densities(&mut self, s: &mut ModelState) -> Result<()>;

    /// Read everything in READIN order.
    fn read_all(&mut self, s: &mut ModelState) -> Result<()> {
        self.read_spectral_data(s)?;
        self.read_rate_constants(s)?;
        self.read_control_params(s)?;
        self.read_atmosphere(s)?;
        self.read_ozone_profile(s)?;
        self.read_jh2o(s)?;
        // Fort02.x: initial mixing ratios for DIURN/TPATH mode (skipped in CTM mode).
        // Fortran bread.f reads fort02.x BEFORE the 1.049/1.066/280K test-scaling block,
        // so initial densities are rescaled by the pre-1.049 DM.
        if s.nd216 <= 0 {
            let _ = self.read_initial_densities(s); // optional — fort02.x may not exist
        }
        // Apply calibration scaling (bread.f lines 396–402): AFTER fort02.x initial densities.
        // Fortran: do3ref/do3int *= 1.066; dm/do2int *= 1.049
        // NOTE: the original Fortran test block also set t = 280.0 (a flat override for
        // sensitivity testing only). That line was ported unconditionally, clobbering the
        // real US STD temperature profile loaded from fort13.x. It is removed here.
        let nc = s.nc;
        for j in 0..nc {
            s.do3ref[j] *= 1.066;
            s.do3int[j] *= 1.066;
            s.dm[j]     *= 1.049;
            s.do2int[j] *= 1.049;
        }
        Ok(())
    }
}

// ── FortranReader ─────────────────────────────────────────────────────────────

// ── Embedded default data ─────────────────────────────────────────────────────

static EMBEDDED_FORT01:  &[u8] = include_bytes!("../data/fort01.x");
static EMBEDDED_FORT02:  &[u8] = include_bytes!("../data/fort02.x");
static EMBEDDED_FORT10:  &[u8] = include_bytes!("../data/fort10_cam06.x");
static EMBEDDED_FORT11:  &[u8] = include_bytes!("../data/fort11_jpl09.x");
static EMBEDDED_FORT13:  &[u8] = include_bytes!("../data/fort13.x");
static EMBEDDED_FORT14:  &[u8] = include_bytes!("../data/fort14.x");
static EMBEDDED_JH2O:    &[u8] = include_bytes!("../data/J_H2O_SZA0.dat");

// ── FortranReader ─────────────────────────────────────────────────────────────

/// Reads Fortran fixed-format files from a directory, with optional per-file
/// byte overrides (used by [`FortranReader::embedded`] to serve compiled-in data).
pub struct FortranReader {
    pub input_dir: std::path::PathBuf,
    overrides: HashMap<&'static str, &'static [u8]>,
}

impl FortranReader {
    /// File-based reader: all files are read from `input_dir`.
    pub fn new(input_dir: impl AsRef<Path>) -> Self {
        Self {
            input_dir: input_dir.as_ref().to_owned(),
            overrides: HashMap::new(),
        }
    }

    /// Embedded reader: core science data is served from compiled-in bytes.
    /// Falls back to `input_dir` (empty) for any file not in the override set,
    /// which is the right behaviour for optional files (fort03, fort05, fort51).
    pub fn embedded() -> Self {
        let mut overrides: HashMap<&'static str, &'static [u8]> = HashMap::new();
        overrides.insert("fort01.x",        EMBEDDED_FORT01);
        overrides.insert("fort02.x",        EMBEDDED_FORT02);
        overrides.insert("fort10_cam06.x",  EMBEDDED_FORT10);
        overrides.insert("fort11_jpl09.x",  EMBEDDED_FORT11);
        overrides.insert("fort13.x",        EMBEDDED_FORT13);
        overrides.insert("fort14.x",        EMBEDDED_FORT14);
        overrides.insert("J_H2O_SZA0.dat",  EMBEDDED_JH2O);
        Self { input_dir: std::path::PathBuf::new(), overrides }
    }

    fn open(&self, name: &str) -> Result<Box<dyn BufRead>> {
        if let Some(bytes) = self.overrides.get(name) {
            return Ok(Box::new(BufReader::new(Cursor::new(*bytes))));
        }
        let path = self.input_dir.join(name);
        let f = std::fs::File::open(&path)
            .with_context(|| format!("cannot open {}", path.display()))?;
        Ok(Box::new(BufReader::new(f)))
    }
}

// ── Fixed-format parsing helpers ─────────────────────────────────────────────

/// Parse a line of up to `n` whitespace-separated f64 values.
fn parse_f64s(line: &str) -> Vec<f64> {
    line.split_whitespace()
        .filter_map(|t| t.parse::<f64>().ok())
        .collect()
}

/// Parse fixed-width chunks of `width` chars starting at `offset`, as i32.
fn parse_fixed_i32s(line: &str, offset: usize, width: usize, count: usize) -> Vec<i32> {
    let bytes = line.as_bytes();
    let mut out = Vec::with_capacity(count);
    for i in 0..count {
        let start = offset + i * width;
        let end = (start + width).min(bytes.len());
        if start >= bytes.len() {
            out.push(0);
            continue;
        }
        let chunk = std::str::from_utf8(&bytes[start..end]).unwrap_or("0").trim();
        out.push(chunk.parse::<i32>().unwrap_or(0));
    }
    out
}

/// Read lines from `reader` collecting all f64 tokens until `target` count reached.
/// Handles values that span multiple lines (Fortran implied-DO loops).
fn read_f64_array(reader: &mut impl BufRead, target: usize) -> Result<Vec<f64>> {
    let mut vals: Vec<f64> = Vec::with_capacity(target);
    while vals.len() < target {
        let mut line = String::new();
        if reader.read_line(&mut line)? == 0 {
            bail!("unexpected EOF reading f64 array (got {}, need {})", vals.len(), target);
        }
        for tok in line.split_whitespace() {
            if let Ok(v) = tok.parse::<f64>() {
                vals.push(v);
                if vals.len() == target {
                    break;
                }
            }
        }
    }
    Ok(vals)
}

/// Read exactly one non-empty text line.
fn next_line(reader: &mut impl BufRead) -> Result<String> {
    let mut line = String::new();
    reader.read_line(&mut line)?;
    Ok(line.trim_end_matches('\n').trim_end_matches('\r').to_owned())
}

/// Skip one line.
fn skip_line(reader: &mut impl BufRead) -> Result<()> {
    let mut buf = String::new();
    reader.read_line(&mut buf)?;
    Ok(())
}

// ── fort10_cam06.x — spectral data ───────────────────────────────────────────

impl ModelReader for FortranReader {
    /// Translates the fort10 read section of READIN (bread.f lines 22–87).
    fn read_spectral_data(&mut self, s: &mut ModelState) -> Result<()> {
        let mut r = self.open("fort10_cam06.x")?;

        // Line 1: TITLE0
        s.title0 = next_line(&mut r)?;
        eprintln!("{}", s.title0);

        // Line 2: `10X,14I5` — NQQQ NWWW NW1 NW2 NWSRB NSR NODF
        // Format: 10 chars skipped, then 14 five-char integers.
        // Actual file: "NW--QQQQQQ   40   77    1   77    0   15    6"
        let line2 = next_line(&mut r)?;
        let ints = parse_fixed_i32s(&line2, 10, 5, 14);
        let nqqq  = ints[0] as usize;
        let nwww  = ints[1] as usize;
        s.nw1    = ints[2] as usize;
        s.nw2    = ints[3] as usize;
        s.nwsrb  = ints[4] as usize;
        s.nsr    = ints[5] as usize;
        s.nodf   = ints[6] as usize;
        s.njval  = nqqq + 4;

        // Line 3: `8E10.3` → QMI1M, QMIWVL
        let line3 = next_line(&mut r)?;
        let fv = parse_f64s(&line3);
        s.qmi1m  = *fv.get(0).unwrap_or(&0.0);
        s.qmiwvl = *fv.get(1).unwrap_or(&0.0);

        // Line 4: title/blank
        skip_line(&mut r)?;

        // WBIN(1..nwww+1) — wavelength bin edges
        let wbin = read_f64_array(&mut r, nwww + 1)?;
        for (i, v) in wbin.iter().enumerate() {
            s.wbin[i] = *v;
        }
        for iw in 0..nwww {
            s.wl[iw] = 0.5 * (s.wbin[iw] + s.wbin[iw + 1]);
        }

        // blank line
        skip_line(&mut r)?;

        // FL(1..nwww) — solar flux
        let fl = read_f64_array(&mut r, nwww)?;
        for (i, v) in fl.iter().enumerate() {
            s.fl[i] = *v;
        }

        // blank line
        skip_line(&mut r)?;

        // NSR lines: `6F10.1,I3` → ODF(I,L), ISR(L)
        // Format: nodf floats of 10 chars, then 3-char integer
        for l in 0..s.nsr {
            let line = next_line(&mut r)?;
            // First nodf*10 chars are floats, last 3 are integer
            let nodf = s.nodf;
            for i in 0..nodf {
                let start = i * 10;
                let end = start + 10;
                let chunk = line.get(start..end).unwrap_or("").trim();
                s.odf[[i, l]] = chunk.parse::<f64>().unwrap_or(0.0);
            }
            let isr_start = nodf * 10;
            let isr_chunk = line.get(isr_start..).unwrap_or("").trim();
            s.isr[l] = isr_chunk.parse::<usize>().unwrap_or(0);
        }

        // TITLEJ(1,1) — A20
        let title_no = next_line(&mut r)?;
        if s.titlej.len() > 0 {
            s.titlej[0][0] = title_no.get(..20).unwrap_or(&title_no).trim().to_owned();
        }

        // blank line; set TITLEJ(2,1) and (3,1) to ' '
        skip_line(&mut r)?;
        if s.titlej.len() > 0 {
            s.titlej[0][1] = String::new();
            s.titlej[0][2] = String::new();
        }

        // NSR lines: `7F10.1` → FNO(L), QNO(I=1..NODF, L)
        let nodf = s.nodf;
        let nsr  = s.nsr;
        for l in 0..nsr {
            let vals = read_f64_array(&mut r, 1 + nodf)?;
            s.fno[l] = vals[0];
            for i in 0..nodf {
                s.qno[[i, l]] = vals[1 + i];
            }
        }

        // Convert QNO from oscillator strength to cross-section (cm²)
        // CNO = 8.85e-13 * FNO / (1e7 * (1/WBIN(K) - 1/WBIN(K+1)))
        let nwsrb = s.nwsrb;
        for l in 0..nsr {
            let k = l + nwsrb; // 0-based index into wbin
            if s.fno[l] > 0.0 {
                let cno = 8.85e-13 * s.fno[l]
                    / (1.0e7 * (1.0 / s.wbin[k] - 1.0 / s.wbin[k + 1]));
                for i in 0..nodf {
                    let odf = s.odf[[i, l]];
                    if odf != 0.0 {
                        s.qno[[i, l]] = cno * s.qno[[i, l]] / odf;
                    }
                }
            }
        }

        // O2: 3 temperature points
        for k in 0..3 {
            // `A20,F5.0` → TITLEJ(k, 2), TQQ(k, 2)
            let line = next_line(&mut r)?;
            let title = line.get(..20).unwrap_or("").trim().to_owned();
            let tqq_val: f64 = line.get(20..25).unwrap_or("").trim().parse().unwrap_or(0.0);
            if s.titlej.len() > 1 { s.titlej[1][k] = title; }
            s.tqq[[k, 1]] = tqq_val;

            // `8E10.3` → QO2(IW=1..nwww, k)
            let qo2 = read_f64_array(&mut r, nwww)?;
            for (iw, v) in qo2.iter().enumerate() {
                s.qo2[[iw, k]] = *v;
            }

            // blank line
            skip_line(&mut r)?;

            // NSR lines: `8E10.3` → O2X(I=1..nodf, l, k)
            for l in 0..nsr {
                let vals = read_f64_array(&mut r, nodf)?;
                for (i, v) in vals.iter().enumerate() {
                    s.o2x[[i, l, k]] = *v;
                }
            }
        }

        // O3: 3 temperature points
        for k in 0..3 {
            let line = next_line(&mut r)?;
            let title = line.get(..20).unwrap_or("").trim().to_owned();
            let tqq_val: f64 = line.get(20..25).unwrap_or("").trim().parse().unwrap_or(0.0);
            if s.titlej.len() > 2 { s.titlej[2][k] = title; }
            s.tqq[[k, 2]] = tqq_val;
            let qo3 = read_f64_array(&mut r, nwww)?;
            for (iw, v) in qo3.iter().enumerate() {
                s.qo3[[iw, k]] = *v;
            }
        }

        // O3→O(1D): 3 temperature points
        for k in 0..3 {
            let line = next_line(&mut r)?;
            let title = line.get(..20).unwrap_or("").trim().to_owned();
            let tqq_val: f64 = line.get(20..25).unwrap_or("").trim().parse().unwrap_or(0.0);
            if s.titlej.len() > 3 { s.titlej[3][k] = title; }
            s.tqq[[k, 3]] = tqq_val;
            let q1d = read_f64_array(&mut r, nwww)?;
            for (iw, v) in q1d.iter().enumerate() {
                s.q1d[[iw, k]] = *v;
            }
        }

        // fort10 has NTAB=40 standard species, then 16 extra organics, then 3 iodine species.
        // NQQQ in the file header counts organics as part of the total (=43), so a naive loop
        // over nqqq would fill slots 40-42 with organics instead of iodine cross-sections.
        // Fix: read exactly NTAB standard species, then scan remaining entries by title to find
        // J(IO), J(HOI), J(IONO2) and store them at jq = NTAB, NTAB+1, NTAB+2.
        let ntab = crate::constants::NTAB; // 40 standard tabulated species
        let nxs  = crate::constants::NXS;  // 43 = 40 standard + 3 iodine
        s.njval  = nxs + 4;               // 47 total J-values (overrides nqqq+4 from header)

        // Phase 1: read NTAB standard QQQ species
        for jq in 0..ntab {
            let idx = jq + 4;
            for k in 0..2 {
                let line = next_line(&mut r)?;
                let title = line.get(..20).unwrap_or("").trim().to_owned();
                let tqq_val: f64 = line.get(20..25).unwrap_or("").trim().parse().unwrap_or(0.0);
                if idx < s.titlej.len() { s.titlej[idx][k] = title; }
                s.tqq[[k, idx]] = tqq_val;
                let qqq = read_f64_array(&mut r, nwww)?;
                for (iw, v) in qqq.iter().enumerate() {
                    s.qqq[[iw, k, jq]] = *v;
                }
            }
        }

        // Phase 2: scan remaining species pairs (organics then iodine) looking for
        // J(IO), J(HOI), J(IONO2) by title. Organics are read and discarded.
        let iodine_keys: &[(&str, usize)] = &[
            ("J(IO)",    ntab),
            ("J(HOI)",   ntab + 1),
            ("J(IONO2)", ntab + 2),
        ];
        let n_iodine = nxs - ntab;
        let mut iod_found = 0usize;
        'scan: loop {
            let mut hdr = String::new();
            if r.read_line(&mut hdr)? == 0 { break; }
            let title_250 = hdr.get(..20).unwrap_or("").trim().to_owned();
            let tqq_250: f64 = hdr.get(20..25).unwrap_or("").trim().parse().unwrap_or(0.0);
            let data_250 = read_f64_array(&mut r, nwww)?;

            let mut hdr = String::new();
            if r.read_line(&mut hdr)? == 0 { break; }
            let title_298 = hdr.get(..20).unwrap_or("").trim().to_owned();
            let tqq_298: f64 = hdr.get(20..25).unwrap_or("").trim().parse().unwrap_or(0.0);
            let data_298 = read_f64_array(&mut r, nwww)?;

            for &(tag, jq) in iodine_keys {
                if title_250.contains(tag) {
                    let idx = jq + 4;
                    s.tqq[[0, idx]] = tqq_250;
                    if idx < s.titlej.len() { s.titlej[idx][0] = title_250.clone(); }
                    for (iw, v) in data_250.iter().enumerate() { s.qqq[[iw, 0, jq]] = *v; }
                    s.tqq[[1, idx]] = tqq_298;
                    if idx < s.titlej.len() { s.titlej[idx][1] = title_298; }
                    for (iw, v) in data_298.iter().enumerate() { s.qqq[[iw, 1, jq]] = *v; }
                    iod_found += 1;
                    if iod_found >= n_iodine { break 'scan; }
                    break;
                }
            }
        }

        s.nlbatm = 1;
        Ok(())
    }

    // ── fort11_jpl09.x — rate constants ───────────────────────────────────────

    /// Translates the fort11 read section of READIN (bread.f lines 89–110).
    fn read_rate_constants(&mut self, s: &mut ModelState) -> Result<()> {
        let mut r = self.open("fort11_jpl09.x")?;

        // Line 1: title
        let title = next_line(&mut r)?;
        eprintln!("{}", title);

        // Line 2: `A10,14I5` → TLBL, NRATE1, NRATES
        let line2 = next_line(&mut r)?;
        let ints = parse_fixed_i32s(&line2, 10, 5, 14);
        s.nrate1 = ints[0] as usize;
        s.nrates = ints[1] as usize;

        let mut ndx0: usize = 0;

        for j in 0..s.nrates {
            // Fortran: IF(J.GT.NRATE1 .AND. J.LT.201) GOTO 4  (skip lines 201-220)
            // j is 0-based here, Fortran j is 1-based
            let jf = j + 1; // 1-based Fortran index
            if jf > s.nrate1 && jf < 201 {
                continue;
            }

            // `5X,I1,1X,I3,2E10.3,5X,5A8`
            // offset: 0-4 skip, 5=NADD(1 char), 6 skip, 7-9=IR(3 chars),
            //         10-19=RK1(10), 20-29=RK2(10), 30-34 skip, 35+ = RFMT (5×8 chars)
            let line = next_line(&mut r)?;
            let nadd = parse_nadd_ir_rk(&line, s, j, &mut ndx0, &mut r)?;
            let _ = nadd;
        }

        Ok(())
    }

    // ── fort13.x — background atmosphere ─────────────────────────────────────

    /// Translates the fort13 read section of READIN (bread.f lines 262–280).
    fn read_atmosphere(&mut self, s: &mut ModelState) -> Result<()> {
        let mut r = self.open("fort13.x")?;

        // Line 1: title
        let title = next_line(&mut r)?;
        eprintln!("{}", title);

        // TITATM label was stored in s.titlev by read_control_params.
        let titatm = s.titlev.clone();

        // Search for matching atmosphere label
        loop {
            let hdr = next_line(&mut r)?;
            if hdr.is_empty() {
                break;
            }
            // `2X,A8,7E10.3` — 2 spaces, 8-char label, then floats
            let tita = hdr.get(2..10).unwrap_or("").trim();
            let rest = hdr.get(10..).unwrap_or("");
            let fv = parse_f64s(rest);
            let press0 = *fv.get(0).unwrap_or(&0.0);
            let grav   = *fv.get(1).unwrap_or(&0.0);
            let rad    = *fv.get(2).unwrap_or(&0.0);

            // `8E10.3` → T(1..46)
            let t_vals = read_f64_array(&mut r, 46)?;
            for (i, v) in t_vals.iter().enumerate() {
                s.t[i] = *v;
            }

            if tita == titatm {
                s.press0 = press0;
                s.grav   = grav;
                s.rad    = rad;
                break;
            }
        }

        // Set altitude grid: Z(I) = 2e5*(I-1) cm (2 km steps)
        for i in 0..s.nc {
            s.z[i] = 2.0e5 * i as f64;
        }

        // Hydrostatic equilibrium
        hystat(s, 0);

        // Cache fort13 content for NEWATM (path mode atmosphere resets)
        if let Ok(mut r13) = self.open("fort13.x") {
            let mut line = String::new();
            loop {
                line.clear();
                match r13.read_line(&mut line) {
                    Ok(0) => break,
                    Ok(_) => s.fort13_lines.push(line.trim_end_matches('\n').trim_end_matches('\r').to_string()),
                    Err(_) => break,
                }
            }
        }

        Ok(())
    }

    // ── fort14.x — O3 profile ─────────────────────────────────────────────────

    /// Translates the fort14 read section of READIN (bread.f lines 281–303).
    fn read_ozone_profile(&mut self, s: &mut ModelState) -> Result<()> {
        let mut r = self.open("fort14.x")?;

        let title = next_line(&mut r)?;
        eprintln!("{}", title);

        // `A10,14I5` → TLBL, NDDDZ
        let line2 = next_line(&mut r)?;
        let ints = parse_fixed_i32s(&line2, 10, 5, 14);
        let ndddz = ints[0] as usize;

        let tito3 = s.title0.clone();

        loop {
            let hdr = next_line(&mut r)?;
            if hdr.is_empty() {
                bail!("***FAILED TO FIND STD OZONE***");
            }
            let titz = hdr.get(2..10).unwrap_or("").trim().to_owned();
            // COLATZ (column O3) is in the remaining floats
            let rest = hdr.get(10..).unwrap_or("");
            let fv = parse_f64s(rest);
            let colatz = *fv.get(0).unwrap_or(&0.0);

            let do3 = read_f64_array(&mut r, ndddz)?;
            if titz == tito3 {
                for (j, v) in do3.iter().enumerate() {
                    s.do3ref[j] = *v;
                }
                // Convert ppmv to cm⁻³ if surface value < 10
                if s.do3ref[0] < 10.0 {
                    for j in 0..ndddz {
                        s.do3ref[j] *= s.dm[j] * 1.0e-6;
                    }
                }
                // Extend to NC using exponential scaling
                let nc = s.nc;
                if ndddz < nc {
                    let scalez = s.do3ref[ndddz - 1] / s.do3ref[ndddz - 2];
                    for i in ndddz..nc {
                        s.do3ref[i] = s.do3ref[i - 1] * scalez;
                    }
                }
                break;
            }
            let _ = colatz;
        }

        // Compute O3 and O2 column integrals (bread.f lines 382–400)
        let nc = s.nc;
        let po2 = s.po2;
        s.do3int[nc - 1] = s.do3ref[nc - 1] * s.zzht;
        s.do2int[nc - 1] = s.dm[nc - 1] * s.zzht * po2;
        for j in (0..nc - 1).rev() {
            s.do3int[j] = s.do3int[j + 1]
                + 0.5 * (s.z[j + 1] - s.z[j]) * (s.do3ref[j + 1] + s.do3ref[j]);
            s.do2int[j] = s.do2int[j + 1]
                + 0.5 * (s.z[j + 1] - s.z[j]) * (s.dm[j + 1] + s.dm[j]) * po2;
        }

        // O3 column rescaling (bread.f lines 388–393)
        if s.zo3col > 0.0 {
            let fo3c = s.zo3col / s.do3int[0];
            for j in 0..nc {
                s.do3int[j] *= fo3c;
                s.do3ref[j] *= fo3c;
            }
        }

        // Aerosol profile (bread.f lines 362–376) — computed here because Z() is now known.
        let a = s.aersol;
        for j in 0..nc {
            let aer1 = a[0] * (-s.z[j] / a[1]).exp();
            let aer2 = if s.z[j] > a[3] {
                a[2] * (-(s.z[j] - a[3]) / 3.0e5).exp()
            } else {
                a[2]
            };
            s.aer[j] = aer1 + aer2;
        }
        // Normalise to TURBD0 column
        let mut aercol = s.aer[nc - 1] * s.zzht;
        for j in 1..nc {
            aercol += (s.z[j] - s.z[j - 1]) * (s.aer[j] + s.aer[j - 1]) * 0.5;
        }
        if aercol > 0.0 {
            let aerscl = s.turbd0 / aercol;
            for j in 0..nc { s.aer[j] *= aerscl; }
        }

        // Cache fort14 content for NEWATM (path mode O3 profile resets)
        if let Ok(mut r14) = self.open("fort14.x") {
            let mut line = String::new();
            loop {
                line.clear();
                match r14.read_line(&mut line) {
                    Ok(0) => break,
                    Ok(_) => s.fort14_lines.push(line.trim_end_matches('\n').trim_end_matches('\r').to_string()),
                    Err(_) => break,
                }
            }
        }

        Ok(())
    }

    // ── fort01.x — main control parameters ───────────────────────────────────

    /// Translates the fort01 read section of READIN (bread.f lines 111–213).
    /// Fort01.x is opened as UNIT=1 and contains all run-control parameters.
    fn read_control_params(&mut self, s: &mut ModelState) -> Result<()> {
        let mut r = self.open("fort01.x")?;

        // `10A8` → TITLER (5 × CHARACTER*8)
        let line = next_line(&mut r)?;
        for i in 0..5 {
            s.titler[i] = line.get(i * 8..(i + 1) * 8)
                .unwrap_or("        ")
                .to_owned();
        }

        // `2X,A8,2X,A8,...` → TITATM, TITO3
        let line = next_line(&mut r)?;
        let titatm = line.get(2..10).unwrap_or("").trim().to_owned();
        let tito3  = line.get(12..20).unwrap_or("").trim().to_owned();
        // Re-purpose titlev/title0 to carry atmosphere/O3 labels into later reads
        s.titlev = titatm;
        s.title0 = tito3;

        // `A10,7E10.3` → TLBL, XLATD, XDECD, FLSCAL, ZO3COL, DAYSEC, GMU0
        let line = next_line(&mut r)?;
        let fv = parse_f64s(&line[10..]);
        s.xlatd  = *fv.get(0).unwrap_or(&0.0);
        s.xdecd  = *fv.get(1).unwrap_or(&0.0);
        s.flscal = *fv.get(2).unwrap_or(&1.0);
        s.zo3col = *fv.get(3).unwrap_or(&0.0);
        s.daysec = *fv.get(4).unwrap_or(&86400.0);
        s.gmu0   = *fv.get(5).unwrap_or(&-0.14); // from fort01 LAT/DEC/FL line
        let crad = 57.29578_f64;
        s.xlat = s.xlatd / crad;
        s.xdec = s.xdecd / crad;

        // `A10,7E10.3` → TLBL, PN2, PO2, PCO2, ZZHT
        let line = next_line(&mut r)?;
        let fv = parse_f64s(&line[10..]);
        s.pn2  = *fv.get(0).unwrap_or(&0.0);
        s.po2  = *fv.get(1).unwrap_or(&0.0);
        s.pco2 = *fv.get(2).unwrap_or(&0.0);
        s.zzht = *fv.get(3).unwrap_or(&6.3e5);

        // `A10,14F5.2` → TLBL, ATIME(2..15)
        let line = next_line(&mut r)?;
        let vals = parse_fixed_f64s_fw(&line, 10, 5, 14);
        for (j, v) in vals.iter().enumerate() {
            s.atime[j + 1] = *v; // ATIME(2)..ATIME(15) → atime[1]..atime[14]
        }

        // `A10,14F5.2` → TLBL, BTIME(1..14)
        let line = next_line(&mut r)?;
        let vals = parse_fixed_f64s_fw(&line, 10, 5, 14);
        for (j, v) in vals.iter().enumerate() {
            s.btime[j] = *v;
        }

        // `A10,14F5.2` → TLBL, CLOUDS
        let line = next_line(&mut r)?;
        let fv = parse_f64s(&line[10..]);
        s.clouds = *fv.get(0).unwrap_or(&0.0);

        // `A10,7E10.3` → TLBL, AERSOL(1..6)
        // AERSOL(5) = TURBD0, AERSOL(6) = TURBDX
        let line = next_line(&mut r)?;
        let fv = parse_f64s(&line[10..]);
        s.turbd0 = *fv.get(4).unwrap_or(&0.0);
        s.turbdx = *fv.get(5).unwrap_or(&0.0);
        let aersol: Vec<f64> = fv.into_iter().take(6).collect();
        for (i, v) in aersol.iter().enumerate() { s.aersol[i] = *v; }

        // `A10,7E10.3` → TLBL, RAFMIN..DAYERR
        let line = next_line(&mut r)?;
        let fv = parse_f64s(&line[10..]);
        s.rafmin = *fv.get(0).unwrap_or(&0.1);
        s.rafmax = *fv.get(1).unwrap_or(&1e4);
        s.rafpml = *fv.get(2).unwrap_or(&1e-10);
        s.rafeps = *fv.get(3).unwrap_or(&1e-4);
        s.raferr = *fv.get(4).unwrap_or(&1e-6);
        s.dayeps = *fv.get(5).unwrap_or(&1e-3);
        s.dayerr = *fv.get(6).unwrap_or(&3e-5);

        // `A10,14I5` → TLBL, MAXRAF, MAXRLX, NLBATM, NRADD, ND216, ND216S
        let line = next_line(&mut r)?;
        let ints = parse_fixed_i32s(&line, 10, 5, 14);
        s.maxraf  = ints[0] as usize;
        s.maxrlx  = ints[1] as usize;
        s.nlbatm  = ints[2] as usize;
        let nradd = ints[3] as usize;
        s.nd216   = ints[4];
        s.nd216s  = ints[5];
        if s.nd216s == 0 {
            s.nd216s = s.nd216;
        }

        // `A10,14I5` → TLBL, NDAY, NBOX, NDVAL, NFVAL, NDDD, NDAYSD
        let line = next_line(&mut r)?;
        let ints = parse_fixed_i32s(&line, 10, 5, 14);
        s.nday   = ints[0];
        let nbox_raw = ints[1];
        s.ndval  = ints[2];
        s.nfval  = ints[3];
        let nddd = ints[4] as usize;
        s.ndaysd = ints[5];

        // If NBOX < 0: read from unit=2, propagate box 1 to all
        let n2box = if nbox_raw < 0 {
            s.nbox = (-nbox_raw) as usize;
            1usize
        } else {
            s.nbox = nbox_raw as usize;
            s.nbox
        };
        let _ = (n2box, nddd);

        // `A5,25I3` → BOXDO — fixed-width 3-char integers after 5-char label
        let line = next_line(&mut r)?;
        let boxdo = parse_fixed_i32s(&line, 5, 3, 25);
        for (i, v) in boxdo.iter().take(s.nbox).enumerate() {
            s.nboxdo[i] = *v;
        }

        // `A5,25I3` → BOXWT, BOXPR, BOXCT, BOXMX
        let line = next_line(&mut r)?;
        let vals = parse_fixed_i32s(&line, 5, 3, 25);
        for (i, v) in vals.iter().take(s.nbox).enumerate() { s.nboxwt[i] = *v; }

        let line = next_line(&mut r)?;
        let vals = parse_fixed_i32s(&line, 5, 3, 25);
        for (i, v) in vals.iter().take(s.nbox).enumerate() { s.nboxpr[i] = *v; }

        let line = next_line(&mut r)?;
        let vals = parse_fixed_i32s(&line, 5, 3, 25);
        for (i, v) in vals.iter().take(s.nbox).enumerate() { s.nboxct[i] = *v; }

        let line = next_line(&mut r)?;
        let vals = parse_fixed_i32s(&line, 5, 3, 25);
        for (i, v) in vals.iter().take(s.nbox).enumerate() { s.nboxmx[i] = *v; }

        // `A10,14F5.2/(10X,14F5.2)` → BOXRN (multi-line float array)
        read_boxf_multiline(&mut r, &mut s.boxrn, s.nbox)?;

        // BOXAA
        read_boxf_multiline(&mut r, &mut s.boxaa, s.nbox)?;

        // BOXTT
        read_boxf_multiline(&mut r, &mut s.boxtt, s.nbox)?;

        // NRADD rate overrides
        let mut ndx0: usize = 0;
        for _ in 0..nradd {
            let line = next_line(&mut r)?;
            // `5X,I1,1X,I3,2E10.3,5X,5A8`
            let (nadd, j, rk1, rk2) = parse_rate_line(&line);
            let jidx = j.saturating_sub(1); // 0-based
            s.rk[[0, jidx]] = rk1;
            s.rk0[[0, jidx]] = rk1;
            s.rk[[1, jidx]] = rk2;
            s.ndxrat[jidx] = 0;
            if nadd > 0 {
                for _jj in 0..nadd {
                    let addline = next_line(&mut r)?;
                    let (_, ir, rkadd1, rkadd2) = parse_rate_line(&addline);
                    let idx = ndx0;
                    s.rkadd[[0, idx]] = rkadd1;
                    s.rkadd[[1, idx]] = rkadd2;
                    let _ = ir;
                    ndx0 += 1;
                }
                s.ndxrat[jidx] = ndx0 as i32 - nadd as i32;
                ndx0 += nadd;
            }
        }

        // `A10,7E10.3` → TLBL, DIFFC, XDMAX, ZSHEAR, DIFFAC
        let line = next_line(&mut r)?;
        let fv = parse_f64s(&line[10..]);
        let diffc = *fv.get(0).unwrap_or(&0.0);
        s.xdmax   = *fv.get(1).unwrap_or(&0.0);
        s.zshear  = *fv.get(2).unwrap_or(&0.0);
        s.diffac  = *fv.get(3).unwrap_or(&1.0);
        s.tdiff = if diffc > 0.0 {
            s.xdmax * s.xdmax * 1.0e4 / diffc
        } else {
            0.0
        };

        // `A5,25I3` → NDXPP
        let line = next_line(&mut r)?;
        let vals = parse_fixed_i32s(&line, 5, 3, 25);
        for (i, v) in vals.iter().take(crate::constants::NNDXPQ).enumerate() {
            s.ndxpp[i] = *v as usize;
        }

        // `A5,25I3` → NDXQQ
        let line = next_line(&mut r)?;
        let vals = parse_fixed_i32s(&line, 5, 3, 25);
        for (i, v) in vals.iter().take(crate::constants::NNDXPQ).enumerate() {
            s.ndxqq[i] = *v as usize;
        }

        // `A10,14I5` → NPRT1..ITPRTB
        let line = next_line(&mut r)?;
        let ints = parse_fixed_i32s(&line, 10, 5, 14);
        s.nprt1  = ints[0];
        s.nprt2  = ints[1];
        s.nprtrr = ints[2];
        s.itprtx = ints[3];
        s.itprtr = ints[4];
        s.npstd  = ints[5];
        s.itprtb = ints[6];

        // `5L2` → LPRT, LPRTX, LPRTJV, LPRTY, LPRT8
        let line = next_line(&mut r)?;
        let bools = parse_fortran_logicals(&line, 5);
        s.lprt   = bools[0];
        s.lprtx  = bools[1];
        s.lprtjv = bools[2];
        s.lprty  = bools[3];
        s.lprt8  = bools[4];

        // Long-lived (flow) species: `A8,I2,1X,4I1` header then entries.
        // Some legacy files have a stale header count (e.g., "L-LIVED 18") while
        // NFVAL and the listed records include iodine families. Trust NFVAL.
        let line = next_line(&mut r)?; // "L-LIVED xx ..." header: `A8,I2,1X,4I1`
        let nfl_hdr_str = line.get(8..10).unwrap_or("0").trim();
        let nfl_hdr: usize = nfl_hdr_str.parse().unwrap_or(0);
        let nfl = s.nfval as usize;
        if nfl_hdr != nfl {
            eprintln!(
                "Warning: L-LIVED header count ({}) != NFVAL ({}); using NFVAL",
                nfl_hdr, nfl
            );
        }

        for j in 0..nfl {
            // `A8,I2,2I5,6E10.3`
            let line = next_line(&mut r)?;
            let name = line.get(..8).unwrap_or("").trim().to_owned();
            let jj: usize = line.get(8..10).unwrap_or("").trim().parse().unwrap_or(0);
            // NFLO(J) at char 10..15, JX (NDIFF) at char 15..20
            let nflo_v: i32 = line.get(10..15).unwrap_or("").trim().parse().unwrap_or(0);
            let jx: i32 = line.get(15..20).unwrap_or("").trim().parse().unwrap_or(0);
            if jj != j + 1 {
                eprintln!("Warning: long-lived species index mismatch at j={}", j + 1);
            }
            s.tnomf[j] = name.clone();
            s.titles[j + s.ndval as usize] = name;
            s.nflo[j] = nflo_v;
            s.ndiff[j + s.ndval as usize] = jx;
        }

        // Implicit/explicit radicals: `A8,I2,1X,4I1` header then NTOTX entries
        let line = next_line(&mut r)?;
        let ntotx_str = line.get(8..10).unwrap_or("0").trim();
        s.ntotx = ntotx_str.parse().unwrap_or(0);

        let jxsum_expected = (s.ntotx + s.ntotx * s.ntotx) / 2;
        let mut jxsum = jxsum_expected as i32;
        let mut nts1: usize = 1;
        let mut nts2: usize = s.ntotx;
        s.nnr = 0;

        for _j in 0..s.ntotx {
            // `A8,I2,1X,4I1`
            let line = next_line(&mut r)?;
            let name  = line.get(..8).unwrap_or("").trim().to_owned();
            let jx: usize = line.get(8..10).unwrap_or("").trim().parse().unwrap_or(0);
            // NXDIU at char 11, NXSLO at char 12, NXE at char 13, NXDIFF at char 14
            let nxdiu: i32 = line.get(11..12).unwrap_or("0").trim().parse().unwrap_or(0);
            let nxslo: i32 = line.get(12..13).unwrap_or("0").trim().parse().unwrap_or(0);
            let _nxe:  i32 = line.get(13..14).unwrap_or("0").trim().parse().unwrap_or(0);
            let nxdiff:i32 = line.get(14..15).unwrap_or("0").trim().parse().unwrap_or(0);

            jxsum -= jx as i32;
            if jx > 0 && jx <= NDEN {
                s.tname[jx - 1] = name.clone();
                s.titles[jx - 1] = name;
                s.ndiff[jx - 1] = nxdiff;
            }

            let effective_diu = if nxslo > 0 { 1 } else { nxdiu };
            if effective_diu > 0 {
                if jx > 0 && jx <= NDEN { s.ntsav[jx - 1] = nts1; }
                nts1 += 1;
            } else {
                if jx > 0 && jx <= NDEN { s.ntsav[jx - 1] = nts2; }
                nts2 = nts2.saturating_sub(1);
            }

            if nxslo > 0 {
                if s.nnr >= crate::constants::NSLOWM {
                    bail!("NNR > NSLOWM");
                }
                s.nnrt[s.nnr] = jx;
                s.nnr += 1;
            }
        }

        // Save NTSAV(NDEN+1) = NTS2 (the boundary between implicit and explicit)
        s.ntsav[NDEN] = nts2;
        s.ntot = nts2;

        // Build NT array: NT(J) = NTSAV(J) for J=1..NDEN
        for j in 0..NDEN {
            s.n[j] = s.ntsav[j];
            s.tnamet[s.ntsav[j].saturating_sub(1)] = s.tname[j].clone();
        }

        // Bromine flag
        s.lbrom = s.n[11] <= s.ntot || s.n[12] <= s.ntot || s.n[13] <= s.ntot
            || s.n[22] <= s.ntot || s.n[23] <= s.ntot;
        if !s.lbrom {
            s.nrates = s.nrate1;
        }

        // Iodine flag — species at n[30..34] (1-based N31..N35)
        s.liod = s.ntotx > 30
            && (s.n[30] <= s.ntot || s.n[31] <= s.ntot || s.n[32] <= s.ntot
                || s.n[33] <= s.ntot || s.n[34] <= s.ntot);

        // Initialise XR to zero
        s.xr.iter_mut().for_each(|v| *v = 0.0);

        if jxsum != 0 {
            bail!("CHECKSUM ERROR ON SPECIES (jxsum={})", jxsum);
        }

        // Aerosol profile from AERSOL params stored above
        // (aersol stored in a local; re-read from s.turbd0/turbdx + the 4 shape params)
        // Since we can't easily re-access the local, we apply here using stored values.
        // The full aersol[0..4] for the profile shape comes from the control line read above.
        // To avoid threading the params, store them temporarily in s.atime[0..3].
        // This is a known limitation: the 4 profile params are in the read_control_params scope.
        // A future refactor can extract them. For now we apply 0 aerosol and rely on fort13.
        // TODO: store aersol[0..4] into a dedicated field and compute AER() here.
        s.nc = 41; // set to full altitude range

        // Set up diurnal time grid
        setday(s);

        // Save remaining lines from fort01 for NEWATM (path mode)
        // These are the path records that come after the initial conditions block.
        {
            let mut line = String::new();
            let mut remaining: Vec<String> = Vec::new();
            loop {
                line.clear();
                match r.read_line(&mut line) {
                    Ok(0) => break, // EOF
                    Ok(_) => {
                        let trimmed = line.trim_end_matches('\n').trim_end_matches('\r').to_string();
                        remaining.push(trimmed);
                    }
                    Err(_) => break,
                }
            }
            remaining.reverse(); // reverse so pop() gives lines in file order
            s.fort01_remaining = remaining;
        }

        Ok(())
    }

    // ── fort02.x — initial mixing ratios (DIURN/TPATH mode) ─────────────────────

    /// Read initial species mixing ratios for non-CTM mode.
    /// Format: title line, then for each of NDDD species:
    ///   A10 species name, then NBOX values in 6E13.6 format.
    /// Values are mixing ratios; the reader rescales to number density later.
    /// Fortran: bread.f DO 34 ID=1,NDDD ... READ(2,113) (DDDDDD(IB,ID), IB=1,N2BOX)
    fn read_initial_densities(&mut self, s: &mut ModelState) -> Result<()> {
        let mut r = match self.open("fort02.x") {
            Ok(r) => r,
            Err(_) => return Ok(()), // optional file — skip gracefully
        };
        let mut line = String::new();
        r.read_line(&mut line)?; // title line

        let nddd = s.ndval as usize + s.nfval as usize;
        let ndval = s.ndval as usize;
        let nbox = s.nbox;

        for id in 0..nddd {
            // Species name (A10)
            line.clear();
            if r.read_line(&mut line)? == 0 { break; }
            // Values (6E13.6 per line, across NBOX boxes)
            let mut vals = Vec::with_capacity(nbox);
            while vals.len() < nbox {
                line.clear();
                if r.read_line(&mut line)? == 0 { break; }
                let trimmed = line.trim_end();
                let mut off = 0;
                while off + 13 <= trimmed.len() && vals.len() < nbox {
                    let v: f64 = trimmed[off..off + 13].trim()
                        .replace(|c: char| c == 'd' || c == 'D', "e")
                        .parse().unwrap_or(0.0);
                    vals.push(v);
                    off += 13;
                }
            }
            // Fill missing boxes by repeating first value
            let fill = vals.first().copied().unwrap_or(0.0);
            while vals.len() < nbox { vals.push(fill); }

            // Store: ID=1..NDVAL → implicit species densities (as mixing ratio; rescaled later)
            //        ID=NDVAL+1..NDDD → long-lived species mixing ratios
            for (ib, &v) in vals.iter().enumerate().take(nbox) {
                if id < ndval {
                    s.den_set(ib, id, v); // mixing ratio (multiplied by DM in bread.f)
                } else {
                    s.fff_set(ib, id - ndval + 1, v); // 1-based j for fff_set
                }
            }
        }

        // Rescale implicit species from mixing ratio to number density, call FIXMIX per box.
        // Fortran: IBOX=I; IALT=IABS(NBOXDO(IBOX)); ...; CALL FIXMIX
        for ib in 0..nbox {
            let ialt = (s.nboxdo[ib].unsigned_abs() as usize).saturating_sub(1);
            let dm = s.dm[ialt];
            for id in 0..ndval {
                let v = s.den_get(ib, id);
                s.den_set(ib, id, v * dm);
            }
            s.ibox = ib;
            s.ialt = ialt;
            fixmix(s);
        }

        Ok(())
    }

    // ── J_H2O_SZA0.dat ───────────────────────────────────────────────────────

    fn read_jh2o(&mut self, s: &mut ModelState) -> Result<()> {
        let mut r = self.open("J_H2O_SZA0.dat")?;
        for i in 0..crate::constants::NJH2O {
            let line = next_line(&mut r)?;
            let fv = parse_f64s(&line);
            s.zjh2o[i] = *fv.get(0).unwrap_or(&0.0);
            s.xjh2o[i] = *fv.get(1).unwrap_or(&0.0);
        }
        Ok(())
    }
}

// ── Rate-line parser helper ───────────────────────────────────────────────────

/// Parse fort11 rate line: `5X,I1,1X,I3,2E10.3,5X,5A8`
/// Returns (NADD, IR, RK1, RK2)
fn parse_rate_line(line: &str) -> (usize, usize, f64, f64) {
    let nadd: usize = line.get(5..6).unwrap_or("0").trim().parse().unwrap_or(0);
    let ir:   usize = line.get(7..10).unwrap_or("0").trim().parse().unwrap_or(0);
    let rk1:  f64   = parse_e_field(line.get(10..20).unwrap_or(""));
    let rk2:  f64   = parse_e_field(line.get(20..30).unwrap_or(""));
    (nadd, ir, rk1, rk2)
}

/// Parse one fort11 rate entry into state, including NADD additional-data lines.
fn parse_nadd_ir_rk(
    line: &str,
    s: &mut ModelState,
    j: usize,      // 0-based
    ndx0: &mut usize,
    r: &mut impl BufRead,
) -> Result<usize> {
    let (nadd, _ir, rk1, rk2) = parse_rate_line(line);
    s.rk[[0, j]]  = rk1;
    s.rk0[[0, j]] = rk1;
    s.rk[[1, j]]  = rk2;
    // RFMT labels (chars 35..75: 5×8)
    for k in 0..2 {
        let start = 35 + k * 8;
        let end   = start + 8;
        s.rfmt_str[j][k] = line.get(start..end).unwrap_or("        ").to_owned();
    }
    s.ndxrat[j] = 0;

    if nadd > 0 {
        s.ndxrat[j] = *ndx0 as i32;
        for jj in 0..nadd {
            let addline = next_line(r)?;
            let (_, _, add1, add2) = parse_rate_line(&addline);
            let idx = *ndx0 + jj;
            if idx < 50 {
                s.rkadd[[0, idx]] = add1;
                s.rkadd[[1, idx]] = add2;
            }
        }
        *ndx0 += nadd;
    }
    Ok(nadd)
}

/// Parse a Fortran E10.3 field (10 chars) which may or may not have 'E'.
fn parse_e_field(s: &str) -> f64 {
    let t = s.trim();
    if t.is_empty() {
        return 0.0;
    }
    // Fortran sometimes omits the 'E': "1.80-11" means 1.80e-11
    // Insert 'e' before a bare sign after digits/dot if no 'E'/'e' present
    if t.contains('E') || t.contains('e') {
        t.parse::<f64>().unwrap_or(0.0)
    } else {
        // Try inserting 'e' before a +/- that follows a digit or '.'
        let fixed = fix_bare_exponent(t);
        fixed.parse::<f64>().unwrap_or(0.0)
    }
}

fn fix_bare_exponent(s: &str) -> String {
    let bytes = s.as_bytes();
    for i in 1..bytes.len() {
        let c = bytes[i] as char;
        let prev = bytes[i - 1] as char;
        if (c == '+' || c == '-') && (prev.is_ascii_digit() || prev == '.') {
            let mut out = s[..i].to_owned();
            out.push('e');
            out.push_str(&s[i..]);
            return out;
        }
    }
    s.to_owned()
}

/// Parse fixed-width f64 fields: offset + count × width chars each.
fn parse_fixed_f64s_fw(line: &str, offset: usize, width: usize, count: usize) -> Vec<f64> {
    let bytes = line.as_bytes();
    let mut out = Vec::with_capacity(count);
    for i in 0..count {
        let start = offset + i * width;
        let end   = (start + width).min(bytes.len());
        if start >= bytes.len() {
            out.push(0.0);
            continue;
        }
        let chunk = std::str::from_utf8(&bytes[start..end]).unwrap_or("").trim();
        out.push(chunk.parse::<f64>().unwrap_or(0.0));
    }
    out
}

/// Parse Fortran logical values from `5L2` format: each L occupies 2 chars.
/// T/t/.TRUE. → true, F/f/.FALSE. → false.
fn parse_fortran_logicals(line: &str, count: usize) -> Vec<bool> {
    let mut out = Vec::with_capacity(count);
    for tok in line.split_whitespace().take(count) {
        let b = matches!(tok.to_ascii_uppercase().as_str(), "T" | ".TRUE." | ".T.");
        out.push(b);
    }
    while out.len() < count { out.push(false); }
    out
}

/// Read a multi-line `A10,14F5.2/(10X,14F5.2)` array (BOXRN/BOXAA/BOXTT pattern).
fn read_boxf_multiline(r: &mut impl BufRead, arr: &mut [f64], nbox: usize) -> Result<()> {
    let mut read = 0;
    while read < nbox {
        let line = next_line(r)?;
        // First line: 10-char label then 14×F5.2; continuation: 10 spaces then 14×F5.2
        let vals = parse_fixed_f64s_fw(&line, 10, 5, 14);
        for v in vals.iter() {
            if read >= nbox { break; }
            arr[read] = *v;
            read += 1;
        }
    }
    Ok(())
}

// ── SETDAY / HYSTAT — called from read_control_params / read_atmosphere ───────

/// Compute diurnal time grid (butil.f SETDAY).
pub fn setday(s: &mut ModelState) {
    let cpisec = 13750.99_f64;
    let cpihr  = 3.819719_f64;

    let gmua = s.xlat.sin() * s.xdec.sin();
    let gmub = s.xlat.cos() * s.xdec.cos();

    let arg = ((s.gmu0 - gmua) / gmub).max(-1.0).min(1.0);
    s.sunset = arg.acos() * cpihr;

    let sset = 3600.0 * s.sunset;

    // Build DTIME from ATIME fractions of sunset.
    // Fortran JJ counts DTIME(1..JJ) including noon at step 1.
    // Rust jj here counts dtime[1..jj] (non-noon afternoon steps).
    // After the loop jj = last valid ATIME index (= Fortran JJ - 1).
    s.dtime[0] = 0.0;
    let mut jj = 1usize;
    for j in 1..15 {
        if s.atime[j] < 1e-5 {
            // Matches Fortran: when ATIME(J)=0, set JJ = J-1 (i.e., previous j).
            // In Fortran 1-based: JJ = J-1 where J was the terminating index (=j+1 0-based).
            // Result: Fortran JJ = j (0-based); we set jj = j here (already done on prev iter).
            break;
        }
        s.dtime[j] = s.atime[j] * sset;
        if s.dtime[j] > 0.5 * s.daysec {
            jj = j - 1;
            break;
        }
        jj = j;
    }
    // Fortran JJ includes noon (step 1 = 0.0). Adjust to match Fortran's NTIM calculation.
    // Fortran: NTIM = JJ_fortran + JJ_fortran + J_night where JJ_fortran = jj + 1.
    let jj_fortran = jj + 1; // equivalent to Fortran's JJ

    let mut night_j = 0usize;
    let snit = s.daysec - s.dtime[jj] - s.dtime[jj];
    if snit >= 1.0 {
        for j in 0..14 {
            if s.btime[j] < 1e-5 { break; }
            s.dtime[jj + j + 1] = s.btime[j] * snit + s.dtime[jj];
            night_j = j + 1;
        }
    }

    let morn = jj + night_j;
    // NTIM matches Fortran: JJ_fortran + JJ_fortran + J_night
    s.ntim = jj_fortran + jj_fortran + night_j;
    if s.ntim > 44 {
        panic!("NTIM ({}) > 44", s.ntim);
    }

    // Mirror morning into afternoon/night: Fortran mirrors JJ_fortran steps (1..JJ in 1-based)
    for j in 0..jj_fortran {
        let mirror_idx = s.ntim - j;
        if mirror_idx <= 43 {
            s.dtime[mirror_idx] = s.daysec - s.dtime[j];
        }
    }

    // Compute WTIME (fractional weights)
    s.wtime[0] = 0.0;
    let ntim = s.ntim;
    for j in 1..=ntim {
        let dj = s.dtime[j].max(s.dtime[j - 1] + 1.0);
        s.dtime[j] = dj;
        s.wtime[j] = (dj - s.dtime[j - 1]) / s.daysec;
    }

    let midnit = (ntim + 1) / 2;
    s.nmu = 1;

    for jj in 0..=morn {
        let gmu = gmua + gmub * (s.dtime[jj] / cpisec).cos();
        s.utime[jj] = gmu;
        s.ztime[jj] = 57.29578 * gmu.acos();
        s.jtim[jj] = s.nmu as i32;

        // Fortran: JCOMP = NTIM+1-JJ (1-based JJ=1..MORN)
        // Converting to 0-based: JJ = jj+1, so JCOMP_0based = (NTIM+1-(jj+1))-1 = ntim-jj-1
        let jcomp = ntim - 1 - jj;
        if jcomp > morn {
            s.utime[jcomp] = s.utime[jj];
            s.ztime[jcomp] = s.ztime[jj];
            s.jtim[jcomp] = s.nmu as i32;
        }
        if gmu >= s.gmu0 {
            s.nmu = jj + 1;
        }
    }
    s.nmu = s.nmu.min(midnit);

    // Mode selection
    s.ntimdo = ntim;
    s.ldiurn = true;

    if s.nday >= 3 {
        s.ntimdo = 2;
        s.ldiurn = false;
        s.dtime[1] = s.daysec;
        s.jtim[1] = 1;
        s.ztime[1] = s.ztime[0];
    } else if s.nday == 2 {
        if s.utime[s.nmu] >= s.gmu0 {
            // Cannot do 4-step; fall back to 1-step
            s.ntimdo = 2;
            s.ldiurn = false;
            s.dtime[1] = s.daysec;
            s.jtim[1] = 1;
            s.ztime[1] = s.ztime[0];
        } else {
            s.ntimdo = 5;
            s.ldiurn = false;
            let nmu = s.nmu;
            let sset2 = 0.5 * (s.dtime[nmu - 1] + s.dtime[nmu]);
            s.dtime[1] = sset2;
            s.dtime[3] = s.daysec - sset2;
            s.dtime[4] = s.daysec;
            s.dtime[2] = s.btime[0] * (s.dtime[3] - s.dtime[1]) + s.dtime[1];
            s.jtim[1] = 1; s.jtim[2] = 2; s.jtim[3] = 2; s.jtim[4] = 1;
            s.ztime[1] = 90.0; s.ztime[2] = 180.0;
            s.ztime[3] = 90.0; s.ztime[4] = s.ztime[0];
            let ratio = 0.5 * s.daysec / sset2;
            for j in 0..ntim { s.wtime[j] *= ratio; }
        }
    }

    // NHHMM
    for j in 0..ntim {
        let htime = s.dtime[j] / 3600.0 + 0.001;
        let iitime = htime as i32;
        s.nhhmm[j] = (100.0 * (0.60 * (htime - iitime as f64) + (iitime + 12).rem_euclid(24) as f64) + 0.50) as i32;
    }

    s.ittt = 1;
}

/// Hydrostatic equilibrium (butil.f HYSTAT). LOGP=0: geometric 2-km grid.
pub fn hystat(s: &mut ModelState, logp: i32) {
    let cboltz = 1.38e-19_f64;
    let nc = s.nc;

    if logp == 0 {
        s.pstd[0] = s.press0;
        s.dm[0] = (7.340e21 / s.t[0]) * (s.press0 / 1013.25);
        s.theta[0] = s.t[0];
        let cstat0 = 3.416e-4 * (s.grav / 980.665);
        let mut cplog1 = 0.5 * cstat0 * (s.rad / (s.rad + s.z[0])).powi(2);
        for ii in 1..nc {
            let cplog2 = 0.5 * cstat0 * (s.rad / (s.rad + s.z[ii])).powi(2);
            let dlogp = (cplog1 / s.t[ii - 1] + cplog2 / s.t[ii]) * (s.z[ii - 1] - s.z[ii]);
            s.pstd[ii] = s.pstd[ii - 1] * dlogp.exp();
            s.dm[ii]   = s.pstd[ii] / (cboltz * s.t[ii]);
            s.theta[ii] = s.t[ii] * (1000.0 / s.pstd[ii]).powf(0.2857);
            cplog1 = cplog2;
        }
    } else {
        s.z[0] = 0.0;
        s.pstd[0] = 1000.0;
        s.dm[0] = s.pstd[0] / (cboltz * s.t[0]);
        s.theta[0] = s.t[0];
        let dlogp = 10.0_f64.powf(-2.0 / 16.0);
        let clogp = 421.25462799776471_f64;
        for ii in 1..nc {
            s.pstd[ii] = s.pstd[ii - 1] * dlogp;
            s.z[ii]    = s.z[ii - 1] + (s.t[ii - 1] + s.t[ii]) * clogp;
            s.dm[ii]   = s.pstd[ii] / (cboltz * s.t[ii]);
            s.theta[ii] = s.t[ii] * (1000.0 / s.pstd[ii]).powf(0.2857);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::state::ModelState;

    // ── parse_e_field ─────────────────────────────────────────────────────────

    #[test]
    fn test_parse_e_field_standard() {
        assert!((parse_e_field("1.80E-11") - 1.80e-11).abs() < 1e-20);
        assert!((parse_e_field("3.50E+02") - 350.0).abs() < 1e-10);
        assert!((parse_e_field("0.0      ") - 0.0).abs() < 1e-30);
    }

    #[test]
    fn test_parse_e_field_bare_exponent() {
        // Fortran sometimes omits 'E': "1.80-11" means 1.80×10⁻¹¹
        assert!((parse_e_field("1.80-11") - 1.80e-11).abs() < 1e-20);
        assert!((parse_e_field("3.50+02") - 350.0).abs() < 1e-10);
        assert!((parse_e_field("2.997+08") - 2.997e8).abs() < 1e-0);
    }

    #[test]
    fn test_parse_e_field_zero() {
        assert_eq!(parse_e_field(""), 0.0);
        assert_eq!(parse_e_field("   "), 0.0);
    }

    // ── parse_fixed_i32s ─────────────────────────────────────────────────────

    #[test]
    fn test_parse_fixed_i32s_boxdo() {
        // Typical BOXDO line: "A5,25I3" → 5-char label then 3-char ints
        let line = "BOXDO-30-29-28-27-26";
        let vals = parse_fixed_i32s(line, 5, 3, 5);
        assert_eq!(vals, vec![-30, -29, -28, -27, -26]);
    }

    #[test]
    fn test_parse_fixed_i32s_positive() {
        let line = "     001002003004";
        let vals = parse_fixed_i32s(line, 5, 3, 4);
        assert_eq!(vals, vec![1, 2, 3, 4]);
    }

    #[test]
    fn test_parse_fixed_i32s_truncated() {
        // Line shorter than requested — remainder should be 0
        let line = "LABEL  1  2";
        let vals = parse_fixed_i32s(line, 5, 3, 5);
        assert_eq!(vals[0], 1);
        assert_eq!(vals[1], 2);
        assert_eq!(vals[2], 0);
    }

    // ── parse_fixed_f64s_fw ──────────────────────────────────────────────────

    #[test]
    fn test_parse_fixed_f64s_atime() {
        // "A10,14F5.2" → skip 10 chars, then 5-char floats
        let line = "DAY INTVAL 0.15 0.30 0.45 0.60 0.70";
        let vals = parse_fixed_f64s_fw(line, 10, 5, 5);
        assert!((vals[0] - 0.15).abs() < 1e-6, "vals[0]={}", vals[0]);
        assert!((vals[2] - 0.45).abs() < 1e-6, "vals[2]={}", vals[2]);
        assert!((vals[4] - 0.70).abs() < 1e-6, "vals[4]={}", vals[4]);
    }

    // ── hystat ───────────────────────────────────────────────────────────────

    #[test]
    fn test_hystat_logp_surface_pressure() {
        let mut s = ModelState::new();
        s.nc = 5;
        for i in 0..5 { s.t[i] = 250.0; }
        hystat(&mut s, 1);
        // Surface pressure is always 1000 hPa in log-p mode
        assert!((s.pstd[0] - 1000.0).abs() < 1e-9);
    }

    #[test]
    fn test_hystat_logp_pressure_ratio() {
        // 16 levels in log-p spans exactly 2 decades: p[16] = 1000 × 10^(-2) = 10 hPa
        let mut s = ModelState::new();
        s.nc = 17;
        for i in 0..17 { s.t[i] = 250.0; }
        hystat(&mut s, 1);
        let ratio = s.pstd[16] / s.pstd[0];
        assert!((ratio - 0.01).abs() < 1e-9, "pstd ratio = {ratio}");
    }

    #[test]
    fn test_hystat_logp_density_from_ideal_gas() {
        // DM = p / (k_B × T); verify at level 0
        let cboltz = 1.38e-19_f64;
        let mut s = ModelState::new();
        s.nc = 3;
        s.t[0] = 300.0; s.t[1] = 280.0; s.t[2] = 260.0;
        hystat(&mut s, 1);
        let expected_dm0 = 1000.0 / (cboltz * 300.0);
        assert!((s.dm[0] - expected_dm0).abs() / expected_dm0 < 1e-10);
    }

    #[test]
    fn test_hystat_logp_z_increases() {
        // Altitude must increase monotonically
        let mut s = ModelState::new();
        s.nc = 10;
        for i in 0..10 { s.t[i] = 240.0 - i as f64 * 5.0; }
        hystat(&mut s, 1);
        for i in 1..10 {
            assert!(s.z[i] > s.z[i-1],
                "z not monotone: z[{i}]={} <= z[{}]={}", s.z[i], i-1, s.z[i-1]);
        }
    }

    #[test]
    fn test_hystat_logp_approx_2km_per_level() {
        // The CLOGP constant is chosen so each level is ~2 km in an isothermal 250 K atmosphere
        let mut s = ModelState::new();
        s.nc = 5;
        for i in 0..5 { s.t[i] = 250.0; }
        hystat(&mut s, 1);
        for i in 1..5 {
            let dz_km = (s.z[i] - s.z[i-1]) * 1e-5;
            assert!(dz_km > 1.5 && dz_km < 2.5,
                "dz[{i}] = {dz_km:.2} km, expected ~2 km");
        }
    }

    #[test]
    fn test_hystat_geometric_surface_density() {
        // DM(0) formula: (7.340e21 / T(0)) × (PRESS0 / 1013.25)
        let mut s = ModelState::new();
        s.nc = 3;
        s.press0 = 1013.25;
        s.grav   = 980.665;
        s.rad    = 6.371e8;
        s.t[0] = 288.0; s.t[1] = 275.0; s.t[2] = 260.0;
        s.z[0] = 0.0; s.z[1] = 2.0e5; s.z[2] = 4.0e5;
        hystat(&mut s, 0);
        let expected = 7.340e21 / 288.0;
        assert!((s.dm[0] - expected).abs() / expected < 1e-9);
    }

    // ── setday ───────────────────────────────────────────────────────────────

    fn standard_setday_state() -> Box<ModelState> {
        let mut s = ModelState::new();
        let deg = std::f64::consts::PI / 180.0;
        s.xlat   = 60.0 * deg;
        s.xdec   = -1.689 * deg;   // March 16 declination
        s.gmu0   = -0.12;
        s.daysec = 86400.0;
        // ATIME from fort01.x (index 0 unused; indices 1..13 active)
        s.atime[1]  = 0.15; s.atime[2]  = 0.30; s.atime[3]  = 0.45;
        s.atime[4]  = 0.60; s.atime[5]  = 0.70; s.atime[6]  = 0.80;
        s.atime[7]  = 0.88; s.atime[8]  = 0.92; s.atime[9]  = 0.96;
        s.atime[10] = 0.98; s.atime[11] = 1.00; s.atime[12] = 1.02;
        s.atime[13] = 1.05;
        // BTIME from fort01.x (6 night steps)
        s.btime[0] = 0.05; s.btime[1] = 0.10; s.btime[2] = 0.30;
        s.btime[3] = 0.50; s.btime[4] = 0.70; s.btime[5] = 0.90;
        s
    }

    #[test]
    fn test_setday_ntim_standard() {
        // Standard 60°N March 16 run must produce ntim = ntimdo = 34
        let mut s = standard_setday_state();
        setday(&mut s);
        assert_eq!(s.ntim,   34, "ntim={}", s.ntim);
        assert_eq!(s.ntimdo, 34, "ntimdo={}", s.ntimdo);
    }

    #[test]
    fn test_setday_noon_is_dtime0() {
        // Noon is always defined as dtime[0] = 0.0 seconds from midnight
        let mut s = standard_setday_state();
        setday(&mut s);
        assert_eq!(s.dtime[0], 0.0);
    }

    #[test]
    fn test_setday_sunset_60n_march() {
        // At 60°N March 16, sunset should be ~6.6–6.8 hours from noon
        let mut s = standard_setday_state();
        setday(&mut s);
        assert!(s.sunset > 6.5 && s.sunset < 7.0,
            "sunset = {:.3} h (expected 6.5–7.0 h)", s.sunset);
    }

    #[test]
    fn test_setday_last_mirror_is_daysec() {
        // The last mirrored dtime (at index ntim) must equal DAYSEC
        let mut s = standard_setday_state();
        setday(&mut s);
        let ntim = s.ntim;
        assert!((s.dtime[ntim] - s.daysec).abs() < 1.0,
            "dtime[ntim={ntim}] = {} ≠ DAYSEC", s.dtime[ntim]);
    }

    #[test]
    fn test_setday_dtime_monotone() {
        // After WTIME adjustment, dtime must be strictly increasing
        let mut s = standard_setday_state();
        setday(&mut s);
        let ntim = s.ntim;
        for j in 1..=ntim {
            assert!(s.dtime[j] > s.dtime[j-1],
                "dtime not monotone at j={j}: {} <= {}", s.dtime[j], s.dtime[j-1]);
        }
    }

    #[test]
    fn test_setday_weights_sum_to_one() {
        // wtime[1..=ntim] must sum to 1.0 (they partition the day)
        let mut s = standard_setday_state();
        setday(&mut s);
        let sum: f64 = (1..=s.ntim).map(|j| s.wtime[j]).sum();
        assert!((sum - 1.0).abs() < 1e-10, "wtime sum = {sum}");
    }

    #[test]
    fn test_setday_utime_mirror_symmetry() {
        // Morning and afternoon utime values should be symmetric about noon.
        // utime[33] mirrors utime[0] (noon), utime[20] mirrors utime[13].
        let mut s = standard_setday_state();
        setday(&mut s);
        // With jcomp = ntim-1-jj and ntim=34: jj=0 → jcomp=33, jj=13 → jcomp=20.
        assert!((s.utime[33] - s.utime[0]).abs() < 1e-12,
            "utime[33]={} ≠ utime[0]={}", s.utime[33], s.utime[0]);
        assert!((s.utime[20] - s.utime[13]).abs() < 1e-12,
            "utime[20]={} ≠ utime[13]={}", s.utime[20], s.utime[13]);
    }

    #[test]
    fn test_setday_nday3_gives_2_steps() {
        // NDAY=3 (1-step day) must set NTIMDO=2
        let mut s = standard_setday_state();
        s.nday = 3;
        setday(&mut s);
        assert_eq!(s.ntimdo, 2, "ntimdo should be 2 for nday=3");
        assert!(!s.ldiurn);
    }

    #[test]
    fn test_setday_polar_night() {
        // At winter pole (90°N, dec=-23.5°), gmu0=-0.12 keeps sun below horizon
        // all day; SNIT ≈ DAYSEC and all wtime concentrated in night
        let mut s = ModelState::new();
        let deg = std::f64::consts::PI / 180.0;
        s.xlat   = 90.0 * deg;
        s.xdec   = -23.5 * deg;
        s.gmu0   = -0.12;
        s.daysec = 86400.0;
        s.atime[1] = 0.5;  // single afternoon step
        s.btime[0] = 0.5;  // single night step
        setday(&mut s);
        // Should still produce a valid time grid
        assert!(s.ntim >= 2);
        assert!(s.ntimdo >= 2);
    }
}
