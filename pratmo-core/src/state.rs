use ndarray::{Array2, Array3, Array4};

use crate::constants::*;

/// All COMMON blocks from bcomm.h unified into one heap-allocated struct.
/// Use `Box::new(ModelState::default())` — the struct is ~3 MB.
///
/// Array index conventions:
///   - All indices are 0-based (Fortran 1-based → subtract 1)
///   - Multi-dim arrays follow Fortran column-major logical ordering where noted
pub struct ModelState {
    // ── COMMON / TITLES / ───────────────────────────────────────────────────
    pub title0: String,                         // CHARACTER*78
    pub titlev: String,                         // CHARACTER*78
    pub titlej: Vec<[String; 3]>,               // CHARACTER*20 TITLEJ(3, NXS+4) → [NJVAL][3]

    // ── COMMON/CCATMO/ ──────────────────────────────────────────────────────
    /// Temperature profile (K); T(46) in Fortran — 46 to allow interpolation headroom
    pub t: [f64; 46],
    pub refo3: [f64; NL],       // reference O3 profile (cm⁻³)
    pub pstd: [f64; NL],        // standard pressure profile (mb)
    pub dm: [f64; NL],          // air number density (cm⁻³)
    pub do3ref: [f64; NL],      // local O3 reference density (cm⁻³)
    pub z: [f64; NL],           // altitude (cm)
    pub theta: [f64; NL],       // potential temperature (K)
    pub aer: [f64; NL],         // aerosol extinction (km⁻¹)
    pub do3int: [f64; NL],      // integrated O3 column above (cm⁻²·atm)
    /// WTAU(NL,NL): optical depth table [from_level][to_level]
    pub wtau: Array2<f64>,
    pub dnoy_ref: [f64; NL],    // reference NOy profile (cm⁻³)
    pub rflect: f64,            // surface albedo
    pub sza: f64,               // solar zenith angle (degrees)
    pub u0: f64,                // cos(SZA)
    pub tanht: f64,             // tangent height (km)
    pub nlbatm: usize,          // number of atmospheric levels used
    pub nc: usize,              // number of active altitude levels
    pub ncdim: usize,           // declared dimension for nc
    pub nbdim: usize,           // declared dimension for nb
    pub nbox: usize,            // number of boxes in current run
    pub nd216: i32,             // run mode flag (>0: CTM, 0: DIURN/TPATH, <0: DERIVS)
    pub nd216s: i32,            // saved nd216

    // ── COMMON/CCETC/ ───────────────────────────────────────────────────────
    pub press0: f64,            // surface pressure (mb)
    pub pn2: f64,               // N2 partial pressure
    pub po2: f64,               // O2 partial pressure
    pub pco2: f64,              // CO2 partial pressure
    pub rad: f64,               // Earth radius (cm)
    pub grav: f64,              // gravity (cm/s²)
    pub zzht: f64,              // scale height (cm)
    pub xlat: f64,              // latitude (radians)
    pub xdec: f64,              // solar declination (radians)
    pub xlatd: f64,             // latitude (degrees)
    pub xdecd: f64,             // solar declination (degrees)
    pub atime: [f64; 15],       // time array A
    pub btime: [f64; 15],       // time array B
    pub clouds: f64,            // cloud fraction
    pub turbd0: f64,            // turbidity parameter 0
    pub turbdx: f64,            // turbidity parameter x
    pub flscal: f64,            // solar flux scale factor

    // ── COMMON/CCSCAT/ ──────────────────────────────────────────────────────
    pub tau: [f64; 42],         // total optical depth per layer
    pub piray: [f64; 42],       // Rayleigh scattering phase integral
    pub piaer: [f64; 42],       // aerosol scattering phase integral
    pub fltau: [f64; 42],       // filtered tau
    pub ttfbot: f64,            // total tau at bottom
    pub ntt: usize,             // number of tau levels

    // ── COMMON/CCWORK/ ──────────────────────────────────────────────────────
    /// Previous time step species densities for 30 implicit species
    pub xnold: [f64; NDEN],
    pub deltt: f64,             // current time step (s)
    pub gmu: f64,               // cos(SZA) at current time step
    pub ittt: i32,              // iteration counter
    pub ialt: usize,            // current altitude level index (0-based)
    pub ibox: usize,            // current box index (0-based)

    // ── COMMON/CCDEN/ ───────────────────────────────────────────────────────
    // Species number densities (cm⁻³), one per box.
    // Fortran: DNO(NB) ... DBRCL(NB) laid out consecutively.
    // Equivalence: DDDDDD(NB, NSPEC) overlaps starting at DNO(1).
    pub dno: [f64; NB],
    pub dno2: [f64; NB],
    pub dno3: [f64; NB],
    pub dn2o5: [f64; NB],
    pub dhno3: [f64; NB],
    pub dh: [f64; NB],
    pub doh: [f64; NB],
    pub dho2: [f64; NB],
    pub dh2o2: [f64; NB],
    pub do_: [f64; NB],         // DO is a Rust keyword, renamed do_
    pub do3: [f64; NB],
    pub dbro: [f64; NB],
    pub dbr: [f64; NB],
    pub dhbr: [f64; NB],
    pub dhno2: [f64; NB],
    pub dhcl: [f64; NB],
    pub dcl: [f64; NB],
    pub dcl2: [f64; NB],
    pub dclo: [f64; NB],
    pub dclno3: [f64; NB],
    pub dhno4: [f64; NB],
    pub dhocl: [f64; NB],
    pub dbrno3: [f64; NB],
    pub dhobr: [f64; NB],
    pub dh2co: [f64; NB],
    pub droo: [f64; NB],
    pub drooh: [f64; NB],
    pub doclo: [f64; NB],
    pub dcl2o2: [f64; NB],
    pub dbrcl: [f64; NB],
    // Long-lived species mixing ratios (dimensionless) — FFFFFF(NB,18)
    pub fo3: [f64; NB],
    pub fn2o: [f64; NB],
    pub fnoy: [f64; NB],
    pub fch4: [f64; NB],
    pub fco: [f64; NB],
    pub fclx: [f64; NB],
    pub fcf2cl: [f64; NB],
    pub fcfcl3: [f64; NB],
    pub fccl4: [f64; NB],
    pub fch3cl: [f64; NB],
    pub fmecl: [f64; NB],       // CH3CCl3
    pub fh2: [f64; NB],
    pub fh2o: [f64; NB],
    pub fnh3: [f64; NB],
    pub fc5h8: [f64; NB],
    pub fbrx: [f64; NB],
    pub fch3br: [f64; NB],
    pub focs: [f64; NB],
    pub fxxx: [f64; NB],
    /// Species name labels (CHARACTER*8)
    pub titler: [String; 5],
    pub titles: Vec<String>,    // TITLES(NSPEC+3)
    pub ndval: i32,
    pub nfval: i32,

    // ── COMMON/CCWVL/ ───────────────────────────────────────────────────────
    /// Wavelength bin boundaries (nm); WBIN(NQ+1)
    pub wbin: [f64; NQ + 1],
    /// Wavelength bin centres (nm); WL(NQ)
    pub wl: [f64; NQ],
    /// Solar flux at each bin; FL(NQ)
    pub fl: [f64; NQ],
    /// O2 absorption cross-sections; QO2(NQ, 3) — 3 temperature points
    pub qo2: Array2<f64>,       // shape [NQ, 3]
    /// O2 Schumann-Runge band data; O2X(6, 15, 3)
    pub o2x: Array3<f64>,       // shape [6, 15, 3]
    /// O2 SRB oscillator strengths; ODF(6, 15)
    pub odf: Array2<f64>,       // shape [6, 15]
    /// NO photolysis flux coefficients; FNO(15)
    pub fno: [f64; 15],
    /// NO cross-sections; QNO(6, 15)
    pub qno: Array2<f64>,       // shape [6, 15]
    /// O3 cross-sections; QO3(NQ, 3)
    pub qo3: Array2<f64>,       // shape [NQ, 3]
    /// O3→O(1D) cross-sections; Q1D(NQ, 3)
    pub q1d: Array2<f64>,       // shape [NQ, 3]
    /// Additional species cross-sections; QQQ(NQ, 2, NXS) — 2 temperature points
    pub qqq: Array3<f64>,       // shape [NQ, 2, NXS]
    /// Temperature interpolation points for each J-value; TQQ(3, NXS+4)
    pub tqq: Array2<f64>,       // shape [3, NJVAL]
    pub qmi1m: f64,
    pub qmiwvl: f64,
    pub nw1: usize,             // first wavelength bin for fast-J
    pub nw2: usize,             // last wavelength bin for fast-J
    pub nwsrb: usize,           // number of SRB bins
    pub nsr: usize,             // number of SRB sub-intervals
    pub nodf: usize,            // number of ODF points
    pub isr: [usize; 15],       // SRB bin indices

    // ── COMMON/CCJVL/ ───────────────────────────────────────────────────────
    // J-value profiles (s⁻¹) at each altitude level.
    // Equivalence: VVVVVV(NL, NXS+4) overlaps starting at VNO(1).
    pub vno: [f64; NL],
    pub vo2: [f64; NL],
    pub vo3: [f64; NL],
    pub vo3d: [f64; NL],        // O3 → O(1D)
    pub vh2coa: [f64; NL],
    pub vh2cob: [f64; NL],
    pub vh2o2: [f64; NL],
    pub vrooh: [f64; NL],
    pub vno2: [f64; NL],
    pub vno3x: [f64; NL],
    pub vno3l: [f64; NL],
    pub vn2o5: [f64; NL],
    pub vhno2: [f64; NL],
    pub vhno3: [f64; NL],
    pub vhno4: [f64; NL],
    pub vclno3: [f64; NL],
    pub vcl2: [f64; NL],
    pub vhocl: [f64; NL],
    pub voclo: [f64; NL],
    pub vcl2o2: [f64; NL],
    pub vclo: [f64; NL],
    pub vbro: [f64; NL],
    pub vbrno3: [f64; NL],
    pub vhobr: [f64; NL],
    pub vn2o: [f64; NL],
    pub vcfcl3: [f64; NL],
    pub vf2cl2: [f64; NL],
    pub vf113: [f64; NL],
    pub vf114: [f64; NL],
    pub vf115: [f64; NL],
    pub vccl4: [f64; NL],
    pub vch3cl: [f64; NL],
    pub vmecf: [f64; NL],
    pub vch3br: [f64; NL],
    pub vh1211: [f64; NL],
    pub vh1301: [f64; NL],
    pub vh2402: [f64; NL],
    pub vh22: [f64; NL],
    pub vh123: [f64; NL],
    pub vh141b: [f64; NL],
    pub vchbr3: [f64; NL],
    pub vch3i: [f64; NL],
    pub vcf3i: [f64; NL],
    pub vocs: [f64; NL],
    /// Actinic flux per wavelength bin at each level; FFF(NQ, NL)
    pub fff: Array2<f64>,       // shape [NQ, NL]
    pub njval: usize,           // number of J-values computed

    // ── COMMON/CCNNN/ ───────────────────────────────────────────────────────
    // Species index mapping: N1..N30 are the 0-based indices for each of the
    // 30 implicit species in the Newton-Raphson solve (NTOT = N30).
    pub n: [usize; 30],         // N1..N30 mapped to n[0]..n[29]
    pub ntot: usize,
    pub ntotx: usize,
    pub nnrt: [usize; NSLOWM],
    pub nnr: usize,
    pub ntsav: [usize; 31],
    pub nflo: [i32; 20],        // explicit integration flags per species
    pub tname: [String; 30],    // species names for Newton-Raphson set
    pub tnamet: [String; 30],   // transient species names
    pub tnomf: [String; 20],    // long-lived species names

    // ── COMMON/CCRTS/ ───────────────────────────────────────────────────────
    // Equivalence: RCOLUM(430) = XR(30) ++ R(250) ++ RP(30) ++ RL(30)
    //              ++ RPF(30) ++ RLF(30) ++ RQF(30)  (total = 430)
    pub xr: [f64; 30],          // species densities passed to solver
    pub r: [f64; NR],           // reaction rates (cm⁻³ s⁻¹ or s⁻¹)
    pub rp: [f64; 30],          // production rates
    pub rl: [f64; 30],          // loss rates
    pub rpf: [f64; 30],         // production rates (family)
    pub rlf: [f64; 30],         // loss rates (family)
    pub rqf: [f64; 30],         // quasi-steady-state rates (family)

    // ── COMMON/DMEANS/ ──────────────────────────────────────────────────────
    /// Daily mean diagnostic storage; PMEAN(490)
    pub pmean: [f64; NPMEAN],
    /// Species densities at each time step; XNOFT(30, 44)
    pub xnoft: Array2<f64>,     // shape [NDEN, NXNOFT]
    /// Per-box mean diagnostics; PPMEAN(50+NNDXPQ, NB)
    pub ppmean: Array2<f64>,    // shape [50+NNDXPQ, NB]
    /// Per-box species at each time step; XXNOFT(30, 44, NB)
    pub xxnoft: Array3<f64>,    // shape [NDEN, NXNOFT, NB]
    /// Stored J-values; STORJV(NXS+4, 16, NB)
    pub storjv: Array3<f64>,    // shape [NJVAL, 16, NB]
    pub ndxpp: [usize; NNDXPQ],
    pub ndxqq: [usize; NNDXPQ],

    // ── COMMON/CCKINE/ ──────────────────────────────────────────────────────
    /// Arrhenius parameters; RK(2, 250) — (A, Ea) or similar
    pub rk: Array2<f64>,        // shape [2, NR]
    /// Rate format codes; RFMT(2, 250) stored as integers encoded in f64
    pub rfmt: Array2<f64>,      // shape [2, NR]  (CHARACTER*8 in Fortran → encoded)
    pub rfmt_str: Vec<[String; 2]>, // actual CHARACTER*8 values; RFMT(2,250)
    /// Additional rate parameters; RKADD(2, 50)
    pub rkadd: Array2<f64>,     // shape [2, 50]
    /// Rate index mapping; NDXRAT(250)
    pub ndxrat: [i32; NR],
    pub nrates: usize,          // total number of rates
    pub nrate1: usize,          // index of first photolysis rate
    /// Zero-pressure rate parameters; rk0(2, 250)
    pub rk0: Array2<f64>,       // shape [2, NR]

    // ── COMMON/CCRATE/ ──────────────────────────────────────────────────────
    /// Evaluated rate constants at current T, M; RATEK(250)
    pub ratek: [f64; NR],
    pub ztemp: f64,             // current temperature (K)
    pub zdnum: f64,             // current number density (cm⁻³)
    pub zalt: f64,              // current altitude (cm)
    pub izalt: usize,           // current altitude level index

    // ── COMMON/CCSOLV/ ──────────────────────────────────────────────────────
    /// Jacobian matrix for Newton-Raphson; A(30, 30)
    pub a_mat: Array2<f64>,     // shape [NDEN, NDEN]
    /// Pivot indices from LU; IPA(30)
    pub ipa: [usize; NDEN],
    pub astore: [f64; NDEN],
    pub aextra: [f64; NDEN],

    // ── COMMON/CDAILY/ ──────────────────────────────────────────────────────
    /// Diurnal time array (s from midnight); DTIME(44)
    pub dtime: [f64; NXNOFT],
    /// Time weights; WTIME(44)
    pub wtime: [f64; NXNOFT],
    /// UTC time array; UTIME(44)
    pub utime: [f64; NXNOFT],
    /// Zenith angle time array; ZTIME(44)
    pub ztime: [f64; NXNOFT],
    pub gmu0: f64,              // noon cos(SZA)
    pub daysec: f64,            // seconds per day
    pub sunset: f64,            // sunset time (s from midnight)
    pub jtim: [i32; NXNOFT],    // time step type flags
    pub nhhmm: [i32; NXNOFT],   // time labels (HHMM format)
    pub ntim: usize,            // number of diurnal time steps
    pub nmu: usize,             // number of cos(SZA) points
    pub ntimdo: usize,          // actual number of time steps used

    // ── COMMON/CCRAL/ ───────────────────────────────────────────────────────
    pub rafmin: f64,            // min time step fraction
    pub rafmax: f64,            // max time step fraction
    pub rafpml: f64,            // P-L convergence tolerance
    pub rafeps: f64,            // species convergence epsilon
    pub raferr: f64,            // species convergence error limit
    pub dayeps: f64,            // daily convergence epsilon
    pub dayerr: f64,            // daily convergence error
    pub maxraf: usize,          // max NR iterations
    pub maxrlx: usize,          // max relaxation iterations

    // ── COMMON/CCPARM/ ──────────────────────────────────────────────────────
    pub boxrn: [f64; NB],       // box run number / identifier
    pub boxaa: [f64; NB],       // box surface albedo
    pub boxtt: [f64; NB],       // box temperature offset (K)
    pub nboxpr: [i32; NB],      // print flags per box
    pub nboxdo: [i32; NB],      // altitude level for each box (1-based; negative = use DAILY)
    pub nboxmx: [i32; NB],      // max integration days per box
    pub nboxwt: [i32; NB],      // weight flag per box
    pub nboxct: [i32; NB],      // count flag per box
    pub nprt1: i32,
    pub nprt2: i32,
    pub nprtrr: i32,
    pub itprtx: i32,
    pub itprtr: i32,
    pub itprtb: i32,
    pub nday: i32,              // diurnal flag (0=noon avg, 1=full 24h)
    pub ndaysd: i32,            // number of integration days
    pub npstd: i32,

    // ── COMMON/CCDIFF/ ──────────────────────────────────────────────────────
    pub tdiff: f64,
    pub xdmax: f64,
    pub zshear: f64,
    pub diffac: f64,
    pub ndiff: [i32; NSPEC],

    // ── COMMON/CCL/ ─────────────────────────────────────────────────────────
    pub lprt: bool,
    pub lprtx: bool,
    pub lprts: bool,
    pub ldiurn: bool,
    pub ljzer: bool,
    pub lbrom: bool,
    pub lsvjac: bool,
    pub lsvday: bool,
    pub lend: bool,
    pub lprtjv: bool,
    pub lresol: bool,
    pub lprty: bool,
    pub lprt8: bool,

    // ── COMMON/CHRIS/ ───────────────────────────────────────────────────────
    /// CTM diagnostic storage; dpl(12, 18, NB, NTAB)
    pub dpl: Array4<f64>,       // shape [12, 18, NB, NTAB]
    pub asa: [f64; NL],         // aerosol surface area density (μm²/cm³)
    /// Density output; denout(12, 18, 3, 4, NB)
    pub denout: ndarray::Array5<f64>, // shape [12, 18, 3, 4, NB]
    /// AM/PM ratio means; rampm(12, 18, 9, NB)
    pub rampm: Array4<f64>,     // shape [12, 18, 9, NB]
    pub aersf: f64,             // aerosol scale factor
    pub smxsf: [f64; NQ],       // spectral solar flux scale factors
    pub ssf: [f64; NQ],         // spectral scale factors
    pub coclo: [f64; 44],       // ClO diurnal cycle
    pub cbro: [f64; 44],        // BrO diurnal cycle
    pub cno2: [f64; 44],        // NO2 diurnal cycle
    pub cno3: [f64; 44],        // NO3 diurnal cycle
    pub cclo: [f64; 44],        // Cl diurnal cycle
    pub chono: [f64; 44],       // HONOx diurnal cycle
    pub xpo3: f64,              // O3 production
    pub xlo3: f64,              // O3 loss
    pub do2int: [f64; NL],      // integrated O2 column above (cm⁻²·atm)
    /// O3 tracer profiles; o3tr(12, 55)
    pub o3tr: Array2<f64>,      // shape [12, 55]
    pub ztr: [f64; 55],         // tracer altitude grid (cm)
    /// Water albedo table; watalb(11, 34)
    pub watalb: Array2<f64>,    // shape [11, 34]
    pub wvalb: [f64; 34],       // wavelength grid for albedo table
    pub szalb: [f64; 11],       // SZA grid for albedo table
    /// SZA-dependent albedo; ztalb(44, 34)
    pub ztalb: Array2<f64>,     // shape [44, 34]
    pub dn2oref: [f64; NL],     // reference N2O profile (cm⁻³)

    // ── COMMON/CHRIS2/ ──────────────────────────────────────────────────────
    /// Input temperature climatology; TINP(12, 18, 46)
    pub tinp: Array3<f64>,      // shape [12, 18, 46]
    /// Input O3 climatology; DO3INP(12, 18, NL)
    pub do3inp: Array3<f64>,    // shape [12, 18, NL]
    /// Input N2O climatology; DN2OINP(12, 18, NL)
    pub dn2oinp: Array3<f64>,   // shape [12, 18, NL]
    /// Input NOy climatology; DNOyINP(12, 18, NL)
    pub dnoyi_np: Array3<f64>,  // shape [12, 18, NL]
    pub zin: [f64; 101],        // altitude grid for tracer input (cm)
    /// Tracer mixing ratio input; xin(6, 101)
    pub xin: Array2<f64>,       // shape [6, 101]
    pub pin_t: [f64; 101],      // pressure grid for T input (mb)
    pub tfull: [f64; 101],      // full T profile
    pub fbrx_ref: [f64; NL],    // reference Bry mixing ratio profile
    pub zjh2o: [f64; NJH2O],    // altitude grid for H2O photolysis (cm)
    pub xjh2o: [f64; NJH2O],    // H2O photolysis rate table (s⁻¹)
    pub xjdo: f64,              // H2O photolysis rate at current box

    // ── COMMON/ICHRIS/ ──────────────────────────────────────────────────────
    pub nwalb: usize,           // number of wavelength points in albedo table
    pub izalb: usize,           // number of SZA points in albedo table
    pub nztr: usize,            // number of tracer altitude levels
    pub iwn2o: usize,           // N2O input weighting flag
    pub iwnoy: usize,           // NOy input weighting flag
    pub iwbry: usize,           // Bry input weighting flag

    // ── COMMON/CCHRIS/ ──────────────────────────────────────────────────────
    pub cinpdir: String,        // CHARACTER*2 input directory prefix

    // ── COMMON/NICKLLOYD/ ───────────────────────────────────────────────────
    pub raxloop: f64,           // diagnostic loop counter
    pub radcount: f64,          // diagnostic iteration counter

    // ── Scratch / cross-read communication ──────────────────────────────────
    /// Aerosol profile shape parameters read from fort01: AERSOL(1..6).
    /// Used to compute AER() after Z() is known from HYSTAT.
    pub aersol: [f64; 6],
    /// Column O3 override from fort01 (ZO3COL); applied after O3 profile read.
    pub zo3col: f64,

    // ── Output file handles (Fortran units 7, 8, 9) ─────────────────────────
    pub out_unit7: Option<std::io::BufWriter<std::fs::File>>,  // PUNCH output
    pub out_unit8: Option<std::io::BufWriter<std::fs::File>>,  // PRTPTH species
    pub out_unit9: Option<std::io::BufWriter<std::fs::File>>,  // PRTPTH rates

    // ── Cached file content for NEWATM ──────────────────────────────────────
    /// Remaining fort01 records after initial READIN (reversed; pop from back = read forward)
    pub fort01_remaining: Vec<String>,
    /// All lines of fort13.x (atmosphere profiles); used by NEWATM
    pub fort13_lines: Vec<String>,
    /// All lines of fort14.x (ozone profiles); used by NEWATM
    pub fort14_lines: Vec<String>,
}

impl ModelState {
    pub fn new() -> Box<Self> {
        Box::new(Self {
            title0: String::new(),
            titlev: String::new(),
            titlej: vec![[String::new(), String::new(), String::new()]; NJVAL],

            t: [0.0; 46],
            refo3: [0.0; NL],
            pstd: [0.0; NL],
            dm: [0.0; NL],
            do3ref: [0.0; NL],
            z: [0.0; NL],
            theta: [0.0; NL],
            aer: [0.0; NL],
            do3int: [0.0; NL],
            wtau: Array2::zeros((NL, NL)),
            dnoy_ref: [0.0; NL],
            rflect: 0.0,
            sza: 0.0,
            u0: 0.0,
            tanht: 0.0,
            nlbatm: 0,
            nc: 0,
            ncdim: NL,
            nbdim: NB,
            nbox: 0,
            nd216: 0,
            nd216s: 0,

            press0: 0.0,
            pn2: 0.0,
            po2: 0.0,
            pco2: 0.0,
            rad: 0.0,
            grav: 0.0,
            zzht: 0.0,
            xlat: 0.0,
            xdec: 0.0,
            xlatd: 0.0,
            xdecd: 0.0,
            atime: [0.0; 15],
            btime: [0.0; 15],
            clouds: 0.0,
            turbd0: 0.0,
            turbdx: 0.0,
            flscal: 0.0,

            tau: [0.0; 42],
            piray: [0.0; 42],
            piaer: [0.0; 42],
            fltau: [0.0; 42],
            ttfbot: 0.0,
            ntt: 0,

            xnold: [0.0; NDEN],
            deltt: 0.0,
            gmu: 0.0,
            ittt: 0,
            ialt: 0,
            ibox: 0,

            dno: [0.0; NB],
            dno2: [0.0; NB],
            dno3: [0.0; NB],
            dn2o5: [0.0; NB],
            dhno3: [0.0; NB],
            dh: [0.0; NB],
            doh: [0.0; NB],
            dho2: [0.0; NB],
            dh2o2: [0.0; NB],
            do_: [0.0; NB],
            do3: [0.0; NB],
            dbro: [0.0; NB],
            dbr: [0.0; NB],
            dhbr: [0.0; NB],
            dhno2: [0.0; NB],
            dhcl: [0.0; NB],
            dcl: [0.0; NB],
            dcl2: [0.0; NB],
            dclo: [0.0; NB],
            dclno3: [0.0; NB],
            dhno4: [0.0; NB],
            dhocl: [0.0; NB],
            dbrno3: [0.0; NB],
            dhobr: [0.0; NB],
            dh2co: [0.0; NB],
            droo: [0.0; NB],
            drooh: [0.0; NB],
            doclo: [0.0; NB],
            dcl2o2: [0.0; NB],
            dbrcl: [0.0; NB],
            fo3: [0.0; NB],
            fn2o: [0.0; NB],
            fnoy: [0.0; NB],
            fch4: [0.0; NB],
            fco: [0.0; NB],
            fclx: [0.0; NB],
            fcf2cl: [0.0; NB],
            fcfcl3: [0.0; NB],
            fccl4: [0.0; NB],
            fch3cl: [0.0; NB],
            fmecl: [0.0; NB],
            fh2: [0.0; NB],
            fh2o: [0.0; NB],
            fnh3: [0.0; NB],
            fc5h8: [0.0; NB],
            fbrx: [0.0; NB],
            fch3br: [0.0; NB],
            focs: [0.0; NB],
            fxxx: [0.0; NB],
            titler: [
                String::new(), String::new(), String::new(),
                String::new(), String::new(),
            ],
            titles: vec![String::new(); NSPEC + 3],
            ndval: 0,
            nfval: 0,

            wbin: [0.0; NQ + 1],
            wl: [0.0; NQ],
            fl: [0.0; NQ],
            qo2: Array2::zeros((NQ, 3)),
            o2x: Array3::zeros((6, 15, 3)),
            odf: Array2::zeros((6, 15)),
            fno: [0.0; 15],
            qno: Array2::zeros((6, 15)),
            qo3: Array2::zeros((NQ, 3)),
            q1d: Array2::zeros((NQ, 3)),
            qqq: Array3::zeros((NQ, 2, NXS)),
            tqq: Array2::zeros((3, NJVAL)),
            qmi1m: 0.0,
            qmiwvl: 0.0,
            nw1: 0,
            nw2: 0,
            nwsrb: 0,
            nsr: 0,
            nodf: 0,
            isr: [0; 15],

            vno: [0.0; NL],
            vo2: [0.0; NL],
            vo3: [0.0; NL],
            vo3d: [0.0; NL],
            vh2coa: [0.0; NL],
            vh2cob: [0.0; NL],
            vh2o2: [0.0; NL],
            vrooh: [0.0; NL],
            vno2: [0.0; NL],
            vno3x: [0.0; NL],
            vno3l: [0.0; NL],
            vn2o5: [0.0; NL],
            vhno2: [0.0; NL],
            vhno3: [0.0; NL],
            vhno4: [0.0; NL],
            vclno3: [0.0; NL],
            vcl2: [0.0; NL],
            vhocl: [0.0; NL],
            voclo: [0.0; NL],
            vcl2o2: [0.0; NL],
            vclo: [0.0; NL],
            vbro: [0.0; NL],
            vbrno3: [0.0; NL],
            vhobr: [0.0; NL],
            vn2o: [0.0; NL],
            vcfcl3: [0.0; NL],
            vf2cl2: [0.0; NL],
            vf113: [0.0; NL],
            vf114: [0.0; NL],
            vf115: [0.0; NL],
            vccl4: [0.0; NL],
            vch3cl: [0.0; NL],
            vmecf: [0.0; NL],
            vch3br: [0.0; NL],
            vh1211: [0.0; NL],
            vh1301: [0.0; NL],
            vh2402: [0.0; NL],
            vh22: [0.0; NL],
            vh123: [0.0; NL],
            vh141b: [0.0; NL],
            vchbr3: [0.0; NL],
            vch3i: [0.0; NL],
            vcf3i: [0.0; NL],
            vocs: [0.0; NL],
            fff: Array2::zeros((NQ, NL)),
            njval: NJVAL,

            n: [0; 30],
            ntot: 0,
            ntotx: 0,
            nnrt: [0; NSLOWM],
            nnr: 0,
            ntsav: [0; 31],
            nflo: [0; 20],
            tname: core::array::from_fn(|_| String::new()),
            tnamet: core::array::from_fn(|_| String::new()),
            tnomf: core::array::from_fn(|_| String::new()),

            xr: [0.0; 30],
            r: [0.0; NR],
            rp: [0.0; 30],
            rl: [0.0; 30],
            rpf: [0.0; 30],
            rlf: [0.0; 30],
            rqf: [0.0; 30],

            pmean: [0.0; NPMEAN],
            xnoft: Array2::zeros((NDEN, NXNOFT)),
            ppmean: Array2::zeros((50 + NNDXPQ, NB)),
            xxnoft: Array3::zeros((NDEN, NXNOFT, NB)),
            storjv: Array3::zeros((NJVAL, 16, NB)),
            ndxpp: [0; NNDXPQ],
            ndxqq: [0; NNDXPQ],

            rk: Array2::zeros((2, NR)),
            rfmt: Array2::zeros((2, NR)),
            rfmt_str: vec![[String::new(), String::new()]; NR],
            rkadd: Array2::zeros((2, 50)),
            ndxrat: [0; NR],
            nrates: 0,
            nrate1: 0,
            rk0: Array2::zeros((2, NR)),

            ratek: [0.0; NR],
            ztemp: 0.0,
            zdnum: 0.0,
            zalt: 0.0,
            izalt: 0,

            a_mat: Array2::zeros((NDEN, NDEN)),
            ipa: [0; NDEN],
            astore: [0.0; NDEN],
            aextra: [0.0; NDEN],

            dtime: [0.0; NXNOFT],
            wtime: [0.0; NXNOFT],
            utime: [0.0; NXNOFT],
            ztime: [0.0; NXNOFT],
            gmu0: 0.0,
            daysec: 0.0,
            sunset: 0.0,
            jtim: [0; NXNOFT],
            nhhmm: [0; NXNOFT],
            ntim: 0,
            nmu: 0,
            ntimdo: 0,

            rafmin: 0.0,
            rafmax: 0.0,
            rafpml: 0.0,
            rafeps: 0.0,
            raferr: 0.0,
            dayeps: 0.0,
            dayerr: 0.0,
            maxraf: 0,
            maxrlx: 0,

            boxrn: [0.0; NB],
            boxaa: [0.0; NB],
            boxtt: [0.0; NB],
            nboxpr: [0; NB],
            nboxdo: [0i32; NB],
            nboxmx: [0; NB],
            nboxwt: [0; NB],
            nboxct: [0; NB],
            nprt1: 0,
            nprt2: 0,
            nprtrr: 0,
            itprtx: 0,
            itprtr: 0,
            itprtb: 0,
            nday: 0,
            ndaysd: 0,
            npstd: 0,

            tdiff: 0.0,
            xdmax: 0.0,
            zshear: 0.0,
            diffac: 0.0,
            ndiff: [0; NSPEC],

            lprt: false,
            lprtx: false,
            lprts: false,
            ldiurn: false,
            ljzer: false,
            lbrom: false,
            lsvjac: false,
            lsvday: false,
            lend: false,
            lprtjv: false,
            lresol: false,
            lprty: false,
            lprt8: false,

            dpl: Array4::zeros((12, 18, NB, NTAB)),
            asa: [0.0; NL],
            denout: ndarray::Array5::zeros((12, 18, 3, 4, NB)),
            rampm: Array4::zeros((12, 18, 9, NB)),
            aersf: 0.0,
            smxsf: [0.0; NQ],
            ssf: [1.0; NQ],
            coclo: [0.0; 44],
            cbro: [0.0; 44],
            cno2: [0.0; 44],
            cno3: [0.0; 44],
            cclo: [0.0; 44],
            chono: [0.0; 44],
            xpo3: 0.0,
            xlo3: 0.0,
            do2int: [0.0; NL],
            o3tr: Array2::zeros((12, 55)),
            ztr: [0.0; 55],
            watalb: Array2::zeros((11, 34)),
            wvalb: [0.0; 34],
            szalb: [0.0; 11],
            ztalb: Array2::zeros((44, 34)),
            dn2oref: [0.0; NL],

            tinp: Array3::zeros((12, 18, 46)),
            do3inp: Array3::zeros((12, 18, NL)),
            dn2oinp: Array3::zeros((12, 18, NL)),
            dnoyi_np: Array3::zeros((12, 18, NL)),
            zin: [0.0; 101],
            xin: Array2::zeros((6, 101)),
            pin_t: [0.0; 101],
            tfull: [0.0; 101],
            fbrx_ref: [0.0; NL],
            zjh2o: [0.0; NJH2O],
            xjh2o: [0.0; NJH2O],
            xjdo: 0.0,

            nwalb: 0,
            izalb: 0,
            nztr: 0,
            iwn2o: 0,
            iwnoy: 0,
            iwbry: 0,

            cinpdir: String::from(".\\"),

            raxloop: 0.0,
            radcount: 0.0,
            aersol: [0.0; 6],
            zo3col: 0.0,

            out_unit7: None,
            out_unit8: None,
            out_unit9: None,
            fort01_remaining: Vec::new(),
            fort13_lines: Vec::new(),
            fort14_lines: Vec::new(),
        })
    }

    // ── Convenience accessors for the CCDEN species arrays ──────────────────
    // These provide a 2D view compatible with Fortran's DDDDDD(IB, ISPEC) access.
    // The 30 implicit species in order (0-based):
    //   0=NO, 1=NO2, 2=NO3, 3=N2O5, 4=HNO3, 5=H, 6=OH, 7=HO2, 8=H2O2,
    //   9=O, 10=O3, 11=BrO, 12=Br, 13=HBr, 14=HNO2, 15=HCl, 16=Cl,
    //   17=Cl2, 18=ClO, 19=ClONO2, 20=HNO4, 21=HOCl, 22=BrONO2, 23=HOBr,
    //   24=H2CO, 25=CH3O2, 26=CH3O2H, 27=OClO, 28=Cl2O2, 29=BrCl
    pub fn den_get(&self, ib: usize, ispec: usize) -> f64 {
        match ispec {
            0 => self.dno[ib],   1 => self.dno2[ib],  2 => self.dno3[ib],
            3 => self.dn2o5[ib], 4 => self.dhno3[ib], 5 => self.dh[ib],
            6 => self.doh[ib],   7 => self.dho2[ib],  8 => self.dh2o2[ib],
            9 => self.do_[ib],   10 => self.do3[ib],  11 => self.dbro[ib],
            12 => self.dbr[ib],  13 => self.dhbr[ib], 14 => self.dhno2[ib],
            15 => self.dhcl[ib], 16 => self.dcl[ib],  17 => self.dcl2[ib],
            18 => self.dclo[ib], 19 => self.dclno3[ib],20 => self.dhno4[ib],
            21 => self.dhocl[ib],22 => self.dbrno3[ib],23 => self.dhobr[ib],
            24 => self.dh2co[ib],25 => self.droo[ib], 26 => self.drooh[ib],
            27 => self.doclo[ib],28 => self.dcl2o2[ib],29 => self.dbrcl[ib],
            _ => panic!("den_get: species index {} out of range", ispec),
        }
    }

    pub fn den_set(&mut self, ib: usize, ispec: usize, val: f64) {
        match ispec {
            0 => self.dno[ib] = val,   1 => self.dno2[ib] = val,
            2 => self.dno3[ib] = val,  3 => self.dn2o5[ib] = val,
            4 => self.dhno3[ib] = val, 5 => self.dh[ib] = val,
            6 => self.doh[ib] = val,   7 => self.dho2[ib] = val,
            8 => self.dh2o2[ib] = val, 9 => self.do_[ib] = val,
            10 => self.do3[ib] = val,  11 => self.dbro[ib] = val,
            12 => self.dbr[ib] = val,  13 => self.dhbr[ib] = val,
            14 => self.dhno2[ib] = val,15 => self.dhcl[ib] = val,
            16 => self.dcl[ib] = val,  17 => self.dcl2[ib] = val,
            18 => self.dclo[ib] = val, 19 => self.dclno3[ib] = val,
            20 => self.dhno4[ib] = val,21 => self.dhocl[ib] = val,
            22 => self.dbrno3[ib] = val,23 => self.dhobr[ib] = val,
            24 => self.dh2co[ib] = val,25 => self.droo[ib] = val,
            26 => self.drooh[ib] = val,27 => self.doclo[ib] = val,
            28 => self.dcl2o2[ib] = val,29 => self.dbrcl[ib] = val,
            _ => panic!("den_set: species index {} out of range", ispec),
        }
    }

    /// Get FFFFFF(ib, j) — long-lived mixing ratio j (1-based) at box ib (0-based).
    /// FFFFFF(NB,18) is EQUIVALENCED to FO3(1) in Fortran.
    pub fn fff_get(&self, ib: usize, j: usize) -> f64 {
        match j {
            1  => self.fo3[ib],    2  => self.fn2o[ib],  3  => self.fnoy[ib],
            4  => self.fch4[ib],   5  => self.fco[ib],   6  => self.fclx[ib],
            7  => self.fcf2cl[ib], 8  => self.fcfcl3[ib],9  => self.fccl4[ib],
            10 => self.fch3cl[ib], 11 => self.fmecl[ib], 12 => self.fh2[ib],
            13 => self.fh2o[ib],   14 => self.fnh3[ib],  15 => self.fc5h8[ib],
            16 => self.fbrx[ib],   17 => self.fch3br[ib],18 => self.focs[ib],
            _ => 0.0,
        }
    }

    /// Set FFFFFF(ib, j) — long-lived mixing ratio j (1-based) at box ib (0-based).
    pub fn fff_set(&mut self, ib: usize, j: usize, val: f64) {
        match j {
            1  => self.fo3[ib] = val,    2  => self.fn2o[ib] = val,
            3  => self.fnoy[ib] = val,   4  => self.fch4[ib] = val,
            5  => self.fco[ib] = val,    6  => self.fclx[ib] = val,
            7  => self.fcf2cl[ib] = val, 8  => self.fcfcl3[ib] = val,
            9  => self.fccl4[ib] = val,  10 => self.fch3cl[ib] = val,
            11 => self.fmecl[ib] = val,  12 => self.fh2[ib] = val,
            13 => self.fh2o[ib] = val,   14 => self.fnh3[ib] = val,
            15 => self.fc5h8[ib] = val,  16 => self.fbrx[ib] = val,
            17 => self.fch3br[ib] = val, 18 => self.focs[ib] = val,
            _ => {}
        }
    }

    /// Get J-value profile for the k-th J-value (0-based) at altitude level il (0-based).
    /// Equivalent to VVVVVV(IL, K) in Fortran.
    pub fn jval_get(&self, il: usize, k: usize) -> f64 {
        match k {
            0 => self.vno[il],    1 => self.vo2[il],    2 => self.vo3[il],
            3 => self.vo3d[il],   4 => self.vh2coa[il], 5 => self.vh2cob[il],
            6 => self.vh2o2[il],  7 => self.vrooh[il],  8 => self.vno2[il],
            9 => self.vno3x[il],  10 => self.vno3l[il], 11 => self.vn2o5[il],
            12 => self.vhno2[il], 13 => self.vhno3[il], 14 => self.vhno4[il],
            15 => self.vclno3[il],16 => self.vcl2[il],  17 => self.vhocl[il],
            18 => self.voclo[il], 19 => self.vcl2o2[il],20 => self.vclo[il],
            21 => self.vbro[il],  22 => self.vbrno3[il],23 => self.vhobr[il],
            24 => self.vn2o[il],  25 => self.vcfcl3[il],26 => self.vf2cl2[il],
            27 => self.vf113[il], 28 => self.vf114[il], 29 => self.vf115[il],
            30 => self.vccl4[il], 31 => self.vch3cl[il],32 => self.vmecf[il],
            33 => self.vch3br[il],34 => self.vh1211[il],35 => self.vh1301[il],
            36 => self.vh2402[il],37 => self.vh22[il],  38 => self.vh123[il],
            39 => self.vh141b[il],40 => self.vchbr3[il],41 => self.vch3i[il],
            42 => self.vcf3i[il], 43 => self.vocs[il],
            _ => panic!("jval_get: J-value index {} out of range", k),
        }
    }

    pub fn jval_set(&mut self, il: usize, k: usize, val: f64) {
        match k {
            0 => self.vno[il] = val,    1 => self.vo2[il] = val,
            2 => self.vo3[il] = val,    3 => self.vo3d[il] = val,
            4 => self.vh2coa[il] = val, 5 => self.vh2cob[il] = val,
            6 => self.vh2o2[il] = val,  7 => self.vrooh[il] = val,
            8 => self.vno2[il] = val,   9 => self.vno3x[il] = val,
            10 => self.vno3l[il] = val, 11 => self.vn2o5[il] = val,
            12 => self.vhno2[il] = val, 13 => self.vhno3[il] = val,
            14 => self.vhno4[il] = val, 15 => self.vclno3[il] = val,
            16 => self.vcl2[il] = val,  17 => self.vhocl[il] = val,
            18 => self.voclo[il] = val, 19 => self.vcl2o2[il] = val,
            20 => self.vclo[il] = val,  21 => self.vbro[il] = val,
            22 => self.vbrno3[il] = val,23 => self.vhobr[il] = val,
            24 => self.vn2o[il] = val,  25 => self.vcfcl3[il] = val,
            26 => self.vf2cl2[il] = val,27 => self.vf113[il] = val,
            28 => self.vf114[il] = val, 29 => self.vf115[il] = val,
            30 => self.vccl4[il] = val, 31 => self.vch3cl[il] = val,
            32 => self.vmecf[il] = val, 33 => self.vch3br[il] = val,
            34 => self.vh1211[il] = val,35 => self.vh1301[il] = val,
            36 => self.vh2402[il] = val,37 => self.vh22[il] = val,
            38 => self.vh123[il] = val, 39 => self.vh141b[il] = val,
            40 => self.vchbr3[il] = val,41 => self.vch3i[il] = val,
            42 => self.vcf3i[il] = val, 43 => self.vocs[il] = val,
            _ => panic!("jval_set: J-value index {} out of range", k),
        }
    }
}
