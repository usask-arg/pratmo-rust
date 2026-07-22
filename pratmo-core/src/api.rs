/// High-level programmatic API for the PRATMO photochemical box model.
///
/// Typical usage:
/// ```ignore
/// let model = PratmoModel::new("path/to/data");
/// let cfg = DiurnConfig { latitude_deg: 0.0, julian_day: 120, ..Default::default() };
/// let out = model.run_diurn(&cfg)?;
/// let o3_grid = out.species_grid(|s| s.o3); // shape [boxes, timesteps]
/// ```
use std::path::{Path, PathBuf};

use anyhow::{bail, Result};
use ndarray::Array2;

use crate::{
    constants::{LEGACY_NL, NATM, NB, NDEN, NL, NXNOFT},
    ctm::ctmlfq_in_memory,
    diurnal::{diurn, diurn_parallel_boxes},
    path::tpath,
    reader::{setday, setday_cpp, setday_elapsed_hours, FortranReader, ModelReader},
    solver::{fixrat, rplace, splace},
    state::ModelState,
};

// ── Species types ──────────────────────────────────────────────────────────────

/// Number densities (cm⁻³) for the 40 implicit (Newton-Raphson) species.
#[derive(Debug, Clone, Default)]
pub struct ImplicitSpecies {
    pub no: f64,
    pub no2: f64,
    pub no3: f64,
    pub n2o5: f64,
    pub hno3: f64,
    pub h: f64,
    pub oh: f64,
    pub ho2: f64,
    pub h2o2: f64,
    pub o: f64,
    pub o3: f64,
    pub bro: f64,
    pub br: f64,
    pub hbr: f64,
    pub hno2: f64,
    pub hcl: f64,
    pub cl: f64,
    pub cl2: f64,
    pub clo: f64,
    pub clono2: f64,
    pub hno4: f64,
    pub hocl: f64,
    pub brono2: f64,
    pub hobr: f64,
    pub h2co: f64,
    pub ch3o2: f64,
    pub ch3o2h: f64,
    pub oclo: f64,
    pub cl2o2: f64,
    pub brcl: f64,
    // Iodine family
    pub i: f64,
    pub io: f64,
    pub hoi: f64,
    pub iono2: f64,
    pub hi: f64,
    pub oio: f64,
    pub i2: f64,
    pub i2o2: f64,
    pub i2o3: f64,
    pub i2o4: f64,
}

impl ImplicitSpecies {
    fn from_state(s: &ModelState, ib: usize) -> Self {
        Self {
            no: s.dno[ib],
            no2: s.dno2[ib],
            no3: s.dno3[ib],
            n2o5: s.dn2o5[ib],
            hno3: s.dhno3[ib],
            h: s.dh[ib],
            oh: s.doh[ib],
            ho2: s.dho2[ib],
            h2o2: s.dh2o2[ib],
            o: s.do_[ib],
            o3: s.do3[ib],
            bro: s.dbro[ib],
            br: s.dbr[ib],
            hbr: s.dhbr[ib],
            hno2: s.dhno2[ib],
            hcl: s.dhcl[ib],
            cl: s.dcl[ib],
            cl2: s.dcl2[ib],
            clo: s.dclo[ib],
            clono2: s.dclno3[ib],
            hno4: s.dhno4[ib],
            hocl: s.dhocl[ib],
            brono2: s.dbrno3[ib],
            hobr: s.dhobr[ib],
            h2co: s.dh2co[ib],
            ch3o2: s.droo[ib],
            ch3o2h: s.drooh[ib],
            oclo: s.doclo[ib],
            cl2o2: s.dcl2o2[ib],
            brcl: s.dbrcl[ib],
            i: s.di_[ib],
            io: s.dio[ib],
            hoi: s.dhoi[ib],
            iono2: s.diono2[ib],
            hi: s.dhi[ib],
            oio: s.doio[ib],
            i2: s.di2[ib],
            i2o2: s.di2o2[ib],
            i2o3: s.di2o3[ib],
            i2o4: s.di2o4[ib],
        }
    }

    /// Extract from the xxnoft time-series array at time-step kt for box ib.
    /// Species not in the NR active set fall back to the box's final converged value.
    fn from_timeseries(s: &ModelState, ib: usize, kt: usize) -> Self {
        // s.n[k] is the 1-based NR slot for physical species k (0-based).
        // xxnoft[slot-1, kt, ib] holds the value if slot <= ntotx.
        let ntotx = s.ntotx;
        let get = |k: usize| -> f64 {
            let slot = s.n[k]; // 1-based
            if slot > 0 && slot <= ntotx {
                s.xxnoft[[slot - 1, kt, ib]]
            } else {
                s.den_get(ib, k)
            }
        };
        Self {
            no: get(0),
            no2: get(1),
            no3: get(2),
            n2o5: get(3),
            hno3: get(4),
            h: get(5),
            oh: get(6),
            ho2: get(7),
            h2o2: get(8),
            o: get(9),
            o3: get(10),
            bro: get(11),
            br: get(12),
            hbr: get(13),
            hno2: get(14),
            hcl: get(15),
            cl: get(16),
            cl2: get(17),
            clo: get(18),
            clono2: get(19),
            hno4: get(20),
            hocl: get(21),
            brono2: get(22),
            hobr: get(23),
            h2co: get(24),
            ch3o2: get(25),
            ch3o2h: get(26),
            oclo: get(27),
            cl2o2: get(28),
            brcl: get(29),
            i: get(30),
            io: get(31),
            hoi: get(32),
            iono2: get(33),
            hi: get(34),
            oio: get(35),
            i2: get(36),
            i2o2: get(37),
            i2o3: get(38),
            i2o4: get(39),
        }
    }
}

/// Dimensionless mixing ratios for the 19 long-lived species.
/// Used both as initial conditions (input) and daily-mean output.
#[derive(Debug, Clone, Default)]
pub struct LongLivedMixingRatios {
    pub o3: f64,
    pub n2o: f64,
    pub noy: f64,
    pub ch4: f64,
    pub co: f64,
    pub clx: f64,
    pub cf2cl2: f64, // CFC-12
    pub cfcl3: f64,  // CFC-11
    pub ccl4: f64,
    pub ch3cl: f64,
    pub ch3ccl3: f64, // MeCl / CH₃CCl₃
    pub h2: f64,
    pub h2o: f64,
    pub nh3: f64,
    pub c5h8: f64, // isoprene
    pub brx: f64,  // total Bry
    pub ch3br: f64,
    pub ocs: f64,
    pub iodx: f64, // total Iy
}

impl LongLivedMixingRatios {
    fn from_state(s: &ModelState, ib: usize) -> Self {
        Self {
            o3: s.fo3[ib],
            n2o: s.fn2o[ib],
            noy: s.fnoy[ib],
            ch4: s.fch4[ib],
            co: s.fco[ib],
            clx: s.fclx[ib],
            cf2cl2: s.fcf2cl[ib],
            cfcl3: s.fcfcl3[ib],
            ccl4: s.fccl4[ib],
            ch3cl: s.fch3cl[ib],
            ch3ccl3: s.fmecl[ib],
            h2: s.fh2[ib],
            h2o: s.fh2o[ib],
            nh3: s.fnh3[ib],
            c5h8: s.fc5h8[ib],
            brx: s.fbrx[ib],
            ch3br: s.fch3br[ib],
            ocs: s.focs[ib],
            iodx: s.fiodx[ib],
        }
    }

    fn apply_to_state(&self, s: &mut ModelState, ib: usize) {
        s.fo3[ib] = self.o3;
        s.fn2o[ib] = self.n2o;
        s.fnoy[ib] = self.noy;
        s.fch4[ib] = self.ch4;
        s.fco[ib] = self.co;
        s.fclx[ib] = self.clx;
        s.fcf2cl[ib] = self.cf2cl2;
        s.fcfcl3[ib] = self.cfcl3;
        s.fccl4[ib] = self.ccl4;
        s.fch3cl[ib] = self.ch3cl;
        s.fmecl[ib] = self.ch3ccl3;
        s.fh2[ib] = self.h2;
        s.fh2o[ib] = self.h2o;
        s.fnh3[ib] = self.nh3;
        s.fc5h8[ib] = self.c5h8;
        s.fbrx[ib] = self.brx;
        s.fch3br[ib] = self.ch3br;
        s.focs[ib] = self.ocs;
        s.fiodx[ib] = self.iodx;
    }
}

/// Photolysis rates (s⁻¹) for all 52 J-value channels.
#[derive(Debug, Clone, Default)]
pub struct JValues {
    pub no: f64,
    pub o2: f64,
    pub o3: f64,
    pub o3_o1d: f64,
    pub h2co_a: f64,
    pub h2co_b: f64,
    pub h2o2: f64,
    pub rooh: f64,
    pub no2: f64,
    pub no3_x: f64,
    pub no3_l: f64,
    pub n2o5: f64,
    pub hno2: f64,
    pub hno3: f64,
    pub hno4: f64,
    pub clono2: f64,
    pub cl2: f64,
    pub hocl: f64,
    pub oclo: f64,
    pub cl2o2: f64,
    pub clo: f64,
    pub bro: f64,
    pub brono2: f64,
    pub hobr: f64,
    pub n2o: f64,
    pub cfc11: f64, // CFCl₃
    pub cfc12: f64, // CF₂Cl₂
    pub cfc113: f64,
    pub cfc114: f64,
    pub cfc115: f64,
    pub ccl4: f64,
    pub ch3cl: f64,
    pub ch3ccl3: f64, // MeCF
    pub ch3br: f64,
    pub h1211: f64,
    pub h1301: f64,
    pub h2402: f64,
    pub hcfc22: f64,
    pub hcfc123: f64,
    pub hcfc141b: f64,
    pub chbr3: f64,
    pub ch3i: f64,
    pub cf3i: f64,
    pub ocs: f64,
    pub io: f64,
    pub hoi: f64,
    pub iono2: f64,
    pub oio: f64,
    pub i2: f64,
    pub i2o2: f64,
    pub i2o3: f64,
    pub i2o4: f64,
}

impl JValues {
    #[cfg(test)]
    fn from_state_alt(s: &ModelState, il: usize) -> Self {
        Self::from_getter(|index| s.jval_get(il, index))
    }

    fn from_state_box_daily_mean(s: &ModelState, ib: usize) -> Self {
        Self::from_getter(|channel| {
            (1..s.ntimdo)
                .map(|time_index| {
                    let angle_index = (s.jtim[time_index] - 1).max(0) as usize;
                    let weight = (s.dtime[time_index] - s.dtime[time_index - 1]) / s.daysec;
                    weight * s.storjv[[channel, angle_index, ib]]
                })
                .sum()
        })
    }

    fn from_getter(mut get: impl FnMut(usize) -> f64) -> Self {
        Self {
            no: get(0),
            o2: get(1),
            o3: get(2),
            o3_o1d: get(3),
            h2co_a: get(4),
            h2co_b: get(5),
            h2o2: get(6),
            rooh: get(7),
            no2: get(8),
            no3_x: get(9),
            no3_l: get(10),
            n2o5: get(11),
            hno2: get(12),
            hno3: get(13),
            hno4: get(14),
            clono2: get(15),
            cl2: get(16),
            hocl: get(17),
            oclo: get(18),
            cl2o2: get(19),
            clo: get(20),
            bro: get(21),
            brono2: get(22),
            hobr: get(23),
            n2o: get(24),
            cfc11: get(25),
            cfc12: get(26),
            cfc113: get(27),
            cfc114: get(28),
            cfc115: get(29),
            ccl4: get(30),
            ch3cl: get(31),
            ch3ccl3: get(32),
            ch3br: get(33),
            h1211: get(34),
            h1301: get(35),
            h2402: get(36),
            hcfc22: get(37),
            hcfc123: get(38),
            hcfc141b: get(39),
            chbr3: get(40),
            ch3i: get(41),
            cf3i: get(42),
            ocs: get(43),
            io: get(44),
            hoi: get(45),
            iono2: get(46),
            oio: get(47),
            i2: get(48),
            i2o2: get(49),
            i2o3: get(50),
            i2o4: get(51),
        }
    }
}

// ── Config types ───────────────────────────────────────────────────────────────

/// Per-box configuration for a CTM run.
#[derive(Debug, Clone)]
pub struct CtmBoxSpec {
    /// 1-based standard pressure level index (1 = surface, NL = top).
    pub altitude_level: u8,
    /// Generic aerosol surface area density in um2/cm3.
    pub aerosol_surface_area_um2_cm3: f64,
    /// Sea-salt surface area density in um2/cm3 (iodine recycling only).
    pub sea_salt_surface_area_um2_cm3: f64,
    pub temp_offset_k: f64,
}

/// Configuration for a CTM climatological run.
/// This overrides the geographic/temporal fields from fort01.x/boxin_gui.dat.
#[derive(Debug, Clone)]
pub struct CtmConfig {
    pub latitude_deg: f64,
    pub julian_day: u16, // 1..=366
    pub integration_days: u32,
    pub boxes: Vec<CtmBoxSpec>,
    pub bromine: bool,
    pub iodine: bool,
    pub solar_flux_scale: f64,
}

impl Default for CtmConfig {
    fn default() -> Self {
        Self {
            latitude_deg: 60.0,
            julian_day: 75,
            integration_days: 40,
            boxes: Vec::new(),
            bromine: false,
            iodine: true,
            solar_flux_scale: 1.0,
        }
    }
}

/// Per-box configuration for a DIURN run.
#[derive(Debug, Clone)]
pub struct DiurnBoxSpec {
    /// 1-based standard pressure level index (1 = surface, NL = top).
    pub altitude_level: u8,
    /// Optional exact chemistry-box altitude in km. When supplied with a
    /// custom atmosphere, chemistry is evaluated at this altitude while the
    /// actinic flux is linearly interpolated between the surrounding
    /// radiative levels.
    pub altitude_km: Option<f64>,
    /// Generic aerosol surface area density in um2/cm3.
    pub aerosol_surface_area_um2_cm3: f64,
    /// Sea-salt surface area density in um2/cm3 (iodine recycling only).
    pub sea_salt_surface_area_um2_cm3: f64,
    pub temp_offset_k: f64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum O3InputKind {
    MixingRatio,
    NumberDensityCm3,
}

/// Custom vertical atmosphere for DIURN runs.
///
/// Pressures are in mb, temperatures in K. O3 values are either dimensionless
/// mixing ratios or number densities depending on `o3_kind`. If `altitude_km`
/// is omitted, altitudes are estimated hydrostatically from pressure/temperature.
#[derive(Debug, Clone)]
pub struct CustomAtmosphereProfile {
    pub pressure_mb: Vec<f64>,
    pub temperature_k: Vec<f64>,
    pub altitude_km: Option<Vec<f64>>,
    pub o3: Vec<f64>,
    pub o3_kind: O3InputKind,
    /// Optional sulfate aerosol surface-area profile (um2/cm3). This is used
    /// both by heterogeneous chemistry and, when enabled, radiative transfer.
    pub aerosol_surface_area_um2_cm3: Option<Vec<f64>>,
}

/// Configuration for a diurnal cycle run.
///
/// `initial_mixing_ratios` overrides the long-lived species mixing ratios per box.
/// If `None`, the values from fort02.x (if present) or fort01.x defaults are used.
#[derive(Debug, Clone)]
pub struct DiurnConfig {
    pub latitude_deg: f64,
    pub julian_day: u16,
    pub integration_days: u32,
    pub boxes: Vec<DiurnBoxSpec>,
    pub bromine: bool,
    pub iodine: bool,
    pub parallel_boxes: bool,
    /// Match the later C++ box model's fixed-mu grid and daily endpoint
    /// convergence algorithm.
    pub cpp_compatibility: bool,
    /// Optional explicit integration coordinate in elapsed hours after local
    /// noon. When present, this overrides the legacy or C++ generated grid.
    pub elapsed_time_hours: Option<Vec<f64>>,
    pub solar_flux_scale: f64,
    /// Lambertian lower-boundary reflectivity used by the photolysis solver.
    pub surface_albedo: f64,
    /// Enable sulfate-aerosol and sea-salt heterogeneous chemistry.
    pub heterogeneous_chemistry: bool,
    /// Include sulfate aerosol extinction and scattering in photolysis.
    pub radiative_aerosol: bool,
    /// Optional custom pressure/temperature/O3 grid. If provided and `boxes` is
    /// empty, one DIURN box is created for each profile level.
    pub atmosphere: Option<CustomAtmosphereProfile>,
    /// One `LongLivedMixingRatios` per box; must match `boxes.len()` if provided.
    pub initial_mixing_ratios: Option<Vec<LongLivedMixingRatios>>,
}

impl Default for DiurnConfig {
    fn default() -> Self {
        Self {
            latitude_deg: 0.0,
            julian_day: 120,
            integration_days: 20,
            boxes: Vec::new(),
            bromine: false,
            iodine: true,
            parallel_boxes: false,
            cpp_compatibility: false,
            elapsed_time_hours: None,
            solar_flux_scale: 1.0,
            surface_albedo: 0.20,
            heterogeneous_chemistry: true,
            radiative_aerosol: false,
            atmosphere: None,
            initial_mixing_ratios: None,
        }
    }
}

#[derive(Debug, Clone)]
pub struct No2ConstrainedDiurnConfig {
    pub diurn: DiurnConfig,
    pub observed_no2_cm3: Vec<f64>,
    pub target_hhmm: i32,
    pub iterations: usize,
}

#[derive(Debug, Clone)]
pub struct No2ConstrainedDiurnOutput {
    pub output: DiurnOutput,
    pub noy_scale: Vec<f64>,
    pub modeled_no2_cm3: Vec<f64>,
}

// ── Output types ───────────────────────────────────────────────────────────────

/// Run diagnostics available from both CTM and DIURN modes.
#[derive(Debug, Clone, Default)]
pub struct Diagnostics {
    pub raxloop: f64,
    pub radcount: f64,
    pub newraf_nonconvergence_count: usize,
    pub rafday_nonconvergence_count: usize,
    /// Maximum final RAFDAY relative-correction metric across calls and boxes.
    pub rafday_max_final_relative_correction: f64,
    /// Maximum number of RAFDAY Newton corrections across calls and boxes.
    pub rafday_max_correction_iterations: usize,
}

/// Snapshot of a single box's state (daily mean or final converged value).
#[derive(Debug, Clone)]
pub struct BoxSnapshot {
    /// 0-based box index.
    pub box_index: usize,
    pub altitude_km: f64,
    pub pressure_mb: f64,
    pub temperature_k: f64,
    /// Air number density (cm⁻³).
    pub air_density_cm3: f64,
    /// Implicit species number densities (cm⁻³).
    pub implicit: ImplicitSpecies,
    /// Long-lived species mixing ratios (dimensionless).
    pub long_lived: LongLivedMixingRatios,
    /// Photolysis rates (s⁻¹).
    pub jvalues: JValues,
}

/// One time-step in a diurnal time series.
#[derive(Debug, Clone)]
pub struct DiurnTimeStep {
    /// Monotonic seconds since the noon start of the integrated 24-hour orbit.
    pub elapsed_seconds: f64,
    /// Local clock label. Both orbit endpoints are noon (`1200`).
    pub time_hhmm: i32,
    pub implicit: ImplicitSpecies,
}

/// Full diurnal time series for one box.
#[derive(Debug, Clone)]
pub struct DiurnBoxTimeSeries {
    pub box_index: usize,
    pub altitude_km: f64,
    pub pressure_mb: f64,
    pub steps: Vec<DiurnTimeStep>,
}

/// Output from a CTM climatological run.
#[derive(Debug, Clone)]
pub struct CtmOutput {
    /// Final converged snapshot for each box.
    pub boxes: Vec<BoxSnapshot>,
    pub diagnostics: Diagnostics,
}

/// Output from a diurnal cycle run.
#[derive(Debug, Clone)]
pub struct DiurnOutput {
    /// Daily-mean snapshot for each box (equivalent to fort08.x content).
    pub boxes: Vec<BoxSnapshot>,
    /// Full diurnal time series for each box (equivalent to fort07.x content).
    pub time_series: Vec<DiurnBoxTimeSeries>,
    pub diagnostics: Diagnostics,
}

impl DiurnOutput {
    /// Extract one implicit species as an `(n_boxes × n_timesteps)` array.
    ///
    /// Example: `output.species_grid(|s| s.o3)` → O₃ as `Array2<f64>`.
    pub fn species_grid(&self, f: impl Fn(&ImplicitSpecies) -> f64) -> Array2<f64> {
        let nboxes = self.time_series.len();
        if nboxes == 0 {
            return Array2::zeros((0, 0));
        }
        let ntimes = self
            .time_series
            .iter()
            .map(|ts| ts.steps.len())
            .max()
            .unwrap_or(0);
        let mut arr = Array2::zeros((nboxes, ntimes));
        for (ib, ts) in self.time_series.iter().enumerate() {
            for (it, step) in ts.steps.iter().enumerate() {
                arr[[ib, it]] = f(&step.implicit);
            }
        }
        arr
    }

    /// Extract one J-value as an `(n_boxes × n_timesteps)` array.
    ///
    /// J-values are constant within the day for each box (daily-mean level),
    /// so each row is filled with the box's daily-mean value.
    ///
    /// Example: `output.jvalue_grid(|j| j.no2)` → J(NO₂) as `Array2<f64>`.
    pub fn jvalue_grid(&self, f: impl Fn(&JValues) -> f64) -> Array2<f64> {
        let nboxes = self.boxes.len();
        let ntimes = self
            .time_series
            .iter()
            .map(|ts| ts.steps.len())
            .max()
            .unwrap_or(0);
        let mut arr = Array2::zeros((nboxes, ntimes));
        for (ib, snap) in self.boxes.iter().enumerate() {
            let val = f(&snap.jvalues);
            for it in 0..ntimes {
                arr[[ib, it]] = val;
            }
        }
        arr
    }
}

fn hhmm_to_minutes(hhmm: i32) -> Option<i32> {
    let hours = hhmm.div_euclid(100);
    let minutes = hhmm.rem_euclid(100);
    ((0..24).contains(&hours) && minutes < 60).then_some(hours * 60 + minutes)
}

fn time_distance_hhmm(a: i32, b: i32) -> i32 {
    // Parity mode preserves a legacy rounding quirk that can encode 02:00 as
    // 0160. Carry such model minutes rather than panicking in the structured
    // observation-matching API.
    let a_minutes = hhmm_to_minutes(a)
        .unwrap_or_else(|| (a.div_euclid(100) * 60 + a.rem_euclid(100)).rem_euclid(24 * 60));
    let b_minutes = hhmm_to_minutes(b).expect("validated target must be a valid HHMM value");
    let direct = (a_minutes - b_minutes).abs();
    direct.min(24 * 60 - direct)
}

fn no2_at_hhmm(out: &DiurnOutput, target_hhmm: i32) -> Result<Vec<f64>> {
    out.time_series
        .iter()
        .map(|ts| {
            ts.steps
                .iter()
                .min_by_key(|step| time_distance_hhmm(step.time_hhmm, target_hhmm))
                .map(|step| step.implicit.no2)
                .ok_or_else(|| anyhow::anyhow!("DIURN output has an empty time series"))
        })
        .collect()
}

impl CtmOutput {
    /// Extract one implicit species as a `(n_boxes,)` altitude profile.
    ///
    /// Example: `output.species_profile(|s| s.o3)`.
    pub fn species_profile(&self, f: impl Fn(&ImplicitSpecies) -> f64) -> Vec<f64> {
        self.boxes.iter().map(|b| f(&b.implicit)).collect()
    }

    /// Extract one J-value as a `(n_boxes,)` altitude profile.
    pub fn jvalue_profile(&self, f: impl Fn(&JValues) -> f64) -> Vec<f64> {
        self.boxes.iter().map(|b| f(&b.jvalues)).collect()
    }

    /// Extract one long-lived mixing ratio as a `(n_boxes,)` altitude profile.
    pub fn long_lived_profile(&self, f: impl Fn(&LongLivedMixingRatios) -> f64) -> Vec<f64> {
        self.boxes.iter().map(|b| f(&b.long_lived)).collect()
    }
}

// ── PratmoModel ────────────────────────────────────────────────────────────────

/// Entry point for the PRATMO photochemical box model.
///
/// Construct with [`PratmoModel::with_defaults()`] for a self-contained run using
/// compiled-in science data, or [`PratmoModel::new(data_dir)`] to load from the
/// original Fortran-format files (useful for backwards-compatibility testing or
/// overriding individual files).
///
/// For file-based runs the directory must contain:
/// `fort01.x`, `fort10_cam06.x`, `fort11_jpl09.x`, `fort13.x`, `fort14.x`, `J_H2O_SZA0.dat`.
/// CTM climatology runs additionally need `fort03_LLM.x`, `fort05.x`, `fort51.x`.
pub struct PratmoModel {
    /// `None` → use compiled-in embedded data; `Some(path)` → read from files.
    data_dir: Option<PathBuf>,
}

impl PratmoModel {
    /// Create a model using only compiled-in default science data.
    /// No data files on disk are required.
    pub fn with_defaults() -> Self {
        Self { data_dir: None }
    }

    /// Create a model that reads science data from `data_dir` at runtime.
    /// Keeps the original Fortran-file workflow intact for backwards-compatibility
    /// testing and for overriding individual data files.
    pub fn new(data_dir: impl AsRef<Path>) -> Self {
        Self {
            data_dir: Some(data_dir.as_ref().to_owned()),
        }
    }

    /// Run the diurnal cycle (DIURN + TPATH) mode.
    ///
    /// Loads base chemistry setup, then overrides geographic/temporal/box
    /// parameters from `cfg`. Returns structured output without writing any files.
    pub fn run_diurn(&self, cfg: &DiurnConfig) -> Result<DiurnOutput> {
        validate_diurn_config(cfg)?;
        let mut s = self.load_base_state()?;
        apply_diurn_config(&mut s, cfg)?;
        validate_supported_mode(&s)?;

        s.out_unit7 = None;
        s.out_unit8 = None;
        s.out_unit9 = None;

        if cfg.parallel_boxes {
            diurn_parallel_boxes(&mut s)?;
        } else {
            diurn(&mut s)?;
        }
        tpath(&mut s)?;

        Ok(extract_diurn_output(&s))
    }

    pub fn run_diurn_no2_constrained(
        &self,
        cfg: &No2ConstrainedDiurnConfig,
    ) -> Result<No2ConstrainedDiurnOutput> {
        validate_diurn_config(&cfg.diurn)?;
        if hhmm_to_minutes(cfg.target_hhmm).is_none() {
            bail!(
                "target_hhmm must be a valid local-time HHMM value from 0000 to 2359, got {}",
                cfg.target_hhmm
            );
        }
        let nbox = diurn_config_nbox(&cfg.diurn);
        if cfg.observed_no2_cm3.len() != nbox {
            bail!(
                "observed_no2_cm3 length ({}) must match DIURN box count ({})",
                cfg.observed_no2_cm3.len(),
                nbox
            );
        }
        for (ib, &obs) in cfg.observed_no2_cm3.iter().enumerate() {
            if !(obs.is_finite() && obs >= 0.0) {
                bail!("observed_no2_cm3[{ib}] must be finite and non-negative");
            }
        }

        let mut base = self.load_base_state()?;
        apply_diurn_config(&mut base, &cfg.diurn)?;
        let mut base_mr: Vec<LongLivedMixingRatios> = (0..nbox)
            .map(|ib| LongLivedMixingRatios::from_state(&base, ib))
            .collect();
        if let Some(ref init) = cfg.diurn.initial_mixing_ratios {
            for (dst, src) in base_mr.iter_mut().zip(init.iter()) {
                *dst = src.clone();
            }
        }

        let mut noy_scale = vec![1.0_f64; nbox];
        for _ in 0..cfg.iterations {
            let mut run_cfg = cfg.diurn.clone();
            let mut init = base_mr.clone();
            for (ib, mr) in init.iter_mut().enumerate() {
                mr.noy = base_mr[ib].noy * noy_scale[ib];
            }
            run_cfg.initial_mixing_ratios = Some(init);
            let out = self.run_diurn(&run_cfg)?;
            let modeled_no2 = no2_at_hhmm(&out, cfg.target_hhmm)?;
            for ib in 0..nbox {
                if modeled_no2[ib] > 0.0 && cfg.observed_no2_cm3[ib].is_finite() {
                    noy_scale[ib] *= cfg.observed_no2_cm3[ib] / modeled_no2[ib];
                }
            }
        }

        let mut final_cfg = cfg.diurn.clone();
        let mut init = base_mr.clone();
        for (ib, mr) in init.iter_mut().enumerate() {
            mr.noy = base_mr[ib].noy * noy_scale[ib];
        }
        final_cfg.initial_mixing_ratios = Some(init);
        let final_out = self.run_diurn(&final_cfg)?;
        let modeled_no2 = no2_at_hhmm(&final_out, cfg.target_hhmm)?;

        Ok(No2ConstrainedDiurnOutput {
            output: final_out,
            noy_scale,
            modeled_no2_cm3: modeled_no2,
        })
    }

    /// Run the CTM climatological mode.
    ///
    /// Loads base chemistry setup, then overrides geographic/temporal/box
    /// parameters from `cfg`. Returns structured output without writing any files.
    /// The standard climatology is compiled in for [`PratmoModel::with_defaults`].
    /// A file-based model reads `fort03_LLM.x`, `fort05.x`, and `fort51.x` from
    /// its configured directory instead.
    pub fn run_ctm(&self, cfg: &CtmConfig) -> Result<CtmOutput> {
        validate_ctm_config(cfg)?;
        let mut s = self.load_base_state()?;
        apply_ctm_config(&mut s, cfg);
        validate_supported_mode(&s)?;

        s.out_unit7 = None;
        s.out_unit8 = None;
        s.out_unit9 = None;

        ctmlfq_in_memory(
            &mut s,
            cfg.latitude_deg,
            cfg.julian_day,
            self.data_dir.is_none(),
        )?;

        Ok(extract_ctm_output(&s))
    }

    fn load_base_state(&self) -> Result<Box<ModelState>> {
        let mut s = ModelState::new();
        let mut reader = match &self.data_dir {
            Some(dir) => {
                s.cinpdir = dir.to_string_lossy().into_owned();
                FortranReader::new(dir)
            }
            None => FortranReader::embedded(),
        };
        reader.read_all(&mut s)?;
        Ok(s)
    }
}

fn validate_common_config(
    mode: &str,
    latitude_deg: f64,
    julian_day: u16,
    integration_days: u32,
    solar_flux_scale: f64,
) -> Result<()> {
    if !(latitude_deg.is_finite() && (-90.0..=90.0).contains(&latitude_deg)) {
        bail!("{mode} latitude_deg must be finite and between -90 and 90");
    }
    if !(1..=366).contains(&julian_day) {
        bail!("{mode} julian_day must be between 1 and 366, got {julian_day}");
    }
    if integration_days == 0 || integration_days > i32::MAX as u32 {
        bail!("{mode} integration_days must be between 1 and {}", i32::MAX);
    }
    if !(solar_flux_scale.is_finite() && solar_flux_scale > 0.0) {
        bail!("{mode} solar_flux_scale must be finite and positive");
    }
    Ok(())
}

fn validate_box_spec(
    mode: &str,
    index: usize,
    altitude_level: u8,
    altitude_km: Option<f64>,
    max_altitude_level: usize,
    aerosol_surface_area_um2_cm3: f64,
    sea_salt_surface_area_um2_cm3: f64,
    temp_offset_k: f64,
) -> Result<()> {
    if altitude_level == 0 || altitude_level as usize > max_altitude_level {
        bail!(
            "{mode} boxes[{index}].altitude_level must be between 1 and {max_altitude_level}, got {altitude_level}"
        );
    }
    if let Some(altitude_km) = altitude_km {
        if !(altitude_km.is_finite() && altitude_km >= 0.0) {
            bail!("{mode} boxes[{index}].altitude_km must be finite and non-negative");
        }
    }
    if !(aerosol_surface_area_um2_cm3.is_finite() && aerosol_surface_area_um2_cm3 >= 0.0) {
        bail!("{mode} boxes[{index}].aerosol_surface_area_um2_cm3 must be finite and non-negative");
    }
    if !(sea_salt_surface_area_um2_cm3.is_finite() && sea_salt_surface_area_um2_cm3 >= 0.0) {
        bail!(
            "{mode} boxes[{index}].sea_salt_surface_area_um2_cm3 must be finite and non-negative"
        );
    }
    if !temp_offset_k.is_finite() {
        bail!("{mode} boxes[{index}].temp_offset_k must be finite");
    }
    Ok(())
}

fn validate_initial_mixing_ratios(ratios: &[LongLivedMixingRatios]) -> Result<()> {
    for (index, ratio) in ratios.iter().enumerate() {
        let values = [
            ("o3", ratio.o3),
            ("n2o", ratio.n2o),
            ("noy", ratio.noy),
            ("ch4", ratio.ch4),
            ("co", ratio.co),
            ("clx", ratio.clx),
            ("cf2cl2", ratio.cf2cl2),
            ("cfcl3", ratio.cfcl3),
            ("ccl4", ratio.ccl4),
            ("ch3cl", ratio.ch3cl),
            ("ch3ccl3", ratio.ch3ccl3),
            ("h2", ratio.h2),
            ("h2o", ratio.h2o),
            ("nh3", ratio.nh3),
            ("c5h8", ratio.c5h8),
            ("brx", ratio.brx),
            ("ch3br", ratio.ch3br),
            ("ocs", ratio.ocs),
            ("iodx", ratio.iodx),
        ];
        for (name, value) in values {
            if !(value.is_finite() && value >= 0.0) {
                bail!(
                    "DIURN initial_mixing_ratios[{index}].{name} must be finite and non-negative"
                );
            }
        }
    }
    Ok(())
}

fn validate_diurn_config(cfg: &DiurnConfig) -> Result<()> {
    validate_common_config(
        "DIURN",
        cfg.latitude_deg,
        cfg.julian_day,
        cfg.integration_days,
        cfg.solar_flux_scale,
    )?;
    if !(cfg.surface_albedo.is_finite() && (0.0..=1.0).contains(&cfg.surface_albedo)) {
        bail!("DIURN surface_albedo must be finite and between 0 and 1");
    }
    if let Some(hours) = &cfg.elapsed_time_hours {
        if hours.len() < 2 || hours.len() > NXNOFT {
            bail!("DIURN elapsed_time_hours must contain between 2 and {NXNOFT} points");
        }
        if (hours[0] - 0.0).abs() > 1.0e-10 || (hours[hours.len() - 1] - 24.0).abs() > 1.0e-10 {
            bail!("DIURN elapsed_time_hours must start at 0 and end at 24");
        }
        for (index, &value) in hours.iter().enumerate() {
            if !value.is_finite() || !(0.0..=24.0).contains(&value) {
                bail!("DIURN elapsed_time_hours[{index}] must be finite and between 0 and 24");
            }
            if index > 0 && value <= hours[index - 1] {
                bail!("DIURN elapsed_time_hours must be strictly increasing");
            }
        }
    }
    if cfg.boxes.len() > NB {
        bail!("DIURN supports at most {NB} boxes, got {}", cfg.boxes.len());
    }
    let max_altitude_level = cfg
        .atmosphere
        .as_ref()
        .map_or(LEGACY_NL, |profile| profile.pressure_mb.len());
    if cfg.boxes.is_empty() {
        let profile_levels = cfg
            .atmosphere
            .as_ref()
            .map(|profile| profile.pressure_mb.len())
            .unwrap_or(0);
        if profile_levels == 0 {
            bail!("DIURN requires at least one box or a non-empty custom atmosphere");
        }
        if profile_levels > NB {
            bail!(
                "a custom atmosphere used without explicit boxes may contain at most {NB} levels, got {profile_levels}"
            );
        }
    }
    for (index, box_) in cfg.boxes.iter().enumerate() {
        validate_box_spec(
            "DIURN",
            index,
            box_.altitude_level,
            box_.altitude_km,
            max_altitude_level,
            box_.aerosol_surface_area_um2_cm3,
            box_.sea_salt_surface_area_um2_cm3,
            box_.temp_offset_k,
        )?;
        if let Some(altitude_km) = box_.altitude_km {
            let profile = cfg.atmosphere.as_ref().ok_or_else(|| {
                anyhow::anyhow!("DIURN boxes[{index}].altitude_km requires a custom atmosphere")
            })?;
            let profile_altitude = profile.altitude_km.as_ref().ok_or_else(|| {
                anyhow::anyhow!("DIURN boxes[{index}].altitude_km requires atmosphere.altitude_km")
            })?;
            if let (Some(&bottom), Some(&top)) = (profile_altitude.first(), profile_altitude.last())
            {
                if altitude_km < bottom || altitude_km > top {
                    bail!(
                        "DIURN boxes[{index}].altitude_km must lie within the custom atmosphere ({bottom}..={top} km)"
                    );
                }
            }
        }
    }
    if let Some(ratios) = &cfg.initial_mixing_ratios {
        let nbox = diurn_config_nbox(cfg);
        if ratios.len() != nbox {
            bail!(
                "DIURN initial_mixing_ratios length ({}) must match box count ({nbox})",
                ratios.len()
            );
        }
        validate_initial_mixing_ratios(ratios)?;
    }
    Ok(())
}

fn validate_ctm_config(cfg: &CtmConfig) -> Result<()> {
    validate_common_config(
        "CTM",
        cfg.latitude_deg,
        cfg.julian_day,
        cfg.integration_days,
        cfg.solar_flux_scale,
    )?;
    if cfg.boxes.is_empty() {
        bail!("CTM requires at least one box");
    }
    if cfg.boxes.len() > NB {
        bail!("CTM supports at most {NB} boxes, got {}", cfg.boxes.len());
    }
    for (index, box_) in cfg.boxes.iter().enumerate() {
        validate_box_spec(
            "CTM",
            index,
            box_.altitude_level,
            None,
            LEGACY_NL,
            box_.aerosol_surface_area_um2_cm3,
            box_.sea_salt_surface_area_um2_cm3,
            box_.temp_offset_k,
        )?;
    }
    Ok(())
}

/// Reject Fortran modes that are not silently approximated by the public Rust
/// API.  A warning or no-op would make a parity run look successful while
/// producing a different scientific calculation.
fn validate_supported_mode(s: &ModelState) -> Result<()> {
    if s.nd216 < 0 {
        bail!(
            "DERIVS mode (nd216={}) is not implemented; sensitivity Jacobians are unavailable",
            s.nd216
        );
    }
    if s.npstd > 0 {
        bail!(
            "PZSTD mode (NPSTD={}) is not implemented; pressure-to-standard-Z conversion is unavailable",
            s.npstd
        );
    }
    Ok(())
}

// ── Config application ─────────────────────────────────────────────────────────

fn diurn_config_nbox(cfg: &DiurnConfig) -> usize {
    if let Some(profile) = &cfg.atmosphere {
        if cfg.boxes.is_empty() {
            profile.pressure_mb.len().min(NB)
        } else {
            cfg.boxes.len().min(NB)
        }
    } else {
        cfg.boxes.len().min(NB)
    }
}

fn apply_diurn_config(s: &mut ModelState, cfg: &DiurnConfig) -> Result<()> {
    if let Some(profile) = &cfg.atmosphere {
        apply_custom_atmosphere(s, profile)?;
    }

    let nbox = diurn_config_nbox(cfg);
    s.nbox = nbox;
    s.nd216 = 0;
    s.nd216s = 0;

    // Geographic parameters
    s.xlatd = cfg.latitude_deg;
    s.xlat = cfg.latitude_deg.to_radians();

    // Solar declination from Julian day (OSIRIS formula, same as ctmlfq)
    let pi = std::f64::consts::PI;
    let xjd = 2.0 * pi * cfg.julian_day as f64 / 365.0;
    let decang = 6.918e-3 - 0.399912 * xjd.cos() + 0.070257 * xjd.sin()
        - 6.758e-3 * (2.0 * xjd).cos()
        + 9.07e-4 * (2.0 * xjd).sin()
        - 2.697e-3 * (3.0 * xjd).cos()
        + 1.480e-3 * (3.0 * xjd).sin();
    s.xdecd = decang * 57.29578;
    s.xdec = decang;

    // Earth–Sun distance (monthly table, same as ctmlfq)
    const EDIST: [f64; 12] = [
        0.9837, 0.9875, 0.9945, 1.0032, 1.0109, 1.0158, 1.0165, 1.0128, 1.0057, 0.9970, 0.9892,
        0.9842,
    ];
    let mon = ((cfg.julian_day as i32 - 1) / 30).min(11) as usize;
    s.flscal = cfg.solar_flux_scale / (EDIST[mon] * EDIST[mon]);
    s.clouds = cfg.surface_albedo;
    s.heterogeneous_chemistry = cfg.heterogeneous_chemistry;
    s.cpp_compatibility = cfg.cpp_compatibility;
    s.radiative_aerosol = cfg.radiative_aerosol;

    // Recompute diurnal time grid with updated lat/dec
    if let Some(hours) = &cfg.elapsed_time_hours {
        setday_elapsed_hours(s, hours);
    } else if cfg.cpp_compatibility {
        setday_cpp(s, 48);
    } else {
        setday(s);
    }

    s.nday = 1; // full 24-hour integration
    s.ndaysd = cfg.integration_days as i32;
    s.lbrom = cfg.bromine;

    // Box configuration.
    // fort02.x initial densities were rescaled by dm[ialt_old] during read_all.
    // Rescale to dm[ialt_new] before changing nboxdo so diurn() starts from
    // physically consistent densities at the user's requested altitude.
    let ndval = s.ndval as usize;
    let custom_reference = if cfg.atmosphere.is_some() {
        let ialt_ref = (s.nboxdo[0].unsigned_abs() as usize)
            .saturating_sub(1)
            .min(NL - 1);
        let mut densities = [0.0_f64; NDEN];
        for (id, value) in densities.iter_mut().take(ndval).enumerate() {
            *value = s.den_get(0, id);
        }
        Some((ialt_ref, densities))
    } else {
        None
    };
    for ib in 0..nbox {
        let spec = cfg.boxes.get(ib);
        let ialt_old = (s.nboxdo[ib].unsigned_abs() as usize)
            .saturating_sub(1)
            .min(NL - 1);
        let (ialt_new, flux_lower, flux_upper, flux_upper_weight) = if cfg.atmosphere.is_some() {
            if let Some(altitude_km) = spec.and_then(|b| b.altitude_km) {
                append_interpolated_chemistry_level(s, altitude_km)?
            } else {
                let level = spec.map(|b| b.altitude_level).unwrap_or((ib + 1) as u8) as usize;
                if level == 0 {
                    bail!("DIURN box altitude_level must be 1-based, got 0 for box {ib}");
                }
                let ialt = level - 1;
                if let Some(profile) = &cfg.atmosphere {
                    if ialt >= profile.pressure_mb.len() {
                        bail!(
                            "DIURN box altitude_level {} exceeds custom atmosphere level count {}",
                            level,
                            profile.pressure_mb.len()
                        );
                    }
                }
                let ialt = ialt.min(NL - 1);
                (ialt, ialt, ialt, 0.0)
            }
        } else {
            let ialt = (spec.map(|b| b.altitude_level).unwrap_or((ib + 1) as u8) as usize)
                .saturating_sub(1)
                .min(NL - 1);
            (ialt, ialt, ialt, 0.0)
        };
        if let Some((ialt_ref, densities)) = custom_reference {
            if s.dm[ialt_ref] > 0.0 {
                let scale = s.dm[ialt_new] / s.dm[ialt_ref];
                for (id, value) in densities.iter().take(ndval).enumerate() {
                    s.den_set(ib, id, value * scale);
                }
            }
        } else if ialt_old != ialt_new && s.dm[ialt_old] > 0.0 {
            let scale = s.dm[ialt_new] / s.dm[ialt_old];
            for id in 0..ndval {
                let v = s.den_get(ib, id);
                s.den_set(ib, id, v * scale);
            }
        }
        s.nboxdo[ib] = (ialt_new + 1) as i32;
        s.box_flux_lower[ib] = flux_lower;
        s.box_flux_upper[ib] = flux_upper;
        s.box_flux_upper_weight[ib] = flux_upper_weight;
        s.boxaa[ib] = spec.map(|b| b.aerosol_surface_area_um2_cm3).unwrap_or(0.25);
        s.boxss[ib] = spec.map(|b| b.sea_salt_surface_area_um2_cm3).unwrap_or(0.0);
        s.boxtt[ib] = spec.map(|b| b.temp_offset_k).unwrap_or(0.0);
        s.nboxmx[ib] = cfg.integration_days as i32;
        s.nboxwt[ib] = 1;
        s.nboxpr[ib] = 0;
        s.nboxct[ib] = 0;
        // BOXRN is the rainout switch/rate multiplier, not a box identifier.
        // The later C++ box model leaves H2O2RAINOUT at zero by default, so a
        // structured API run must not synthesize rainout from the box index.
        // Doing so preferentially removed the NOy reservoirs and made an
        // otherwise identical atmosphere depend on its position in `boxes`.
        s.boxrn[ib] = 0.0;
        if cfg.atmosphere.is_some() && s.dm[ialt_new] > 0.0 {
            s.do3[ib] = s.do3ref[ialt_new];
            s.fo3[ib] = s.do3[ib] / s.dm[ialt_new];
        }
    }
    // Zero out any leftover boxes from fort01
    for ib in nbox..NB {
        s.nboxdo[ib] = 0;
        s.nboxwt[ib] = 0;
        s.box_flux_lower[ib] = 0;
        s.box_flux_upper[ib] = 0;
        s.box_flux_upper_weight[ib] = 0.0;
    }

    // Initial mixing ratios override
    if let Some(ref init_mr) = cfg.initial_mixing_ratios {
        for (ib, mr) in init_mr.iter().take(nbox).enumerate() {
            mr.apply_to_state(s, ib);
            let ialt = (s.nboxdo[ib].unsigned_abs() as usize)
                .saturating_sub(1)
                .min(NATM - 1);
            s.do3[ib] = s.fo3[ib] * s.dm[ialt];
        }
    }

    if !cfg.iodine {
        disable_iodine(s, cfg.atmosphere.is_none());
    }

    if cfg.atmosphere.is_some() || cfg.initial_mixing_ratios.is_some() {
        for ib in 0..nbox {
            reconcile_custom_box_implicit_species(s, ib);
        }
    }

    Ok(())
}

fn apply_custom_atmosphere(s: &mut ModelState, profile: &CustomAtmosphereProfile) -> Result<()> {
    let n = profile.pressure_mb.len();
    if n == 0 || n > NL {
        bail!("custom atmosphere must contain 1..={NL} levels, got {n}");
    }
    if profile.temperature_k.len() != n || profile.o3.len() != n {
        bail!("custom atmosphere pressure, temperature, and O3 arrays must have equal length");
    }
    if let Some(altitude_km) = &profile.altitude_km {
        if altitude_km.len() != n {
            bail!("custom atmosphere altitude_km length must match pressure length");
        }
    }
    if let Some(aerosol) = &profile.aerosol_surface_area_um2_cm3 {
        if aerosol.len() != n {
            bail!("custom atmosphere aerosol surface-area length must match pressure length");
        }
    }

    let cboltz = 1.38e-19_f64;
    let po2 = s.po2;
    s.nc = n;
    s.nlev = n;

    for i in 0..n {
        let p = profile.pressure_mb[i];
        let t = profile.temperature_k[i];
        let o3 = profile.o3[i];
        if !(p.is_finite() && p > 0.0) {
            bail!("custom atmosphere pressure at index {i} must be positive");
        }
        if !(t.is_finite() && t > 0.0) {
            bail!("custom atmosphere temperature at index {i} must be positive");
        }
        if !(o3.is_finite() && o3 >= 0.0) {
            bail!("custom atmosphere O3 at index {i} must be non-negative");
        }
        if i > 0 && p >= profile.pressure_mb[i - 1] {
            bail!("custom atmosphere pressure must decrease with altitude");
        }

        s.pstd[i] = p;
        s.t[i] = t;
        s.dm[i] = p / (cboltz * t);
        s.theta[i] = t * (1000.0 / p).powf(0.2857);
        s.do3ref[i] = match profile.o3_kind {
            O3InputKind::MixingRatio => o3 * s.dm[i],
            O3InputKind::NumberDensityCm3 => o3,
        };
        s.refo3[i] = s.do3ref[i];
        let aerosol_area = match &profile.aerosol_surface_area_um2_cm3 {
            Some(profile) => {
                let value = profile[i];
                if !(value.is_finite() && value >= 0.0) {
                    bail!("custom atmosphere aerosol surface area at index {i} must be finite and non-negative");
                }
                value
            }
            None => -1.0,
        };
        s.asa[i] = aerosol_area;
        // Original PRATMO conversion: ASA (um2/cm3) × 10^-8 converts
        // surface area to cm2/cm3, then Qext/4 with Qext=2.7 gives the
        // 300-nm extinction coefficient.
        s.aer[i] = aerosol_area.max(0.0) * 0.25 * 2.7e-8;
    }

    if let Some(altitude_km) = &profile.altitude_km {
        for i in 0..n {
            let zkm = altitude_km[i];
            if !(zkm.is_finite() && zkm >= 0.0) {
                bail!("custom atmosphere altitude_km at index {i} must be finite and non-negative");
            }
            if i > 0 && zkm <= altitude_km[i - 1] {
                bail!("custom atmosphere altitude_km must increase with altitude");
            }
            s.z[i] = zkm * 1.0e5;
        }
    } else {
        s.z[0] = 0.0;
        for i in 1..n {
            // Hypsometric estimate in cm using dry-air gas constant / g.
            let tmean = 0.5 * (s.t[i - 1] + s.t[i]);
            let dz_m = 29.263 * tmean * (s.pstd[i - 1] / s.pstd[i]).ln();
            s.z[i] = s.z[i - 1] + dz_m * 100.0;
        }
    }

    s.do3int[n - 1] = s.do3ref[n - 1] * s.zzht;
    s.do2int[n - 1] = s.dm[n - 1] * s.zzht * po2;
    for j in (0..n - 1).rev() {
        let dz = (s.z[j + 1] - s.z[j]).abs();
        s.do3int[j] = s.do3int[j + 1] + 0.5 * dz * (s.do3ref[j + 1] + s.do3ref[j]);
        s.do2int[j] = s.do2int[j + 1] + 0.5 * dz * (s.dm[j + 1] + s.dm[j]) * po2;
    }

    for i in n..NATM {
        s.pstd[i] = 0.0;
        s.t[i] = 0.0;
        s.dm[i] = 0.0;
        s.theta[i] = 0.0;
        s.z[i] = s.z[n - 1];
        s.do3ref[i] = 0.0;
        s.refo3[i] = 0.0;
        s.do3int[i] = 0.0;
        s.do2int[i] = 0.0;
        s.asa[i] = -1.0;
        s.aer[i] = 0.0;
    }

    Ok(())
}

/// Add one chemistry-only atmospheric level and return its state index plus
/// the two radiative shells used to interpolate actinic flux.
fn append_interpolated_chemistry_level(
    s: &mut ModelState,
    altitude_km: f64,
) -> Result<(usize, usize, usize, f64)> {
    let target_z = altitude_km * 1.0e5;
    let exact_tolerance_cm = 1.0e-4;
    if let Some(level) = (0..s.nc).find(|&i| (s.z[i] - target_z).abs() <= exact_tolerance_cm) {
        return Ok((level, level, level, 0.0));
    }

    let upper = (0..s.nc).find(|&i| s.z[i] > target_z).ok_or_else(|| {
        anyhow::anyhow!("chemistry altitude {altitude_km} km is above the radiative grid")
    })?;
    if upper == 0 {
        bail!("chemistry altitude {altitude_km} km is below the radiative grid");
    }
    let lower = upper - 1;
    let dz = s.z[upper] - s.z[lower];
    if dz <= 0.0 {
        bail!("custom atmosphere altitude grid must be strictly increasing");
    }
    let weight = ((target_z - s.z[lower]) / dz).clamp(0.0, 1.0);

    if s.nlev >= NATM {
        bail!("too many off-grid chemistry levels; atmospheric capacity is {NATM}");
    }
    let level = s.nlev;
    s.nlev += 1;

    let interpolate =
        |lower_value: f64, upper_value: f64| lower_value + weight * (upper_value - lower_value);
    s.z[level] = target_z;
    s.pstd[level] = (interpolate(s.pstd[lower].ln(), s.pstd[upper].ln())).exp();
    s.t[level] = interpolate(s.t[lower], s.t[upper]);
    s.dm[level] = s.pstd[level] / (1.38e-19 * s.t[level]);
    s.theta[level] = s.t[level] * (1000.0 / s.pstd[level]).powf(0.2857);

    let lower_o3_mr = s.do3ref[lower] / s.dm[lower].max(1.0);
    let upper_o3_mr = s.do3ref[upper] / s.dm[upper].max(1.0);
    s.do3ref[level] = interpolate(lower_o3_mr, upper_o3_mr) * s.dm[level];
    s.refo3[level] = s.do3ref[level];
    s.do3int[level] = interpolate(s.do3int[lower], s.do3int[upper]);
    s.do2int[level] = interpolate(s.do2int[lower], s.do2int[upper]);

    if s.asa[lower] >= 0.0 && s.asa[upper] >= 0.0 {
        s.asa[level] = interpolate(s.asa[lower], s.asa[upper]);
        s.aer[level] = s.asa[level] * 0.25 * 2.7e-8;
    } else {
        s.asa[level] = -1.0;
        s.aer[level] = 0.0;
    }

    Ok((level, lower, upper, weight))
}

fn reconcile_custom_box_implicit_species(s: &mut ModelState, ib: usize) {
    let ialt = (s.nboxdo[ib].unsigned_abs() as usize)
        .saturating_sub(1)
        .min(NATM - 1);
    let mut xn = [0.0f64; NDEN];
    rplace(s, &mut xn, ib);
    let old_ialt = s.ialt;
    s.ialt = ialt;
    fixrat(&mut xn, s, ib);
    s.ialt = old_ialt;
    splace(s, &xn, ib);

    for j in 0..s.ntotx {
        for it in 0..s.ntimdo {
            s.xxnoft[[j, it, ib]] = xn[j];
        }
    }
    for j in 0..s.ntotx {
        s.xnoft[[j, 0]] = xn[j];
    }
}

fn disable_iodine(s: &mut ModelState, preserve_structure: bool) {
    if !preserve_structure {
        let old_o3_slot = s.n[10];
        s.liod = false;
        for j in 30..40 {
            s.n[j] = 0;
            s.ntsav[j] = 0;
        }
        // The augmented permutation stores the 29 active base radicals in
        // slots 1..29, iodine in 30..39, and explicit O3 in slot 40. Move O3
        // into slot 30 before compacting; truncating without this move makes
        // the chemistry operate on an uninitialized ozone slot.
        if old_o3_slot > 0 && old_o3_slot <= NDEN {
            let old = old_o3_slot - 1;
            let new = 29usize;
            for it in 0..s.xnoft.ncols() {
                s.xnoft[[new, it]] = s.xnoft[[old, it]];
            }
            for ib in 0..NB {
                for it in 0..s.xxnoft.shape()[1] {
                    s.xxnoft[[new, it, ib]] = s.xxnoft[[old, it, ib]];
                }
            }
        }
        s.n[10] = 30;
        s.ntsav[10] = 30;
        s.ntsav[NDEN] = 29;
        s.tnamet[29] = s.tname[10].clone();
        s.ntotx = 30;
        s.ntot = 29;
        let mut nnr = 0usize;
        for j in 0..s.nnr {
            let species = s.nnrt[j];
            if species <= 30 {
                s.nnrt[nnr] = species;
                nnr += 1;
            }
        }
        for j in nnr..s.nnrt.len() {
            s.nnrt[j] = 0;
        }
        s.nnr = nnr;
        for ib in 0..NB {
            s.fiodx[ib] = 0.0;
            s.di_[ib] = 0.0;
            s.dio[ib] = 0.0;
            s.dhoi[ib] = 0.0;
            s.diono2[ib] = 0.0;
            s.dhi[ib] = 0.0;
            s.doio[ib] = 0.0;
            s.di2[ib] = 0.0;
            s.di2o2[ib] = 0.0;
            s.di2o3[ib] = 0.0;
            s.di2o4[ib] = 0.0;
        }
        return;
    }

    // Keep the 40-species system structurally intact and use a negligible Iy
    // numerical floor. Removing the last ten rows makes the legacy Newton
    // system follow a different, cancellation-limited convergence path and is
    // not a reproducible iodine-off control.
    const IODINE_OFF_FLOOR: f64 = 1.0e-20;
    s.liod = true;
    let mut nnr = 0usize;
    for j in 0..s.nnr {
        let species = s.nnrt[j];
        if species <= 30 {
            s.nnrt[nnr] = species;
            nnr += 1;
        }
    }
    for j in nnr..s.nnrt.len() {
        s.nnrt[j] = 0;
    }
    s.nnr = nnr;
    for ib in 0..NB {
        s.fiodx[ib] = IODINE_OFF_FLOOR;
        let ialt = (s.nboxdo[ib].unsigned_abs() as usize)
            .saturating_sub(1)
            .min(NATM - 1);
        let one_iodine = IODINE_OFF_FLOOR * s.dm[ialt] / 10.0;
        s.di_[ib] = one_iodine;
        s.dio[ib] = one_iodine;
        s.dhoi[ib] = one_iodine;
        s.diono2[ib] = one_iodine;
        s.dhi[ib] = one_iodine;
        s.doio[ib] = one_iodine;
        s.di2[ib] = 0.5 * one_iodine;
        s.di2o2[ib] = 0.5 * one_iodine;
        s.di2o3[ib] = 0.5 * one_iodine;
        s.di2o4[ib] = 0.5 * one_iodine;
    }
}

fn apply_ctm_config(s: &mut ModelState, cfg: &CtmConfig) {
    let nbox = cfg.boxes.len().min(NB);
    s.nbox = nbox;

    // CTM nd216/nd216s are computed inside ctmlfq from boxin_gui.dat (jdaydo × lat).
    // Override the relevant boxin fields: ctmlfq reads s.xlat/xlatd and s.flscal/xdecd
    // only if boxin_gui.dat is absent; otherwise it reads the file.
    // We write the key fields directly so they take effect after boxin_gui.dat is parsed.
    s.xlatd = cfg.latitude_deg;
    s.xlat = cfg.latitude_deg.to_radians();

    let pi = std::f64::consts::PI;
    let xjd = 2.0 * pi * cfg.julian_day as f64 / 365.0;
    let decang = 6.918e-3 - 0.399912 * xjd.cos() + 0.070257 * xjd.sin()
        - 6.758e-3 * (2.0 * xjd).cos()
        + 9.07e-4 * (2.0 * xjd).sin()
        - 2.697e-3 * (3.0 * xjd).cos()
        + 1.480e-3 * (3.0 * xjd).sin();
    s.xdecd = decang * 57.29578;
    s.xdec = decang;

    const EDIST: [f64; 12] = [
        0.9837, 0.9875, 0.9945, 1.0032, 1.0109, 1.0158, 1.0165, 1.0128, 1.0057, 0.9970, 0.9892,
        0.9842,
    ];
    let mon = ((cfg.julian_day as i32 - 1) / 30).min(11) as usize;
    s.flscal = cfg.solar_flux_scale / (EDIST[mon] * EDIST[mon]);

    s.ndaysd = cfg.integration_days as i32;
    s.lbrom = cfg.bromine;

    for (ib, spec) in cfg.boxes.iter().take(nbox).enumerate() {
        s.nboxdo[ib] = spec.altitude_level as i32;
        let level = spec.altitude_level.saturating_sub(1) as usize;
        s.box_flux_lower[ib] = level;
        s.box_flux_upper[ib] = level;
        s.box_flux_upper_weight[ib] = 0.0;
        s.boxaa[ib] = spec.aerosol_surface_area_um2_cm3;
        s.boxss[ib] = spec.sea_salt_surface_area_um2_cm3;
        s.boxtt[ib] = spec.temp_offset_k;
        s.nboxmx[ib] = cfg.integration_days as i32;
        s.nboxwt[ib] = 1;
        s.nboxpr[ib] = 0;
        s.nboxct[ib] = 0;
        // Match the C++ default: rainout is disabled unless a future public
        // configuration field explicitly supplies it.
        s.boxrn[ib] = 0.0;
    }
    for ib in nbox..NB {
        s.nboxdo[ib] = 0;
        s.nboxwt[ib] = 0;
        s.box_flux_lower[ib] = 0;
        s.box_flux_upper[ib] = 0;
        s.box_flux_upper_weight[ib] = 0.0;
    }

    if !cfg.iodine {
        disable_iodine(s, true);
    }
}

// ── Output extraction ──────────────────────────────────────────────────────────

fn extract_box_snapshot(s: &ModelState, ib: usize) -> BoxSnapshot {
    let ialt = (s.nboxdo[ib].unsigned_abs() as usize)
        .saturating_sub(1)
        .min(NATM - 1);
    BoxSnapshot {
        box_index: ib,
        altitude_km: s.z[ialt] * 1e-5,
        pressure_mb: s.pstd[ialt],
        temperature_k: s.t[ialt] + s.boxtt[ib],
        air_density_cm3: s.dm[ialt],
        implicit: ImplicitSpecies::from_state(s, ib),
        long_lived: LongLivedMixingRatios::from_state(s, ib),
        jvalues: JValues::from_state_box_daily_mean(s, ib),
    }
}

fn extract_diurn_timeseries(s: &ModelState, ib: usize) -> DiurnBoxTimeSeries {
    let ialt = (s.nboxdo[ib].unsigned_abs() as usize)
        .saturating_sub(1)
        .min(NATM - 1);
    let ntimdo = s.ntimdo;

    let steps = (0..ntimdo)
        .map(|kt| DiurnTimeStep {
            elapsed_seconds: s.dtime[kt],
            time_hhmm: s.nhhmm[kt],
            implicit: ImplicitSpecies::from_timeseries(s, ib, kt),
        })
        .collect();

    DiurnBoxTimeSeries {
        box_index: ib,
        altitude_km: s.z[ialt] * 1e-5,
        pressure_mb: s.pstd[ialt],
        steps,
    }
}

fn extract_diurn_output(s: &ModelState) -> DiurnOutput {
    let nbox = s.nbox;
    let boxes: Vec<BoxSnapshot> = (0..nbox)
        .filter(|&ib| s.nboxdo[ib] != 0)
        .map(|ib| extract_box_snapshot(s, ib))
        .collect();

    let time_series: Vec<DiurnBoxTimeSeries> = (0..nbox)
        .filter(|&ib| s.nboxdo[ib] != 0)
        .map(|ib| extract_diurn_timeseries(s, ib))
        .collect();

    DiurnOutput {
        boxes,
        time_series,
        diagnostics: Diagnostics {
            raxloop: s.raxloop,
            radcount: s.radcount,
            newraf_nonconvergence_count: s.newraf_nonconvergence_count,
            rafday_nonconvergence_count: s.rafday_nonconvergence_count,
            rafday_max_final_relative_correction: s.rafday_max_final_relative_correction,
            rafday_max_correction_iterations: s.rafday_max_correction_iterations,
        },
    }
}

fn extract_ctm_output(s: &ModelState) -> CtmOutput {
    let nbox = s.nbox;
    let boxes: Vec<BoxSnapshot> = (0..nbox)
        .filter(|&ib| s.nboxdo[ib] != 0)
        .map(|ib| extract_box_snapshot(s, ib))
        .collect();

    CtmOutput {
        boxes,
        diagnostics: Diagnostics {
            raxloop: s.raxloop,
            radcount: s.radcount,
            newraf_nonconvergence_count: s.newraf_nonconvergence_count,
            rafday_nonconvergence_count: s.rafday_nonconvergence_count,
            rafday_max_final_relative_correction: s.rafday_max_final_relative_correction,
            rafday_max_correction_iterations: s.rafday_max_correction_iterations,
        },
    }
}

#[cfg(test)]
mod tests {
    use super::{
        extract_box_snapshot, time_distance_hhmm, validate_ctm_config, validate_diurn_config,
        validate_supported_mode, CtmBoxSpec, CtmConfig, CustomAtmosphereProfile, DiurnBoxSpec,
        DiurnConfig, JValues, LongLivedMixingRatios, O3InputKind,
    };
    use crate::state::ModelState;

    #[test]
    fn unsupported_derivs_mode_is_rejected() {
        let mut state = ModelState::new();
        state.nd216 = -1;
        let error = validate_supported_mode(&state).expect_err("DERIVS must be rejected");
        assert!(error.to_string().contains("DERIVS"));
    }

    #[test]
    fn unsupported_pzstd_mode_is_rejected() {
        let mut state = ModelState::new();
        state.npstd = 1;
        let error = validate_supported_mode(&state).expect_err("PZSTD must be rejected");
        assert!(error.to_string().contains("PZSTD"));
    }

    #[test]
    fn hhmm_distance_wraps_across_midnight() {
        assert_eq!(time_distance_hhmm(2359, 1), 2);
        assert_eq!(time_distance_hhmm(1200, 1200), 0);
        assert_eq!(time_distance_hhmm(160, 200), 0);
    }

    #[test]
    fn empty_box_configs_are_rejected() {
        let ctm_error = validate_ctm_config(&CtmConfig::default())
            .expect_err("an empty CTM box list must be rejected");
        assert!(ctm_error.to_string().contains("at least one box"));

        let diurn_error = validate_diurn_config(&DiurnConfig::default())
            .expect_err("DIURN needs boxes unless it has a custom atmosphere");
        assert!(diurn_error.to_string().contains("at least one box"));
    }

    #[test]
    fn invalid_box_and_mixing_ratio_values_are_rejected() {
        let ctm = CtmConfig {
            boxes: vec![CtmBoxSpec {
                altitude_level: 42,
                aerosol_surface_area_um2_cm3: 0.0,
                sea_salt_surface_area_um2_cm3: 0.0,
                temp_offset_k: 0.0,
            }],
            ..CtmConfig::default()
        };
        let ctm_error = validate_ctm_config(&ctm).expect_err("level 42 exceeds the model grid");
        assert!(ctm_error.to_string().contains("altitude_level"));

        let diurn = DiurnConfig {
            boxes: vec![DiurnBoxSpec {
                altitude_level: 1,
                altitude_km: None,
                aerosol_surface_area_um2_cm3: 0.0,
                sea_salt_surface_area_um2_cm3: 0.0,
                temp_offset_k: 0.0,
            }],
            initial_mixing_ratios: Some(vec![LongLivedMixingRatios {
                o3: f64::NAN,
                ..LongLivedMixingRatios::default()
            }]),
            ..DiurnConfig::default()
        };
        let diurn_error =
            validate_diurn_config(&diurn).expect_err("NaN mixing ratios must be rejected");
        assert!(diurn_error.to_string().contains(".o3"));

        let invalid_albedo = DiurnConfig {
            surface_albedo: 1.01,
            boxes: vec![DiurnBoxSpec {
                altitude_level: 1,
                altitude_km: None,
                aerosol_surface_area_um2_cm3: 0.0,
                sea_salt_surface_area_um2_cm3: 0.0,
                temp_offset_k: 0.0,
            }],
            ..DiurnConfig::default()
        };
        let albedo_error = validate_diurn_config(&invalid_albedo)
            .expect_err("surface albedo above one must be rejected");
        assert!(albedo_error.to_string().contains("surface_albedo"));
    }

    #[test]
    fn custom_diurn_accepts_cpp_81_level_radiative_grid() {
        let cfg = DiurnConfig {
            boxes: vec![DiurnBoxSpec {
                altitude_level: 81,
                altitude_km: None,
                aerosol_surface_area_um2_cm3: 0.0,
                sea_salt_surface_area_um2_cm3: 0.0,
                temp_offset_k: 0.0,
            }],
            atmosphere: Some(CustomAtmosphereProfile {
                pressure_mb: (0..81)
                    .map(|i| 1000.0 * (-0.125 * i as f64).exp())
                    .collect(),
                temperature_k: vec![230.0; 81],
                altitude_km: Some((0..81).map(|i| i as f64).collect()),
                o3: vec![5.0e-6; 81],
                o3_kind: O3InputKind::MixingRatio,
                aerosol_surface_area_um2_cm3: None,
            }),
            ..DiurnConfig::default()
        };
        validate_diurn_config(&cfg).expect("the C++ 0..80 km shell grid must be supported");
    }

    #[test]
    fn explicit_diurn_time_grid_is_validated() {
        let base = DiurnConfig {
            boxes: vec![DiurnBoxSpec {
                altitude_level: 1,
                altitude_km: None,
                aerosol_surface_area_um2_cm3: 0.0,
                sea_salt_surface_area_um2_cm3: 0.0,
                temp_offset_k: 0.0,
            }],
            ..DiurnConfig::default()
        };
        let half_hour = DiurnConfig {
            elapsed_time_hours: Some((0..=48).map(|index| index as f64 * 0.5).collect()),
            ..base.clone()
        };
        validate_diurn_config(&half_hour).expect("the archived 49-point grid must be valid");

        let duplicate = DiurnConfig {
            elapsed_time_hours: Some(vec![0.0, 12.0, 12.0, 24.0]),
            ..base
        };
        let error = validate_diurn_config(&duplicate)
            .expect_err("an explicit time grid must be strictly increasing");
        assert!(error.to_string().contains("strictly increasing"));
    }

    #[test]
    fn iodine_jvalue_fields_map_to_their_own_state_arrays() {
        let mut state = ModelState::new();
        state.vio[0] = 1.0;
        state.vhoi[0] = 2.0;
        state.viono2[0] = 3.0;
        state.voio[0] = 4.0;
        state.vi2[0] = 5.0;
        state.vi2o2[0] = 6.0;
        state.vi2o3[0] = 7.0;
        state.vi2o4[0] = 8.0;

        let values = JValues::from_state_alt(&state, 0);
        assert_eq!(
            [
                values.io,
                values.hoi,
                values.iono2,
                values.oio,
                values.i2,
                values.i2o2,
                values.i2o3,
                values.i2o4,
            ],
            [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
        );
    }

    #[test]
    fn box_snapshot_reports_daily_mean_jvalues() {
        let mut state = ModelState::new();
        state.nboxdo[0] = 15;
        state.ntimdo = 2;
        state.daysec = 86_400.0;
        state.dtime[0] = 0.0;
        state.dtime[1] = state.daysec;
        state.jtim[1] = 1;
        state.storjv[[8, 0, 0]] = 2.5;

        let snapshot = extract_box_snapshot(&state, 0);
        assert_eq!(snapshot.jvalues.no2, 2.5);
    }
}
