use std::collections::HashMap;
use std::path::PathBuf;

use ndarray::Array2;
use numpy::{IntoPyArray, PyArray1, PyArray2};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyTuple;

use pratmo_core::api::{
    BoxSnapshot, CtmBoxSpec, CtmConfig, CtmOutput, CustomAtmosphereProfile, Diagnostics,
    DiurnBoxSpec, DiurnBoxTimeSeries, DiurnConfig, DiurnOutput, DiurnTimeStep, ImplicitSpecies,
    JValues, LongLivedMixingRatios, No2ConstrainedDiurnConfig, No2ConstrainedDiurnOutput,
    O3InputKind, PratmoModel,
};

// ── String-dispatch helpers ────────────────────────────────────────────────────

const IMPLICIT_SPECIES_NAMES: &[&str] = &[
    "no", "no2", "no3", "n2o5", "hno3", "h", "oh", "ho2", "h2o2", "o", "o3", "bro", "br", "hbr",
    "hno2", "hcl", "cl", "cl2", "clo", "clono2", "hno4", "hocl", "brono2", "hobr", "h2co", "ch3o2",
    "ch3o2h", "oclo", "cl2o2", "brcl", "i", "io", "hoi", "iono2", "hi", "oio", "i2", "i2o2",
    "i2o3", "i2o4",
];

const LONG_LIVED_NAMES: &[&str] = &[
    "o3", "n2o", "noy", "ch4", "co", "clx", "cf2cl2", "cfcl3", "ccl4", "ch3cl", "ch3ccl3", "h2",
    "h2o", "nh3", "c5h8", "brx", "ch3br", "ocs", "iodx",
];

const JVALUE_NAMES: &[&str] = &[
    "no", "o2", "o3", "o3_o1d", "h2co_a", "h2co_b", "h2o2", "rooh", "no2", "no3_x", "no3_l",
    "n2o5", "hno2", "hno3", "hno4", "clono2", "cl2", "hocl", "oclo", "cl2o2", "clo", "bro",
    "brono2", "hobr", "n2o", "cfc11", "cfc12", "cfc113", "cfc114", "cfc115", "ccl4", "ch3cl",
    "ch3ccl3", "ch3br", "h1211", "h1301", "h2402", "hcfc22", "hcfc123", "hcfc141b", "chbr3",
    "ch3i", "cf3i", "ocs", "io", "hoi", "iono2", "oio", "i2", "i2o2", "i2o3", "i2o4",
];

fn normalize_field_name(name: &str) -> String {
    name.trim().to_ascii_lowercase()
}

fn implicit_field_by_name(s: &ImplicitSpecies, name: &str) -> Option<f64> {
    match name {
        "no" => Some(s.no),
        "no2" => Some(s.no2),
        "no3" => Some(s.no3),
        "n2o5" => Some(s.n2o5),
        "hno3" => Some(s.hno3),
        "h" => Some(s.h),
        "oh" => Some(s.oh),
        "ho2" => Some(s.ho2),
        "h2o2" => Some(s.h2o2),
        "o" => Some(s.o),
        "o3" => Some(s.o3),
        "bro" => Some(s.bro),
        "br" => Some(s.br),
        "hbr" => Some(s.hbr),
        "hno2" => Some(s.hno2),
        "hcl" => Some(s.hcl),
        "cl" => Some(s.cl),
        "cl2" => Some(s.cl2),
        "clo" => Some(s.clo),
        "clono2" => Some(s.clono2),
        "hno4" => Some(s.hno4),
        "hocl" => Some(s.hocl),
        "brono2" => Some(s.brono2),
        "hobr" => Some(s.hobr),
        "h2co" => Some(s.h2co),
        "ch3o2" => Some(s.ch3o2),
        "ch3o2h" => Some(s.ch3o2h),
        "oclo" => Some(s.oclo),
        "cl2o2" => Some(s.cl2o2),
        "brcl" => Some(s.brcl),
        "i" => Some(s.i),
        "io" => Some(s.io),
        "hoi" => Some(s.hoi),
        "iono2" => Some(s.iono2),
        "hi" => Some(s.hi),
        "oio" => Some(s.oio),
        "i2" => Some(s.i2),
        "i2o2" => Some(s.i2o2),
        "i2o3" => Some(s.i2o3),
        "i2o4" => Some(s.i2o4),
        _ => None,
    }
}

fn jvalue_field_by_name(j: &JValues, name: &str) -> Option<f64> {
    match name {
        "no" => Some(j.no),
        "o2" => Some(j.o2),
        "o3" => Some(j.o3),
        "o3_o1d" => Some(j.o3_o1d),
        "h2co_a" => Some(j.h2co_a),
        "h2co_b" => Some(j.h2co_b),
        "h2o2" => Some(j.h2o2),
        "rooh" => Some(j.rooh),
        "no2" => Some(j.no2),
        "no3_x" => Some(j.no3_x),
        "no3_l" => Some(j.no3_l),
        "n2o5" => Some(j.n2o5),
        "hno2" => Some(j.hno2),
        "hno3" => Some(j.hno3),
        "hno4" => Some(j.hno4),
        "clono2" => Some(j.clono2),
        "cl2" => Some(j.cl2),
        "hocl" => Some(j.hocl),
        "oclo" => Some(j.oclo),
        "cl2o2" => Some(j.cl2o2),
        "clo" => Some(j.clo),
        "bro" => Some(j.bro),
        "brono2" => Some(j.brono2),
        "hobr" => Some(j.hobr),
        "n2o" => Some(j.n2o),
        "cfc11" => Some(j.cfc11),
        "cfc12" => Some(j.cfc12),
        "cfc113" => Some(j.cfc113),
        "cfc114" => Some(j.cfc114),
        "cfc115" => Some(j.cfc115),
        "ccl4" => Some(j.ccl4),
        "ch3cl" => Some(j.ch3cl),
        "ch3ccl3" => Some(j.ch3ccl3),
        "ch3br" => Some(j.ch3br),
        "h1211" => Some(j.h1211),
        "h1301" => Some(j.h1301),
        "h2402" => Some(j.h2402),
        "hcfc22" => Some(j.hcfc22),
        "hcfc123" => Some(j.hcfc123),
        "hcfc141b" => Some(j.hcfc141b),
        "chbr3" => Some(j.chbr3),
        "ch3i" => Some(j.ch3i),
        "cf3i" => Some(j.cf3i),
        "ocs" => Some(j.ocs),
        "io" => Some(j.io),
        "hoi" => Some(j.hoi),
        "iono2" => Some(j.iono2),
        "oio" => Some(j.oio),
        "i2" => Some(j.i2),
        "i2o2" => Some(j.i2o2),
        "i2o3" => Some(j.i2o3),
        "i2o4" => Some(j.i2o4),
        _ => None,
    }
}

fn long_lived_field_by_name(ratios: &LongLivedMixingRatios, name: &str) -> Option<f64> {
    match name {
        "o3" => Some(ratios.o3),
        "n2o" => Some(ratios.n2o),
        "noy" => Some(ratios.noy),
        "ch4" => Some(ratios.ch4),
        "co" => Some(ratios.co),
        "clx" => Some(ratios.clx),
        "cf2cl2" => Some(ratios.cf2cl2),
        "cfcl3" => Some(ratios.cfcl3),
        "ccl4" => Some(ratios.ccl4),
        "ch3cl" => Some(ratios.ch3cl),
        "ch3ccl3" => Some(ratios.ch3ccl3),
        "h2" => Some(ratios.h2),
        "h2o" => Some(ratios.h2o),
        "nh3" => Some(ratios.nh3),
        "c5h8" => Some(ratios.c5h8),
        "brx" => Some(ratios.brx),
        "ch3br" => Some(ratios.ch3br),
        "ocs" => Some(ratios.ocs),
        "iodx" => Some(ratios.iodx),
        _ => None,
    }
}

fn unknown_field_error(kind: &str, requested: &str, valid_names: &[&str]) -> PyErr {
    PyValueError::new_err(format!(
        "Unknown {kind}: '{requested}'. Valid names: {}",
        valid_names.join(", ")
    ))
}

// ── ImplicitSpecies ────────────────────────────────────────────────────────────

#[pyclass(name = "ImplicitSpecies", frozen, from_py_object)]
#[derive(Clone)]
struct PyImplicitSpecies {
    inner: ImplicitSpecies,
}

#[pymethods]
impl PyImplicitSpecies {
    #[getter]
    fn no(&self) -> f64 {
        self.inner.no
    }
    #[getter]
    fn no2(&self) -> f64 {
        self.inner.no2
    }
    #[getter]
    fn no3(&self) -> f64 {
        self.inner.no3
    }
    #[getter]
    fn n2o5(&self) -> f64 {
        self.inner.n2o5
    }
    #[getter]
    fn hno3(&self) -> f64 {
        self.inner.hno3
    }
    #[getter]
    fn h(&self) -> f64 {
        self.inner.h
    }
    #[getter]
    fn oh(&self) -> f64 {
        self.inner.oh
    }
    #[getter]
    fn ho2(&self) -> f64 {
        self.inner.ho2
    }
    #[getter]
    fn h2o2(&self) -> f64 {
        self.inner.h2o2
    }
    #[getter]
    fn o(&self) -> f64 {
        self.inner.o
    }
    #[getter]
    fn o3(&self) -> f64 {
        self.inner.o3
    }
    #[getter]
    fn bro(&self) -> f64 {
        self.inner.bro
    }
    #[getter]
    fn br(&self) -> f64 {
        self.inner.br
    }
    #[getter]
    fn hbr(&self) -> f64 {
        self.inner.hbr
    }
    #[getter]
    fn hno2(&self) -> f64 {
        self.inner.hno2
    }
    #[getter]
    fn hcl(&self) -> f64 {
        self.inner.hcl
    }
    #[getter]
    fn cl(&self) -> f64 {
        self.inner.cl
    }
    #[getter]
    fn cl2(&self) -> f64 {
        self.inner.cl2
    }
    #[getter]
    fn clo(&self) -> f64 {
        self.inner.clo
    }
    #[getter]
    fn clono2(&self) -> f64 {
        self.inner.clono2
    }
    #[getter]
    fn hno4(&self) -> f64 {
        self.inner.hno4
    }
    #[getter]
    fn hocl(&self) -> f64 {
        self.inner.hocl
    }
    #[getter]
    fn brono2(&self) -> f64 {
        self.inner.brono2
    }
    #[getter]
    fn hobr(&self) -> f64 {
        self.inner.hobr
    }
    #[getter]
    fn h2co(&self) -> f64 {
        self.inner.h2co
    }
    #[getter]
    fn ch3o2(&self) -> f64 {
        self.inner.ch3o2
    }
    #[getter]
    fn ch3o2h(&self) -> f64 {
        self.inner.ch3o2h
    }
    #[getter]
    fn oclo(&self) -> f64 {
        self.inner.oclo
    }
    #[getter]
    fn cl2o2(&self) -> f64 {
        self.inner.cl2o2
    }
    #[getter]
    fn brcl(&self) -> f64 {
        self.inner.brcl
    }
    #[getter]
    fn i(&self) -> f64 {
        self.inner.i
    }
    #[getter]
    fn io(&self) -> f64 {
        self.inner.io
    }
    #[getter]
    fn hoi(&self) -> f64 {
        self.inner.hoi
    }
    #[getter]
    fn iono2(&self) -> f64 {
        self.inner.iono2
    }
    #[getter]
    fn hi(&self) -> f64 {
        self.inner.hi
    }
    #[getter]
    fn oio(&self) -> f64 {
        self.inner.oio
    }
    #[getter]
    fn i2(&self) -> f64 {
        self.inner.i2
    }
    #[getter]
    fn i2o2(&self) -> f64 {
        self.inner.i2o2
    }
    #[getter]
    fn i2o3(&self) -> f64 {
        self.inner.i2o3
    }
    #[getter]
    fn i2o4(&self) -> f64 {
        self.inner.i2o4
    }

    /// Return all 40 species as a ``{name: value}`` dict (cm⁻³).
    fn to_dict(&self) -> HashMap<&'static str, f64> {
        [
            ("no", self.inner.no),
            ("no2", self.inner.no2),
            ("no3", self.inner.no3),
            ("n2o5", self.inner.n2o5),
            ("hno3", self.inner.hno3),
            ("h", self.inner.h),
            ("oh", self.inner.oh),
            ("ho2", self.inner.ho2),
            ("h2o2", self.inner.h2o2),
            ("o", self.inner.o),
            ("o3", self.inner.o3),
            ("bro", self.inner.bro),
            ("br", self.inner.br),
            ("hbr", self.inner.hbr),
            ("hno2", self.inner.hno2),
            ("hcl", self.inner.hcl),
            ("cl", self.inner.cl),
            ("cl2", self.inner.cl2),
            ("clo", self.inner.clo),
            ("clono2", self.inner.clono2),
            ("hno4", self.inner.hno4),
            ("hocl", self.inner.hocl),
            ("brono2", self.inner.brono2),
            ("hobr", self.inner.hobr),
            ("h2co", self.inner.h2co),
            ("ch3o2", self.inner.ch3o2),
            ("ch3o2h", self.inner.ch3o2h),
            ("oclo", self.inner.oclo),
            ("cl2o2", self.inner.cl2o2),
            ("brcl", self.inner.brcl),
            ("i", self.inner.i),
            ("io", self.inner.io),
            ("hoi", self.inner.hoi),
            ("iono2", self.inner.iono2),
            ("hi", self.inner.hi),
            ("oio", self.inner.oio),
            ("i2", self.inner.i2),
            ("i2o2", self.inner.i2o2),
            ("i2o3", self.inner.i2o3),
            ("i2o4", self.inner.i2o4),
        ]
        .into_iter()
        .collect()
    }

    fn __repr__(&self) -> String {
        format!(
            "ImplicitSpecies(o3={:.3e}, oh={:.3e}, no={:.3e}, no2={:.3e})",
            self.inner.o3, self.inner.oh, self.inner.no, self.inner.no2
        )
    }
}

// ── LongLivedMixingRatios ──────────────────────────────────────────────────────

#[pyclass(name = "LongLivedMixingRatios", from_py_object)]
#[derive(Clone)]
struct PyLongLivedMixingRatios {
    inner: LongLivedMixingRatios,
}

#[pymethods]
impl PyLongLivedMixingRatios {
    /// Construct long-lived mixing ratios. All arguments are keyword-only and default to 0.0.
    #[new]
    #[pyo3(signature = (
        o3=0.0, n2o=0.0, noy=0.0, ch4=0.0, co=0.0,
        clx=0.0, cf2cl2=0.0, cfcl3=0.0, ccl4=0.0, ch3cl=0.0,
        ch3ccl3=0.0, h2=0.0, h2o=0.0, nh3=0.0, c5h8=0.0,
        brx=0.0, ch3br=0.0, ocs=0.0, iodx=0.0
    ))]
    #[allow(clippy::too_many_arguments)]
    fn new(
        o3: f64,
        n2o: f64,
        noy: f64,
        ch4: f64,
        co: f64,
        clx: f64,
        cf2cl2: f64,
        cfcl3: f64,
        ccl4: f64,
        ch3cl: f64,
        ch3ccl3: f64,
        h2: f64,
        h2o: f64,
        nh3: f64,
        c5h8: f64,
        brx: f64,
        ch3br: f64,
        ocs: f64,
        iodx: f64,
    ) -> Self {
        Self {
            inner: LongLivedMixingRatios {
                o3,
                n2o,
                noy,
                ch4,
                co,
                clx,
                cf2cl2,
                cfcl3,
                ccl4,
                ch3cl,
                ch3ccl3,
                h2,
                h2o,
                nh3,
                c5h8,
                brx,
                ch3br,
                ocs,
                iodx,
            },
        }
    }

    #[getter]
    fn o3(&self) -> f64 {
        self.inner.o3
    }
    #[setter]
    fn set_o3(&mut self, v: f64) {
        self.inner.o3 = v;
    }

    #[getter]
    fn n2o(&self) -> f64 {
        self.inner.n2o
    }
    #[setter]
    fn set_n2o(&mut self, v: f64) {
        self.inner.n2o = v;
    }

    #[getter]
    fn noy(&self) -> f64 {
        self.inner.noy
    }
    #[setter]
    fn set_noy(&mut self, v: f64) {
        self.inner.noy = v;
    }

    #[getter]
    fn ch4(&self) -> f64 {
        self.inner.ch4
    }
    #[setter]
    fn set_ch4(&mut self, v: f64) {
        self.inner.ch4 = v;
    }

    #[getter]
    fn co(&self) -> f64 {
        self.inner.co
    }
    #[setter]
    fn set_co(&mut self, v: f64) {
        self.inner.co = v;
    }

    #[getter]
    fn clx(&self) -> f64 {
        self.inner.clx
    }
    #[setter]
    fn set_clx(&mut self, v: f64) {
        self.inner.clx = v;
    }

    #[getter]
    fn cf2cl2(&self) -> f64 {
        self.inner.cf2cl2
    }
    #[setter]
    fn set_cf2cl2(&mut self, v: f64) {
        self.inner.cf2cl2 = v;
    }

    #[getter]
    fn cfcl3(&self) -> f64 {
        self.inner.cfcl3
    }
    #[setter]
    fn set_cfcl3(&mut self, v: f64) {
        self.inner.cfcl3 = v;
    }

    #[getter]
    fn ccl4(&self) -> f64 {
        self.inner.ccl4
    }
    #[setter]
    fn set_ccl4(&mut self, v: f64) {
        self.inner.ccl4 = v;
    }

    #[getter]
    fn ch3cl(&self) -> f64 {
        self.inner.ch3cl
    }
    #[setter]
    fn set_ch3cl(&mut self, v: f64) {
        self.inner.ch3cl = v;
    }

    #[getter]
    fn ch3ccl3(&self) -> f64 {
        self.inner.ch3ccl3
    }
    #[setter]
    fn set_ch3ccl3(&mut self, v: f64) {
        self.inner.ch3ccl3 = v;
    }

    #[getter]
    fn h2(&self) -> f64 {
        self.inner.h2
    }
    #[setter]
    fn set_h2(&mut self, v: f64) {
        self.inner.h2 = v;
    }

    #[getter]
    fn h2o(&self) -> f64 {
        self.inner.h2o
    }
    #[setter]
    fn set_h2o(&mut self, v: f64) {
        self.inner.h2o = v;
    }

    #[getter]
    fn nh3(&self) -> f64 {
        self.inner.nh3
    }
    #[setter]
    fn set_nh3(&mut self, v: f64) {
        self.inner.nh3 = v;
    }

    #[getter]
    fn c5h8(&self) -> f64 {
        self.inner.c5h8
    }
    #[setter]
    fn set_c5h8(&mut self, v: f64) {
        self.inner.c5h8 = v;
    }

    #[getter]
    fn brx(&self) -> f64 {
        self.inner.brx
    }
    #[setter]
    fn set_brx(&mut self, v: f64) {
        self.inner.brx = v;
    }

    #[getter]
    fn ch3br(&self) -> f64 {
        self.inner.ch3br
    }
    #[setter]
    fn set_ch3br(&mut self, v: f64) {
        self.inner.ch3br = v;
    }

    #[getter]
    fn ocs(&self) -> f64 {
        self.inner.ocs
    }
    #[setter]
    fn set_ocs(&mut self, v: f64) {
        self.inner.ocs = v;
    }

    #[getter]
    fn iodx(&self) -> f64 {
        self.inner.iodx
    }
    #[setter]
    fn set_iodx(&mut self, v: f64) {
        self.inner.iodx = v;
    }

    /// Return all 19 species as a ``{name: value}`` dict (dimensionless mixing ratios).
    fn to_dict(&self) -> HashMap<&'static str, f64> {
        [
            ("o3", self.inner.o3),
            ("n2o", self.inner.n2o),
            ("noy", self.inner.noy),
            ("ch4", self.inner.ch4),
            ("co", self.inner.co),
            ("clx", self.inner.clx),
            ("cf2cl2", self.inner.cf2cl2),
            ("cfcl3", self.inner.cfcl3),
            ("ccl4", self.inner.ccl4),
            ("ch3cl", self.inner.ch3cl),
            ("ch3ccl3", self.inner.ch3ccl3),
            ("h2", self.inner.h2),
            ("h2o", self.inner.h2o),
            ("nh3", self.inner.nh3),
            ("c5h8", self.inner.c5h8),
            ("brx", self.inner.brx),
            ("ch3br", self.inner.ch3br),
            ("ocs", self.inner.ocs),
            ("iodx", self.inner.iodx),
        ]
        .into_iter()
        .collect()
    }

    fn __repr__(&self) -> String {
        format!(
            "LongLivedMixingRatios(o3={:.3e}, n2o={:.3e}, ch4={:.3e})",
            self.inner.o3, self.inner.n2o, self.inner.ch4
        )
    }
}

// ── JValues ────────────────────────────────────────────────────────────────────

#[pyclass(name = "JValues", frozen, from_py_object)]
#[derive(Clone)]
struct PyJValues {
    inner: JValues,
}

#[pymethods]
impl PyJValues {
    #[getter]
    fn no(&self) -> f64 {
        self.inner.no
    }
    #[getter]
    fn o2(&self) -> f64 {
        self.inner.o2
    }
    #[getter]
    fn o3(&self) -> f64 {
        self.inner.o3
    }
    #[getter]
    fn o3_o1d(&self) -> f64 {
        self.inner.o3_o1d
    }
    #[getter]
    fn h2co_a(&self) -> f64 {
        self.inner.h2co_a
    }
    #[getter]
    fn h2co_b(&self) -> f64 {
        self.inner.h2co_b
    }
    #[getter]
    fn h2o2(&self) -> f64 {
        self.inner.h2o2
    }
    #[getter]
    fn rooh(&self) -> f64 {
        self.inner.rooh
    }
    #[getter]
    fn no2(&self) -> f64 {
        self.inner.no2
    }
    #[getter]
    fn no3_x(&self) -> f64 {
        self.inner.no3_x
    }
    #[getter]
    fn no3_l(&self) -> f64 {
        self.inner.no3_l
    }
    #[getter]
    fn n2o5(&self) -> f64 {
        self.inner.n2o5
    }
    #[getter]
    fn hno2(&self) -> f64 {
        self.inner.hno2
    }
    #[getter]
    fn hno3(&self) -> f64 {
        self.inner.hno3
    }
    #[getter]
    fn hno4(&self) -> f64 {
        self.inner.hno4
    }
    #[getter]
    fn clono2(&self) -> f64 {
        self.inner.clono2
    }
    #[getter]
    fn cl2(&self) -> f64 {
        self.inner.cl2
    }
    #[getter]
    fn hocl(&self) -> f64 {
        self.inner.hocl
    }
    #[getter]
    fn oclo(&self) -> f64 {
        self.inner.oclo
    }
    #[getter]
    fn cl2o2(&self) -> f64 {
        self.inner.cl2o2
    }
    #[getter]
    fn clo(&self) -> f64 {
        self.inner.clo
    }
    #[getter]
    fn bro(&self) -> f64 {
        self.inner.bro
    }
    #[getter]
    fn brono2(&self) -> f64 {
        self.inner.brono2
    }
    #[getter]
    fn hobr(&self) -> f64 {
        self.inner.hobr
    }
    #[getter]
    fn n2o(&self) -> f64 {
        self.inner.n2o
    }
    #[getter]
    fn cfc11(&self) -> f64 {
        self.inner.cfc11
    }
    #[getter]
    fn cfc12(&self) -> f64 {
        self.inner.cfc12
    }
    #[getter]
    fn cfc113(&self) -> f64 {
        self.inner.cfc113
    }
    #[getter]
    fn cfc114(&self) -> f64 {
        self.inner.cfc114
    }
    #[getter]
    fn cfc115(&self) -> f64 {
        self.inner.cfc115
    }
    #[getter]
    fn ccl4(&self) -> f64 {
        self.inner.ccl4
    }
    #[getter]
    fn ch3cl(&self) -> f64 {
        self.inner.ch3cl
    }
    #[getter]
    fn ch3ccl3(&self) -> f64 {
        self.inner.ch3ccl3
    }
    #[getter]
    fn ch3br(&self) -> f64 {
        self.inner.ch3br
    }
    #[getter]
    fn h1211(&self) -> f64 {
        self.inner.h1211
    }
    #[getter]
    fn h1301(&self) -> f64 {
        self.inner.h1301
    }
    #[getter]
    fn h2402(&self) -> f64 {
        self.inner.h2402
    }
    #[getter]
    fn hcfc22(&self) -> f64 {
        self.inner.hcfc22
    }
    #[getter]
    fn hcfc123(&self) -> f64 {
        self.inner.hcfc123
    }
    #[getter]
    fn hcfc141b(&self) -> f64 {
        self.inner.hcfc141b
    }
    #[getter]
    fn chbr3(&self) -> f64 {
        self.inner.chbr3
    }
    #[getter]
    fn ch3i(&self) -> f64 {
        self.inner.ch3i
    }
    #[getter]
    fn cf3i(&self) -> f64 {
        self.inner.cf3i
    }
    #[getter]
    fn ocs(&self) -> f64 {
        self.inner.ocs
    }
    #[getter]
    fn io(&self) -> f64 {
        self.inner.io
    }
    #[getter]
    fn hoi(&self) -> f64 {
        self.inner.hoi
    }
    #[getter]
    fn iono2(&self) -> f64 {
        self.inner.iono2
    }
    #[getter]
    fn oio(&self) -> f64 {
        self.inner.oio
    }
    #[getter]
    fn i2(&self) -> f64 {
        self.inner.i2
    }
    #[getter]
    fn i2o2(&self) -> f64 {
        self.inner.i2o2
    }
    #[getter]
    fn i2o3(&self) -> f64 {
        self.inner.i2o3
    }
    #[getter]
    fn i2o4(&self) -> f64 {
        self.inner.i2o4
    }

    /// Return all 52 J-values as a ``{name: value}`` dict (s⁻¹).
    fn to_dict(&self) -> HashMap<&'static str, f64> {
        [
            ("no", self.inner.no),
            ("o2", self.inner.o2),
            ("o3", self.inner.o3),
            ("o3_o1d", self.inner.o3_o1d),
            ("h2co_a", self.inner.h2co_a),
            ("h2co_b", self.inner.h2co_b),
            ("h2o2", self.inner.h2o2),
            ("rooh", self.inner.rooh),
            ("no2", self.inner.no2),
            ("no3_x", self.inner.no3_x),
            ("no3_l", self.inner.no3_l),
            ("n2o5", self.inner.n2o5),
            ("hno2", self.inner.hno2),
            ("hno3", self.inner.hno3),
            ("hno4", self.inner.hno4),
            ("clono2", self.inner.clono2),
            ("cl2", self.inner.cl2),
            ("hocl", self.inner.hocl),
            ("oclo", self.inner.oclo),
            ("cl2o2", self.inner.cl2o2),
            ("clo", self.inner.clo),
            ("bro", self.inner.bro),
            ("brono2", self.inner.brono2),
            ("hobr", self.inner.hobr),
            ("n2o", self.inner.n2o),
            ("cfc11", self.inner.cfc11),
            ("cfc12", self.inner.cfc12),
            ("cfc113", self.inner.cfc113),
            ("cfc114", self.inner.cfc114),
            ("cfc115", self.inner.cfc115),
            ("ccl4", self.inner.ccl4),
            ("ch3cl", self.inner.ch3cl),
            ("ch3ccl3", self.inner.ch3ccl3),
            ("ch3br", self.inner.ch3br),
            ("h1211", self.inner.h1211),
            ("h1301", self.inner.h1301),
            ("h2402", self.inner.h2402),
            ("hcfc22", self.inner.hcfc22),
            ("hcfc123", self.inner.hcfc123),
            ("hcfc141b", self.inner.hcfc141b),
            ("chbr3", self.inner.chbr3),
            ("ch3i", self.inner.ch3i),
            ("cf3i", self.inner.cf3i),
            ("ocs", self.inner.ocs),
            ("io", self.inner.io),
            ("hoi", self.inner.hoi),
            ("iono2", self.inner.iono2),
            ("oio", self.inner.oio),
            ("i2", self.inner.i2),
            ("i2o2", self.inner.i2o2),
            ("i2o3", self.inner.i2o3),
            ("i2o4", self.inner.i2o4),
        ]
        .into_iter()
        .collect()
    }

    fn __repr__(&self) -> String {
        format!(
            "JValues(no2={:.3e}, o3={:.3e}, no={:.3e})",
            self.inner.no2, self.inner.o3, self.inner.no
        )
    }
}

// ── Diagnostics ────────────────────────────────────────────────────────────────

#[pyclass(name = "Diagnostics", frozen, from_py_object)]
#[derive(Clone)]
struct PyDiagnostics {
    inner: Diagnostics,
}

#[pymethods]
impl PyDiagnostics {
    #[getter]
    fn raxloop(&self) -> f64 {
        self.inner.raxloop
    }
    #[getter]
    fn radcount(&self) -> f64 {
        self.inner.radcount
    }
    #[getter]
    fn newraf_nonconvergence_count(&self) -> usize {
        self.inner.newraf_nonconvergence_count
    }
    #[getter]
    fn rafday_nonconvergence_count(&self) -> usize {
        self.inner.rafday_nonconvergence_count
    }

    #[getter]
    fn rafday_max_final_relative_correction(&self) -> f64 {
        self.inner.rafday_max_final_relative_correction
    }

    #[getter]
    fn rafday_max_correction_iterations(&self) -> usize {
        self.inner.rafday_max_correction_iterations
    }

    fn __repr__(&self) -> String {
        format!(
            "Diagnostics(raxloop={}, radcount={}, newraf_nonconvergence_count={}, rafday_nonconvergence_count={}, rafday_max_final_relative_correction={}, rafday_max_correction_iterations={})",
            self.inner.raxloop,
            self.inner.radcount,
            self.inner.newraf_nonconvergence_count,
            self.inner.rafday_nonconvergence_count,
            self.inner.rafday_max_final_relative_correction,
            self.inner.rafday_max_correction_iterations
        )
    }
}

// ── BoxSnapshot ────────────────────────────────────────────────────────────────

#[pyclass(name = "BoxSnapshot", frozen, from_py_object)]
#[derive(Clone)]
struct PyBoxSnapshot {
    inner: BoxSnapshot,
}

#[pymethods]
impl PyBoxSnapshot {
    #[getter]
    fn box_index(&self) -> usize {
        self.inner.box_index
    }
    #[getter]
    fn altitude_km(&self) -> f64 {
        self.inner.altitude_km
    }
    #[getter]
    fn pressure_mb(&self) -> f64 {
        self.inner.pressure_mb
    }
    #[getter]
    fn temperature_k(&self) -> f64 {
        self.inner.temperature_k
    }
    #[getter]
    fn air_density_cm3(&self) -> f64 {
        self.inner.air_density_cm3
    }

    #[getter]
    fn implicit(&self) -> PyImplicitSpecies {
        PyImplicitSpecies {
            inner: self.inner.implicit.clone(),
        }
    }
    #[getter]
    fn long_lived(&self) -> PyLongLivedMixingRatios {
        PyLongLivedMixingRatios {
            inner: self.inner.long_lived.clone(),
        }
    }
    #[getter]
    fn jvalues(&self) -> PyJValues {
        PyJValues {
            inner: self.inner.jvalues.clone(),
        }
    }

    fn __repr__(&self) -> String {
        format!(
            "BoxSnapshot(box={}, alt={:.1}km, p={:.1}mb, T={:.1}K)",
            self.inner.box_index,
            self.inner.altitude_km,
            self.inner.pressure_mb,
            self.inner.temperature_k
        )
    }
}

// ── DiurnTimeStep ──────────────────────────────────────────────────────────────

#[pyclass(name = "DiurnTimeStep", frozen, from_py_object)]
#[derive(Clone)]
struct PyDiurnTimeStep {
    inner: DiurnTimeStep,
}

#[pymethods]
impl PyDiurnTimeStep {
    /// Monotonic seconds since the noon start of the 24-hour orbit.
    #[getter]
    fn elapsed_seconds(&self) -> f64 {
        self.inner.elapsed_seconds
    }

    /// Local clock label. Both orbit endpoints are noon (1200).
    #[getter]
    fn time_hhmm(&self) -> i32 {
        self.inner.time_hhmm
    }

    #[getter]
    fn implicit(&self) -> PyImplicitSpecies {
        PyImplicitSpecies {
            inner: self.inner.implicit.clone(),
        }
    }

    fn __repr__(&self) -> String {
        format!(
            "DiurnTimeStep(elapsed_seconds={}, time_hhmm={:04})",
            self.inner.elapsed_seconds, self.inner.time_hhmm
        )
    }
}

// ── DiurnBoxTimeSeries ─────────────────────────────────────────────────────────

#[pyclass(name = "DiurnBoxTimeSeries", frozen, from_py_object)]
#[derive(Clone)]
struct PyDiurnBoxTimeSeries {
    inner: DiurnBoxTimeSeries,
}

#[pymethods]
impl PyDiurnBoxTimeSeries {
    fn __len__(&self) -> usize {
        self.inner.steps.len()
    }

    #[getter]
    fn box_index(&self) -> usize {
        self.inner.box_index
    }
    #[getter]
    fn altitude_km(&self) -> f64 {
        self.inner.altitude_km
    }
    #[getter]
    fn pressure_mb(&self) -> f64 {
        self.inner.pressure_mb
    }

    #[getter]
    fn steps(&self) -> Vec<PyDiurnTimeStep> {
        self.inner
            .steps
            .iter()
            .map(|s| PyDiurnTimeStep { inner: s.clone() })
            .collect()
    }

    fn __repr__(&self) -> String {
        format!(
            "DiurnBoxTimeSeries(box={}, alt={:.1}km, {} steps)",
            self.inner.box_index,
            self.inner.altitude_km,
            self.inner.steps.len()
        )
    }
}

// ── Config: DiurnBoxSpec / CtmBoxSpec ──────────────────────────────────────────

/// Per-box configuration for a DIURN run. ``altitude_level`` is a 1-based pressure level index.
#[pyclass(name = "DiurnBoxSpec", from_py_object)]
#[derive(Clone)]
struct PyDiurnBoxSpec {
    #[pyo3(get, set)]
    altitude_level: u8,
    #[pyo3(get, set)]
    aerosol_surface_area_um2_cm3: f64,
    #[pyo3(get, set)]
    sea_salt_surface_area_um2_cm3: f64,
    #[pyo3(get, set)]
    temp_offset_k: f64,
}

#[pymethods]
impl PyDiurnBoxSpec {
    #[new]
    #[pyo3(signature = (altitude_level, aerosol_surface_area_um2_cm3=0.0, sea_salt_surface_area_um2_cm3=0.0, temp_offset_k=0.0))]
    fn new(
        altitude_level: u8,
        aerosol_surface_area_um2_cm3: f64,
        sea_salt_surface_area_um2_cm3: f64,
        temp_offset_k: f64,
    ) -> Self {
        Self {
            altitude_level,
            aerosol_surface_area_um2_cm3,
            sea_salt_surface_area_um2_cm3,
            temp_offset_k,
        }
    }

    fn __repr__(&self) -> String {
        format!(
            "DiurnBoxSpec(altitude_level={}, aerosol_surface_area_um2_cm3={}, sea_salt_surface_area_um2_cm3={}, temp_offset_k={})",
            self.altitude_level,
            self.aerosol_surface_area_um2_cm3,
            self.sea_salt_surface_area_um2_cm3,
            self.temp_offset_k
        )
    }
}

/// Per-box configuration for a CTM run. ``altitude_level`` is a 1-based pressure level index.
#[pyclass(name = "CtmBoxSpec", from_py_object)]
#[derive(Clone)]
struct PyCtmBoxSpec {
    #[pyo3(get, set)]
    altitude_level: u8,
    #[pyo3(get, set)]
    aerosol_surface_area_um2_cm3: f64,
    #[pyo3(get, set)]
    sea_salt_surface_area_um2_cm3: f64,
    #[pyo3(get, set)]
    temp_offset_k: f64,
}

#[pymethods]
impl PyCtmBoxSpec {
    #[new]
    #[pyo3(signature = (altitude_level, aerosol_surface_area_um2_cm3=0.0, sea_salt_surface_area_um2_cm3=0.0, temp_offset_k=0.0))]
    fn new(
        altitude_level: u8,
        aerosol_surface_area_um2_cm3: f64,
        sea_salt_surface_area_um2_cm3: f64,
        temp_offset_k: f64,
    ) -> Self {
        Self {
            altitude_level,
            aerosol_surface_area_um2_cm3,
            sea_salt_surface_area_um2_cm3,
            temp_offset_k,
        }
    }

    fn __repr__(&self) -> String {
        format!(
            "CtmBoxSpec(altitude_level={}, aerosol_surface_area_um2_cm3={}, sea_salt_surface_area_um2_cm3={}, temp_offset_k={})",
            self.altitude_level,
            self.aerosol_surface_area_um2_cm3,
            self.sea_salt_surface_area_um2_cm3,
            self.temp_offset_k
        )
    }
}

#[pyclass(name = "CustomAtmosphereProfile", from_py_object)]
#[derive(Clone)]
struct PyCustomAtmosphereProfile {
    #[pyo3(get, set)]
    pressure_mb: Vec<f64>,
    #[pyo3(get, set)]
    temperature_k: Vec<f64>,
    #[pyo3(get, set)]
    altitude_km: Option<Vec<f64>>,
    #[pyo3(get, set)]
    o3: Vec<f64>,
    #[pyo3(get, set)]
    o3_kind: String,
}

#[pymethods]
impl PyCustomAtmosphereProfile {
    #[new]
    #[pyo3(signature = (pressure_mb, temperature_k, o3, o3_kind="mixing_ratio".to_string(), altitude_km=None))]
    fn new(
        pressure_mb: Vec<f64>,
        temperature_k: Vec<f64>,
        o3: Vec<f64>,
        o3_kind: String,
        altitude_km: Option<Vec<f64>>,
    ) -> Self {
        Self {
            pressure_mb,
            temperature_k,
            altitude_km,
            o3,
            o3_kind,
        }
    }

    fn __repr__(&self) -> String {
        format!(
            "CustomAtmosphereProfile({} levels, o3_kind='{}')",
            self.pressure_mb.len(),
            self.o3_kind
        )
    }
}

impl PyCustomAtmosphereProfile {
    fn to_rust(&self) -> PyResult<CustomAtmosphereProfile> {
        let o3_kind = match self.o3_kind.as_str() {
            "mixing_ratio" | "vmr" => O3InputKind::MixingRatio,
            "number_density" | "density" | "cm-3" | "cm3" => O3InputKind::NumberDensityCm3,
            other => {
                return Err(PyValueError::new_err(format!(
                    "Unknown o3_kind '{other}'. Use 'mixing_ratio' or 'number_density'."
                )));
            }
        };
        Ok(CustomAtmosphereProfile {
            pressure_mb: self.pressure_mb.clone(),
            temperature_k: self.temperature_k.clone(),
            altitude_km: self.altitude_km.clone(),
            o3: self.o3.clone(),
            o3_kind,
        })
    }
}

// ── DiurnConfig ────────────────────────────────────────────────────────────────

/// Configuration for a diurnal cycle (DIURN) run.
#[pyclass(name = "DiurnConfig", from_py_object)]
#[derive(Clone)]
struct PyDiurnConfig {
    #[pyo3(get, set)]
    latitude_deg: f64,
    #[pyo3(get, set)]
    julian_day: u16,
    #[pyo3(get, set)]
    integration_days: u32,
    #[pyo3(get, set)]
    boxes: Vec<PyDiurnBoxSpec>,
    #[pyo3(get, set)]
    bromine: bool,
    #[pyo3(get, set)]
    iodine: bool,
    #[pyo3(get, set)]
    parallel_boxes: bool,
    #[pyo3(get, set)]
    solar_flux_scale: f64,
    #[pyo3(get, set)]
    atmosphere: Option<PyCustomAtmosphereProfile>,
    #[pyo3(get, set)]
    initial_mixing_ratios: Option<Vec<PyLongLivedMixingRatios>>,
}

#[pymethods]
impl PyDiurnConfig {
    #[new]
    #[pyo3(signature = (
        latitude_deg=0.0,
        julian_day=120,
        integration_days=20,
        boxes=vec![],
        bromine=false,
        iodine=true,
        parallel_boxes=false,
        solar_flux_scale=1.0,
        atmosphere=None,
        initial_mixing_ratios=None
    ))]
    #[allow(clippy::too_many_arguments)]
    fn new(
        latitude_deg: f64,
        julian_day: u16,
        integration_days: u32,
        boxes: Vec<PyDiurnBoxSpec>,
        bromine: bool,
        iodine: bool,
        parallel_boxes: bool,
        solar_flux_scale: f64,
        atmosphere: Option<PyCustomAtmosphereProfile>,
        initial_mixing_ratios: Option<Vec<PyLongLivedMixingRatios>>,
    ) -> Self {
        Self {
            latitude_deg,
            julian_day,
            integration_days,
            boxes,
            bromine,
            iodine,
            parallel_boxes,
            solar_flux_scale,
            atmosphere,
            initial_mixing_ratios,
        }
    }

    fn __repr__(&self) -> String {
        format!(
            "DiurnConfig(latitude_deg={}, julian_day={}, integration_days={}, boxes={}, bromine={}, iodine={}, parallel_boxes={})",
            self.latitude_deg,
            self.julian_day,
            self.integration_days,
            self.boxes.len(),
            self.bromine,
            self.iodine,
            self.parallel_boxes
        )
    }
}

impl PyDiurnConfig {
    fn to_rust(&self) -> PyResult<DiurnConfig> {
        Ok(DiurnConfig {
            latitude_deg: self.latitude_deg,
            julian_day: self.julian_day,
            integration_days: self.integration_days,
            boxes: self
                .boxes
                .iter()
                .map(|b| DiurnBoxSpec {
                    altitude_level: b.altitude_level,
                    aerosol_surface_area_um2_cm3: b.aerosol_surface_area_um2_cm3,
                    sea_salt_surface_area_um2_cm3: b.sea_salt_surface_area_um2_cm3,
                    temp_offset_k: b.temp_offset_k,
                })
                .collect(),
            bromine: self.bromine,
            iodine: self.iodine,
            parallel_boxes: self.parallel_boxes,
            solar_flux_scale: self.solar_flux_scale,
            atmosphere: self.atmosphere.as_ref().map(|a| a.to_rust()).transpose()?,
            initial_mixing_ratios: self
                .initial_mixing_ratios
                .as_ref()
                .map(|v| v.iter().map(|mr| mr.inner.clone()).collect()),
        })
    }
}

// ── CtmConfig ──────────────────────────────────────────────────────────────────

/// Configuration for a CTM climatological run.
#[pyclass(name = "CtmConfig")]
struct PyCtmConfig {
    #[pyo3(get, set)]
    latitude_deg: f64,
    #[pyo3(get, set)]
    julian_day: u16,
    #[pyo3(get, set)]
    integration_days: u32,
    #[pyo3(get, set)]
    boxes: Vec<PyCtmBoxSpec>,
    #[pyo3(get, set)]
    bromine: bool,
    #[pyo3(get, set)]
    iodine: bool,
    #[pyo3(get, set)]
    solar_flux_scale: f64,
}

#[pymethods]
impl PyCtmConfig {
    #[new]
    #[pyo3(signature = (
        latitude_deg=60.0,
        julian_day=75,
        integration_days=40,
        boxes=vec![],
        bromine=false,
        iodine=true,
        solar_flux_scale=1.0
    ))]
    fn new(
        latitude_deg: f64,
        julian_day: u16,
        integration_days: u32,
        boxes: Vec<PyCtmBoxSpec>,
        bromine: bool,
        iodine: bool,
        solar_flux_scale: f64,
    ) -> Self {
        Self {
            latitude_deg,
            julian_day,
            integration_days,
            boxes,
            bromine,
            iodine,
            solar_flux_scale,
        }
    }

    fn __repr__(&self) -> String {
        format!(
            "CtmConfig(latitude_deg={}, julian_day={}, integration_days={}, boxes={}, bromine={}, iodine={})",
            self.latitude_deg,
            self.julian_day,
            self.integration_days,
            self.boxes.len(),
            self.bromine,
            self.iodine
        )
    }
}

impl PyCtmConfig {
    fn to_rust(&self) -> CtmConfig {
        CtmConfig {
            latitude_deg: self.latitude_deg,
            julian_day: self.julian_day,
            integration_days: self.integration_days,
            boxes: self
                .boxes
                .iter()
                .map(|b| CtmBoxSpec {
                    altitude_level: b.altitude_level,
                    aerosol_surface_area_um2_cm3: b.aerosol_surface_area_um2_cm3,
                    sea_salt_surface_area_um2_cm3: b.sea_salt_surface_area_um2_cm3,
                    temp_offset_k: b.temp_offset_k,
                })
                .collect(),
            bromine: self.bromine,
            iodine: self.iodine,
            solar_flux_scale: self.solar_flux_scale,
        }
    }
}

// ── DiurnOutput ────────────────────────────────────────────────────────────────

#[pyclass(name = "DiurnOutput", frozen)]
struct PyDiurnOutput {
    inner: DiurnOutput,
}

#[pymethods]
impl PyDiurnOutput {
    fn __len__(&self) -> usize {
        self.inner.boxes.len()
    }

    /// Daily-mean snapshot for each box.
    #[getter]
    fn boxes(&self) -> Vec<PyBoxSnapshot> {
        self.inner
            .boxes
            .iter()
            .map(|b| PyBoxSnapshot { inner: b.clone() })
            .collect()
    }

    /// Full diurnal time series for each box.
    #[getter]
    fn time_series(&self) -> Vec<PyDiurnBoxTimeSeries> {
        self.inner
            .time_series
            .iter()
            .map(|ts| PyDiurnBoxTimeSeries { inner: ts.clone() })
            .collect()
    }

    #[getter]
    fn diagnostics(&self) -> PyDiagnostics {
        PyDiagnostics {
            inner: self.inner.diagnostics.clone(),
        }
    }

    /// Box altitudes in km, with shape ``(n_boxes,)``.
    #[getter]
    fn altitude_km<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        PyArray1::from_vec(py, self.inner.boxes.iter().map(|b| b.altitude_km).collect())
    }

    /// Box pressures in mb, with shape ``(n_boxes,)``.
    #[getter]
    fn pressure_mb<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        PyArray1::from_vec(py, self.inner.boxes.iter().map(|b| b.pressure_mb).collect())
    }

    /// Box temperatures in K, with shape ``(n_boxes,)``.
    #[getter]
    fn temperature_k<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        PyArray1::from_vec(
            py,
            self.inner.boxes.iter().map(|b| b.temperature_k).collect(),
        )
    }

    /// Box air number densities in cm⁻³, with shape ``(n_boxes,)``.
    #[getter]
    fn air_density_cm3<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        PyArray1::from_vec(
            py,
            self.inner.boxes.iter().map(|b| b.air_density_cm3).collect(),
        )
    }

    /// Shared DIURN coordinate in seconds since the noon start, shape ``(n_timesteps,)``.
    #[getter]
    fn elapsed_seconds<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        let values = self
            .inner
            .time_series
            .first()
            .map(|ts| ts.steps.iter().map(|step| step.elapsed_seconds).collect())
            .unwrap_or_default();
        PyArray1::from_vec(py, values)
    }

    /// Shared cyclic local-time labels in HHMM form, shape ``(n_timesteps,)``.
    #[getter]
    fn time_hhmm<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<i32>> {
        let values = self
            .inner
            .time_series
            .first()
            .map(|ts| ts.steps.iter().map(|step| step.time_hhmm).collect())
            .unwrap_or_default();
        PyArray1::from_vec(py, values)
    }

    /// Return a daily-mean implicit species profile with shape ``(n_boxes,)``.
    fn species_profile<'py>(
        &self,
        py: Python<'py>,
        species_name: &str,
    ) -> PyResult<Bound<'py, PyArray1<f64>>> {
        let name = normalize_field_name(species_name);
        if implicit_field_by_name(&ImplicitSpecies::default(), &name).is_none() {
            return Err(unknown_field_error(
                "implicit species",
                species_name,
                IMPLICIT_SPECIES_NAMES,
            ));
        }
        let values = self
            .inner
            .boxes
            .iter()
            .map(|box_| implicit_field_by_name(&box_.implicit, &name).unwrap())
            .collect();
        Ok(PyArray1::from_vec(py, values))
    }

    /// Return a daily-mean long-lived mixing-ratio profile with shape ``(n_boxes,)``.
    fn long_lived_profile<'py>(
        &self,
        py: Python<'py>,
        species_name: &str,
    ) -> PyResult<Bound<'py, PyArray1<f64>>> {
        let name = normalize_field_name(species_name);
        if long_lived_field_by_name(&LongLivedMixingRatios::default(), &name).is_none() {
            return Err(unknown_field_error(
                "long-lived species",
                species_name,
                LONG_LIVED_NAMES,
            ));
        }
        let values = self
            .inner
            .boxes
            .iter()
            .map(|box_| long_lived_field_by_name(&box_.long_lived, &name).unwrap())
            .collect();
        Ok(PyArray1::from_vec(py, values))
    }

    /// Return a daily-mean J-value profile with shape ``(n_boxes,)``.
    fn jvalue_profile<'py>(
        &self,
        py: Python<'py>,
        jvalue_name: &str,
    ) -> PyResult<Bound<'py, PyArray1<f64>>> {
        let name = normalize_field_name(jvalue_name);
        if jvalue_field_by_name(&JValues::default(), &name).is_none() {
            return Err(unknown_field_error("J-value", jvalue_name, JVALUE_NAMES));
        }
        let values = self
            .inner
            .boxes
            .iter()
            .map(|box_| jvalue_field_by_name(&box_.jvalues, &name).unwrap())
            .collect();
        Ok(PyArray1::from_vec(py, values))
    }

    /// Return an implicit species as a 2-D numpy array of shape ``(n_boxes, n_timesteps)``.
    ///
    /// Parameters
    /// ----------
    /// species_name:
    ///     A name from ``pratmo.IMPLICIT_SPECIES_NAMES`` (case-insensitive).
    fn species_grid<'py>(
        &self,
        py: Python<'py>,
        species_name: &str,
    ) -> PyResult<Bound<'py, PyArray2<f64>>> {
        let name = normalize_field_name(species_name);
        if implicit_field_by_name(&ImplicitSpecies::default(), &name).is_none() {
            return Err(unknown_field_error(
                "implicit species",
                species_name,
                IMPLICIT_SPECIES_NAMES,
            ));
        }
        let arr: Array2<f64> = self
            .inner
            .species_grid(|s| implicit_field_by_name(s, &name).unwrap());
        Ok(arr.into_pyarray(py))
    }

    /// Return daily-mean J-values repeated into a 2-D array of shape
    /// ``(n_boxes, n_timesteps)``.
    ///
    /// Parameters
    /// ----------
    /// jvalue_name:
    ///     A name from ``pratmo.JVALUE_NAMES`` (case-insensitive).
    fn jvalue_grid<'py>(
        &self,
        py: Python<'py>,
        jvalue_name: &str,
    ) -> PyResult<Bound<'py, PyArray2<f64>>> {
        let name = normalize_field_name(jvalue_name);
        if jvalue_field_by_name(&JValues::default(), &name).is_none() {
            return Err(unknown_field_error("J-value", jvalue_name, JVALUE_NAMES));
        }
        let arr: Array2<f64> = self
            .inner
            .jvalue_grid(|j| jvalue_field_by_name(j, &name).unwrap());
        Ok(arr.into_pyarray(py))
    }

    fn __repr__(&self) -> String {
        format!(
            "DiurnOutput({} boxes, {} time series)",
            self.inner.boxes.len(),
            self.inner.time_series.len()
        )
    }
}

#[pyclass(name = "No2ConstrainedDiurnConfig")]
struct PyNo2ConstrainedDiurnConfig {
    #[pyo3(get, set)]
    diurn: PyDiurnConfig,
    #[pyo3(get, set)]
    observed_no2_cm3: Vec<f64>,
    #[pyo3(get, set)]
    target_hhmm: i32,
    #[pyo3(get, set)]
    iterations: usize,
}

#[pymethods]
impl PyNo2ConstrainedDiurnConfig {
    #[new]
    #[pyo3(signature = (diurn, observed_no2_cm3, target_hhmm, iterations=3))]
    fn new(
        diurn: PyDiurnConfig,
        observed_no2_cm3: Vec<f64>,
        target_hhmm: i32,
        iterations: usize,
    ) -> Self {
        Self {
            diurn,
            observed_no2_cm3,
            target_hhmm,
            iterations,
        }
    }

    fn __repr__(&self) -> String {
        format!(
            "No2ConstrainedDiurnConfig(boxes={}, target_hhmm={:04}, iterations={})",
            self.observed_no2_cm3.len(),
            self.target_hhmm,
            self.iterations
        )
    }
}

impl PyNo2ConstrainedDiurnConfig {
    fn to_rust(&self) -> PyResult<No2ConstrainedDiurnConfig> {
        Ok(No2ConstrainedDiurnConfig {
            diurn: self.diurn.to_rust()?,
            observed_no2_cm3: self.observed_no2_cm3.clone(),
            target_hhmm: self.target_hhmm,
            iterations: self.iterations,
        })
    }
}

#[pyclass(name = "No2ConstrainedDiurnOutput")]
struct PyNo2ConstrainedDiurnOutput {
    inner: No2ConstrainedDiurnOutput,
}

#[pymethods]
impl PyNo2ConstrainedDiurnOutput {
    #[getter]
    fn output(&self) -> PyDiurnOutput {
        PyDiurnOutput {
            inner: self.inner.output.clone(),
        }
    }

    #[getter]
    fn noy_scale(&self) -> Vec<f64> {
        self.inner.noy_scale.clone()
    }

    #[getter]
    fn modeled_no2_cm3(&self) -> Vec<f64> {
        self.inner.modeled_no2_cm3.clone()
    }

    fn __repr__(&self) -> String {
        format!(
            "No2ConstrainedDiurnOutput({} boxes)",
            self.inner.noy_scale.len()
        )
    }
}

// ── CtmOutput ──────────────────────────────────────────────────────────────────

#[pyclass(name = "CtmOutput", frozen)]
struct PyCtmOutput {
    inner: CtmOutput,
}

#[pymethods]
impl PyCtmOutput {
    fn __len__(&self) -> usize {
        self.inner.boxes.len()
    }

    /// Final converged snapshot for each box.
    #[getter]
    fn boxes(&self) -> Vec<PyBoxSnapshot> {
        self.inner
            .boxes
            .iter()
            .map(|b| PyBoxSnapshot { inner: b.clone() })
            .collect()
    }

    #[getter]
    fn diagnostics(&self) -> PyDiagnostics {
        PyDiagnostics {
            inner: self.inner.diagnostics.clone(),
        }
    }

    /// Box altitudes in km, with shape ``(n_boxes,)``.
    #[getter]
    fn altitude_km<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        PyArray1::from_vec(py, self.inner.boxes.iter().map(|b| b.altitude_km).collect())
    }

    /// Box pressures in mb, with shape ``(n_boxes,)``.
    #[getter]
    fn pressure_mb<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        PyArray1::from_vec(py, self.inner.boxes.iter().map(|b| b.pressure_mb).collect())
    }

    /// Box temperatures in K, with shape ``(n_boxes,)``.
    #[getter]
    fn temperature_k<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        PyArray1::from_vec(
            py,
            self.inner.boxes.iter().map(|b| b.temperature_k).collect(),
        )
    }

    /// Box air number densities in cm⁻³, with shape ``(n_boxes,)``.
    #[getter]
    fn air_density_cm3<'py>(&self, py: Python<'py>) -> Bound<'py, PyArray1<f64>> {
        PyArray1::from_vec(
            py,
            self.inner.boxes.iter().map(|b| b.air_density_cm3).collect(),
        )
    }

    /// Return an implicit species as a 1-D numpy array of shape ``(n_boxes,)``.
    ///
    /// Parameters
    /// ----------
    /// species_name:
    ///     See :meth:`DiurnOutput.species_grid` for valid names.
    fn species_profile<'py>(
        &self,
        py: Python<'py>,
        species_name: &str,
    ) -> PyResult<Bound<'py, PyArray1<f64>>> {
        let name = normalize_field_name(species_name);
        if implicit_field_by_name(&ImplicitSpecies::default(), &name).is_none() {
            return Err(unknown_field_error(
                "implicit species",
                species_name,
                IMPLICIT_SPECIES_NAMES,
            ));
        }
        let v: Vec<f64> = self
            .inner
            .species_profile(|s| implicit_field_by_name(s, &name).unwrap());
        Ok(PyArray1::from_vec(py, v))
    }

    /// Return a long-lived mixing ratio as a 1-D array of shape ``(n_boxes,)``.
    fn long_lived_profile<'py>(
        &self,
        py: Python<'py>,
        species_name: &str,
    ) -> PyResult<Bound<'py, PyArray1<f64>>> {
        let name = normalize_field_name(species_name);
        if long_lived_field_by_name(&LongLivedMixingRatios::default(), &name).is_none() {
            return Err(unknown_field_error(
                "long-lived species",
                species_name,
                LONG_LIVED_NAMES,
            ));
        }
        let values = self
            .inner
            .boxes
            .iter()
            .map(|box_| long_lived_field_by_name(&box_.long_lived, &name).unwrap())
            .collect();
        Ok(PyArray1::from_vec(py, values))
    }

    /// Return a J-value as a 1-D numpy array of shape ``(n_boxes,)``.
    ///
    /// Parameters
    /// ----------
    /// jvalue_name:
    ///     See :meth:`DiurnOutput.jvalue_grid` for valid names.
    fn jvalue_profile<'py>(
        &self,
        py: Python<'py>,
        jvalue_name: &str,
    ) -> PyResult<Bound<'py, PyArray1<f64>>> {
        let name = normalize_field_name(jvalue_name);
        if jvalue_field_by_name(&JValues::default(), &name).is_none() {
            return Err(unknown_field_error("J-value", jvalue_name, JVALUE_NAMES));
        }
        let v: Vec<f64> = self
            .inner
            .jvalue_profile(|j| jvalue_field_by_name(j, &name).unwrap());
        Ok(PyArray1::from_vec(py, v))
    }

    fn __repr__(&self) -> String {
        format!("CtmOutput({} boxes)", self.inner.boxes.len())
    }
}

// ── PratmoModel ────────────────────────────────────────────────────────────────

/// Entry point for the PRATMO photochemical box model.
///
/// Examples
/// --------
/// >>> model = PratmoModel.with_defaults()
/// >>> cfg = DiurnConfig(latitude_deg=0.0, julian_day=120, integration_days=5,
/// ...                   boxes=[DiurnBoxSpec(altitude_level=25)])
/// >>> out = model.run_diurn(cfg)
/// >>> o3 = out.species_grid("o3")   # shape (1, n_timesteps)
#[pyclass(name = "PratmoModel")]
struct PyPratmoModel {
    inner: PratmoModel,
}

#[pymethods]
impl PyPratmoModel {
    /// Create a model using only compiled-in embedded science data. No files needed.
    #[staticmethod]
    fn with_defaults() -> Self {
        Self {
            inner: PratmoModel::with_defaults(),
        }
    }

    /// Create a model that reads science data from *data_dir* at runtime.
    ///
    /// The directory must contain ``fort01.x``, ``fort10_cam06.x``,
    /// ``fort11_jpl09.x``, ``fort13.x``, ``fort14.x``, ``J_H2O_SZA0.dat``.
    /// CTM runs additionally need ``fort03_LLM.x``, ``fort05.x``, ``fort51.x``.
    #[staticmethod]
    fn from_data_dir(data_dir: PathBuf) -> Self {
        Self {
            inner: PratmoModel::new(data_dir),
        }
    }

    /// Run the diurnal cycle (DIURN + TPATH) mode.
    fn run_diurn(&self, cfg: &PyDiurnConfig) -> PyResult<PyDiurnOutput> {
        let rust_cfg = cfg.to_rust()?;
        self.inner
            .run_diurn(&rust_cfg)
            .map(|out| PyDiurnOutput { inner: out })
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    fn run_diurn_no2_constrained(
        &self,
        cfg: &PyNo2ConstrainedDiurnConfig,
    ) -> PyResult<PyNo2ConstrainedDiurnOutput> {
        let rust_cfg = cfg.to_rust()?;
        self.inner
            .run_diurn_no2_constrained(&rust_cfg)
            .map(|out| PyNo2ConstrainedDiurnOutput { inner: out })
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    /// Run the CTM climatological mode.
    fn run_ctm(&self, cfg: &PyCtmConfig) -> PyResult<PyCtmOutput> {
        self.inner
            .run_ctm(&cfg.to_rust())
            .map(|out| PyCtmOutput { inner: out })
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    fn __repr__(&self) -> &'static str {
        "PratmoModel()"
    }
}

// ── Module ─────────────────────────────────────────────────────────────────────

#[pymodule]
fn _pratmo(m: &Bound<'_, PyModule>) -> PyResult<()> {
    let py = m.py();
    m.add(
        "IMPLICIT_SPECIES_NAMES",
        PyTuple::new(py, IMPLICIT_SPECIES_NAMES.iter().copied())?,
    )?;
    m.add(
        "LONG_LIVED_NAMES",
        PyTuple::new(py, LONG_LIVED_NAMES.iter().copied())?,
    )?;
    m.add(
        "JVALUE_NAMES",
        PyTuple::new(py, JVALUE_NAMES.iter().copied())?,
    )?;
    m.add_class::<PyPratmoModel>()?;
    m.add_class::<PyDiurnConfig>()?;
    m.add_class::<PyCtmConfig>()?;
    m.add_class::<PyDiurnBoxSpec>()?;
    m.add_class::<PyCtmBoxSpec>()?;
    m.add_class::<PyCustomAtmosphereProfile>()?;
    m.add_class::<PyNo2ConstrainedDiurnConfig>()?;
    m.add_class::<PyNo2ConstrainedDiurnOutput>()?;
    m.add_class::<PyDiurnOutput>()?;
    m.add_class::<PyCtmOutput>()?;
    m.add_class::<PyBoxSnapshot>()?;
    m.add_class::<PyDiurnBoxTimeSeries>()?;
    m.add_class::<PyDiurnTimeStep>()?;
    m.add_class::<PyImplicitSpecies>()?;
    m.add_class::<PyLongLivedMixingRatios>()?;
    m.add_class::<PyJValues>()?;
    m.add_class::<PyDiagnostics>()?;
    Ok(())
}
