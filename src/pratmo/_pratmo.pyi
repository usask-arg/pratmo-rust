from __future__ import annotations

from typing import Optional

import numpy as np

class ImplicitSpecies:
    """Number densities (cm⁻³) for the 35 implicit Newton-Raphson species."""

    @property
    def no(self) -> float: ...
    @property
    def no2(self) -> float: ...
    @property
    def no3(self) -> float: ...
    @property
    def n2o5(self) -> float: ...
    @property
    def hno3(self) -> float: ...
    @property
    def h(self) -> float: ...
    @property
    def oh(self) -> float: ...
    @property
    def ho2(self) -> float: ...
    @property
    def h2o2(self) -> float: ...
    @property
    def o(self) -> float: ...
    @property
    def o3(self) -> float: ...
    @property
    def bro(self) -> float: ...
    @property
    def br(self) -> float: ...
    @property
    def hbr(self) -> float: ...
    @property
    def hno2(self) -> float: ...
    @property
    def hcl(self) -> float: ...
    @property
    def cl(self) -> float: ...
    @property
    def cl2(self) -> float: ...
    @property
    def clo(self) -> float: ...
    @property
    def clono2(self) -> float: ...
    @property
    def hno4(self) -> float: ...
    @property
    def hocl(self) -> float: ...
    @property
    def brono2(self) -> float: ...
    @property
    def hobr(self) -> float: ...
    @property
    def h2co(self) -> float: ...
    @property
    def ch3o2(self) -> float: ...
    @property
    def ch3o2h(self) -> float: ...
    @property
    def oclo(self) -> float: ...
    @property
    def cl2o2(self) -> float: ...
    @property
    def brcl(self) -> float: ...
    @property
    def i(self) -> float: ...
    @property
    def io(self) -> float: ...
    @property
    def hoi(self) -> float: ...
    @property
    def iono2(self) -> float: ...
    @property
    def hi(self) -> float: ...
    def to_dict(self) -> dict[str, float]:
        """Return all 35 species as a ``{name: value}`` dict (cm⁻³)."""
        ...

class LongLivedMixingRatios:
    """Dimensionless mixing ratios for the 18 long-lived species."""

    def __init__(
        self,
        *,
        o3: float = 0.0,
        n2o: float = 0.0,
        noy: float = 0.0,
        ch4: float = 0.0,
        co: float = 0.0,
        clx: float = 0.0,
        cf2cl2: float = 0.0,
        cfcl3: float = 0.0,
        ccl4: float = 0.0,
        ch3cl: float = 0.0,
        ch3ccl3: float = 0.0,
        h2: float = 0.0,
        h2o: float = 0.0,
        nh3: float = 0.0,
        c5h8: float = 0.0,
        brx: float = 0.0,
        ch3br: float = 0.0,
        ocs: float = 0.0,
        iodx: float = 0.0,
    ) -> None: ...
    o3: float
    n2o: float
    noy: float
    ch4: float
    co: float
    clx: float
    cf2cl2: float
    cfcl3: float
    ccl4: float
    ch3cl: float
    ch3ccl3: float
    h2: float
    h2o: float
    nh3: float
    c5h8: float
    brx: float
    ch3br: float
    ocs: float
    iodx: float
    def to_dict(self) -> dict[str, float]:
        """Return all 19 species as a ``{name: value}`` dict."""
        ...

class JValues:
    """Photolysis rates (s⁻¹) for all 47 J-value channels."""

    @property
    def no(self) -> float: ...
    @property
    def o2(self) -> float: ...
    @property
    def o3(self) -> float: ...
    @property
    def o3_o1d(self) -> float: ...
    @property
    def h2co_a(self) -> float: ...
    @property
    def h2co_b(self) -> float: ...
    @property
    def h2o2(self) -> float: ...
    @property
    def rooh(self) -> float: ...
    @property
    def no2(self) -> float: ...
    @property
    def no3_x(self) -> float: ...
    @property
    def no3_l(self) -> float: ...
    @property
    def n2o5(self) -> float: ...
    @property
    def hno2(self) -> float: ...
    @property
    def hno3(self) -> float: ...
    @property
    def hno4(self) -> float: ...
    @property
    def clono2(self) -> float: ...
    @property
    def cl2(self) -> float: ...
    @property
    def hocl(self) -> float: ...
    @property
    def oclo(self) -> float: ...
    @property
    def cl2o2(self) -> float: ...
    @property
    def clo(self) -> float: ...
    @property
    def bro(self) -> float: ...
    @property
    def brono2(self) -> float: ...
    @property
    def hobr(self) -> float: ...
    @property
    def n2o(self) -> float: ...
    @property
    def cfc11(self) -> float: ...
    @property
    def cfc12(self) -> float: ...
    @property
    def cfc113(self) -> float: ...
    @property
    def cfc114(self) -> float: ...
    @property
    def cfc115(self) -> float: ...
    @property
    def ccl4(self) -> float: ...
    @property
    def ch3cl(self) -> float: ...
    @property
    def ch3ccl3(self) -> float: ...
    @property
    def ch3br(self) -> float: ...
    @property
    def h1211(self) -> float: ...
    @property
    def h1301(self) -> float: ...
    @property
    def h2402(self) -> float: ...
    @property
    def hcfc22(self) -> float: ...
    @property
    def hcfc123(self) -> float: ...
    @property
    def hcfc141b(self) -> float: ...
    @property
    def chbr3(self) -> float: ...
    @property
    def ch3i(self) -> float: ...
    @property
    def cf3i(self) -> float: ...
    @property
    def ocs(self) -> float: ...
    @property
    def io(self) -> float: ...
    @property
    def hoi(self) -> float: ...
    @property
    def iono2(self) -> float: ...
    def to_dict(self) -> dict[str, float]:
        """Return all 47 J-values as a ``{name: value}`` dict (s⁻¹)."""
        ...

class Diagnostics:
    """Run diagnostics from a DIURN or CTM run."""

    @property
    def raxloop(self) -> float: ...
    @property
    def radcount(self) -> float: ...

class BoxSnapshot:
    """Snapshot of a single box's state (daily mean or final converged value)."""

    @property
    def box_index(self) -> int: ...
    @property
    def altitude_km(self) -> float: ...
    @property
    def pressure_mb(self) -> float: ...
    @property
    def temperature_k(self) -> float: ...
    @property
    def air_density_cm3(self) -> float: ...
    @property
    def implicit(self) -> ImplicitSpecies: ...
    @property
    def long_lived(self) -> LongLivedMixingRatios: ...
    @property
    def jvalues(self) -> JValues: ...

class DiurnTimeStep:
    """One time-step in a diurnal time series."""

    @property
    def time_hhmm(self) -> int:
        """Time in HHMM integer format (e.g. 1430 = 14:30 UTC)."""
        ...
    @property
    def implicit(self) -> ImplicitSpecies: ...

class DiurnBoxTimeSeries:
    """Full diurnal time series for one box."""

    @property
    def box_index(self) -> int: ...
    @property
    def altitude_km(self) -> float: ...
    @property
    def pressure_mb(self) -> float: ...
    @property
    def steps(self) -> list[DiurnTimeStep]: ...

class DiurnBoxSpec:
    """Per-box configuration for a DIURN run."""

    altitude_level: int
    """1-based standard pressure level index (1 = surface, 41 = top)."""
    albedo: float
    temp_offset_k: float
    def __init__(
        self,
        altitude_level: int,
        albedo: float = 0.0,
        temp_offset_k: float = 0.0,
    ) -> None: ...

class CtmBoxSpec:
    """Per-box configuration for a CTM run."""

    altitude_level: int
    """1-based standard pressure level index (1 = surface, 41 = top)."""
    albedo: float
    temp_offset_k: float
    def __init__(
        self,
        altitude_level: int,
        albedo: float = 0.0,
        temp_offset_k: float = 0.0,
    ) -> None: ...

class DiurnConfig:
    """Configuration for a diurnal cycle (DIURN) run."""

    latitude_deg: float
    julian_day: int
    integration_days: int
    boxes: list[DiurnBoxSpec]
    bromine: bool
    iodine: bool
    parallel_boxes: bool
    solar_flux_scale: float
    initial_mixing_ratios: Optional[list[LongLivedMixingRatios]]
    def __init__(
        self,
        *,
        latitude_deg: float = 0.0,
        julian_day: int = 120,
        integration_days: int = 20,
        boxes: list[DiurnBoxSpec] = [],
        bromine: bool = False,
        iodine: bool = True,
        parallel_boxes: bool = False,
        solar_flux_scale: float = 1.0,
        initial_mixing_ratios: Optional[list[LongLivedMixingRatios]] = None,
    ) -> None: ...

class CtmConfig:
    """Configuration for a CTM climatological run."""

    latitude_deg: float
    julian_day: int
    integration_days: int
    boxes: list[CtmBoxSpec]
    bromine: bool
    solar_flux_scale: float
    def __init__(
        self,
        *,
        latitude_deg: float = 60.0,
        julian_day: int = 75,
        integration_days: int = 40,
        boxes: list[CtmBoxSpec] = [],
        bromine: bool = False,
        solar_flux_scale: float = 1.0,
    ) -> None: ...

class DiurnOutput:
    """Output from a diurnal cycle run."""

    @property
    def boxes(self) -> list[BoxSnapshot]:
        """Daily-mean snapshot for each box."""
        ...
    @property
    def time_series(self) -> list[DiurnBoxTimeSeries]:
        """Full diurnal time series for each box."""
        ...
    @property
    def diagnostics(self) -> Diagnostics: ...
    def species_grid(self, species_name: str) -> np.ndarray:
        """Return an implicit species as a 2-D array of shape ``(n_boxes, n_timesteps)``.

        Parameters
        ----------
        species_name:
            One of: no, no2, no3, n2o5, hno3, h, oh, ho2, h2o2, o, o3, bro,
            br, hbr, hno2, hcl, cl, cl2, clo, clono2, hno4, hocl, brono2,
            hobr, h2co, ch3o2, ch3o2h, oclo, cl2o2, brcl
        """
        ...
    def jvalue_grid(self, jvalue_name: str) -> np.ndarray:
        """Return a J-value as a 2-D array of shape ``(n_boxes, n_timesteps)``.

        Parameters
        ----------
        jvalue_name:
            One of: no, o2, o3, o3_o1d, h2co_a, h2co_b, h2o2, rooh, no2,
            no3_x, no3_l, n2o5, hno2, hno3, hno4, clono2, cl2, hocl, oclo,
            cl2o2, clo, bro, brono2, hobr, n2o, cfc11, cfc12, cfc113, cfc114,
            cfc115, ccl4, ch3cl, ch3ccl3, ch3br, h1211, h1301, h2402, hcfc22,
            hcfc123, hcfc141b, chbr3, ch3i, cf3i, ocs
        """
        ...

class CtmOutput:
    """Output from a CTM climatological run."""

    @property
    def boxes(self) -> list[BoxSnapshot]:
        """Final converged snapshot for each box."""
        ...
    @property
    def diagnostics(self) -> Diagnostics: ...
    def species_profile(self, species_name: str) -> np.ndarray:
        """Return an implicit species as a 1-D array of shape ``(n_boxes,)``."""
        ...
    def jvalue_profile(self, jvalue_name: str) -> np.ndarray:
        """Return a J-value as a 1-D array of shape ``(n_boxes,)``."""
        ...

class PratmoModel:
    """Entry point for the PRATMO photochemical box model."""

    @staticmethod
    def with_defaults() -> PratmoModel:
        """Create a model using compiled-in embedded science data. No files needed."""
        ...
    @staticmethod
    def from_data_dir(data_dir: str) -> PratmoModel:
        """Create a model that loads science data from *data_dir* at runtime."""
        ...
    def run_diurn(self, cfg: DiurnConfig) -> DiurnOutput:
        """Run the diurnal cycle (DIURN + TPATH) mode."""
        ...
    def run_ctm(self, cfg: CtmConfig) -> CtmOutput:
        """Run the CTM climatological mode."""
        ...
