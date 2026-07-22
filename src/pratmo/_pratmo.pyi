from __future__ import annotations

from typing import Final, Optional

import numpy as np

IMPLICIT_SPECIES_NAMES: Final[tuple[str, ...]]
"""Names accepted by implicit-species profile and grid methods."""

LONG_LIVED_NAMES: Final[tuple[str, ...]]
"""Names accepted by long-lived mixing-ratio profile methods."""

JVALUE_NAMES: Final[tuple[str, ...]]
"""Names accepted by J-value profile and grid methods."""

class ImplicitSpecies:
    """Number densities (cm⁻³) for the 40 implicit Newton-Raphson species."""

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
    @property
    def oio(self) -> float: ...
    @property
    def i2(self) -> float: ...
    @property
    def i2o2(self) -> float: ...
    @property
    def i2o3(self) -> float: ...
    @property
    def i2o4(self) -> float: ...
    def to_dict(self) -> dict[str, float]:
        """Return all 40 species as a ``{name: value}`` dict (cm⁻³)."""
        ...

class LongLivedMixingRatios:
    """Dimensionless mixing ratios for the 19 long-lived species/families.

    Defaults form a representative, non-singular lower-stratospheric state.
    Use climatological or measured values for quantitative work.
    """

    def __init__(
        self,
        *,
        o3: float = 5.0e-6,
        n2o: float = 300.0e-9,
        noy: float = 10.0e-9,
        ch4: float = 1.7e-6,
        co: float = 50.0e-9,
        clx: float = 3.0e-9,
        cf2cl2: float = 500.0e-12,
        cfcl3: float = 250.0e-12,
        ccl4: float = 90.0e-12,
        ch3cl: float = 550.0e-12,
        ch3ccl3: float = 20.0e-12,
        h2: float = 500.0e-9,
        h2o: float = 5.0e-6,
        nh3: float = 0.0,
        c5h8: float = 0.0,
        brx: float = 20.0e-12,
        ch3br: float = 10.0e-12,
        ocs: float = 500.0e-12,
        iodx: float = 1.0e-12,
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
    """Photolysis rates (s⁻¹) for all 52 J-value channels."""

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
    @property
    def oio(self) -> float: ...
    @property
    def i2(self) -> float: ...
    @property
    def i2o2(self) -> float: ...
    @property
    def i2o3(self) -> float: ...
    @property
    def i2o4(self) -> float: ...
    def to_dict(self) -> dict[str, float]:
        """Return all 52 J-values as a ``{name: value}`` dict (s⁻¹)."""
        ...

class Diagnostics:
    """Run diagnostics from a DIURN or CTM run."""

    @property
    def raxloop(self) -> float: ...
    @property
    def radcount(self) -> float: ...
    @property
    def newraf_nonconvergence_count(self) -> int: ...
    @property
    def rafday_nonconvergence_count(self) -> int: ...
    @property
    def rafday_max_final_relative_correction(self) -> float: ...
    @property
    def rafday_max_correction_iterations(self) -> int: ...

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
    def elapsed_seconds(self) -> float:
        """Monotonic seconds since the noon start of the 24-hour orbit."""
        ...
    @property
    def time_hhmm(self) -> int:
        """Local clock label; both orbit endpoints are noon (1200)."""
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
    def __len__(self) -> int: ...

class DiurnBoxSpec:
    """Per-box configuration for a DIURN run."""

    altitude_level: int
    """1-based level index (legacy grid: 1..41; custom radiative grid: up to 81)."""
    altitude_km: Optional[float]
    """Exact chemistry altitude; flux is interpolated between radiative shells."""
    aerosol_surface_area_um2_cm3: float
    sea_salt_surface_area_um2_cm3: float
    temp_offset_k: float
    def __init__(
        self,
        altitude_level: int,
        aerosol_surface_area_um2_cm3: float = 0.25,
        sea_salt_surface_area_um2_cm3: float = 0.0,
        temp_offset_k: float = 0.0,
        altitude_km: Optional[float] = None,
    ) -> None: ...

class CtmBoxSpec:
    """Per-box configuration for a CTM run."""

    altitude_level: int
    """1-based standard pressure level index (1 = surface, 41 = top)."""
    aerosol_surface_area_um2_cm3: float
    sea_salt_surface_area_um2_cm3: float
    temp_offset_k: float
    def __init__(
        self,
        altitude_level: int,
        aerosol_surface_area_um2_cm3: float = 0.0,
        sea_salt_surface_area_um2_cm3: float = 0.0,
        temp_offset_k: float = 0.0,
    ) -> None: ...

class CustomAtmosphereProfile:
    """Custom pressure/temperature/O3 grid for DIURN runs."""

    pressure_mb: list[float]
    temperature_k: list[float]
    altitude_km: Optional[list[float]]
    o3: list[float]
    o3_kind: str
    aerosol_surface_area_um2_cm3: Optional[list[float]]
    def __init__(
        self,
        pressure_mb: list[float],
        temperature_k: list[float],
        o3: list[float],
        o3_kind: str = "mixing_ratio",
        altitude_km: Optional[list[float]] = None,
        aerosol_surface_area_um2_cm3: Optional[list[float]] = None,
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
    cpp_compatibility: bool
    elapsed_time_hours: Optional[list[float]]
    solar_flux_scale: float
    surface_albedo: float
    heterogeneous_chemistry: bool
    radiative_aerosol: bool
    atmosphere: Optional[CustomAtmosphereProfile]
    initial_mixing_ratios: Optional[list[LongLivedMixingRatios]]
    def __init__(
        self,
        *,
        latitude_deg: float = 0.0,
        julian_day: int = 120,
        integration_days: int = 20,
        boxes: Optional[list[DiurnBoxSpec]] = None,
        bromine: bool = True,
        iodine: bool = False,
        parallel_boxes: bool = False,
        cpp_compatibility: bool = False,
        elapsed_time_hours: Optional[list[float]] = None,
        solar_flux_scale: float = 1.0,
        surface_albedo: float = 0.20,
        heterogeneous_chemistry: bool = True,
        radiative_aerosol: bool = False,
        atmosphere: Optional[CustomAtmosphereProfile] = None,
        initial_mixing_ratios: Optional[list[LongLivedMixingRatios]] = None,
    ) -> None: ...

class CtmConfig:
    """Configuration for a CTM climatological run."""

    latitude_deg: float
    julian_day: int
    integration_days: int
    boxes: list[CtmBoxSpec]
    bromine: bool
    iodine: bool
    solar_flux_scale: float
    def __init__(
        self,
        *,
        latitude_deg: float = 60.0,
        julian_day: int = 75,
        integration_days: int = 40,
        boxes: Optional[list[CtmBoxSpec]] = None,
        bromine: bool = True,
        iodine: bool = False,
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
    @property
    def altitude_km(self) -> np.ndarray:
        """Box altitudes in km, with shape ``(n_boxes,)``."""
        ...
    @property
    def pressure_mb(self) -> np.ndarray:
        """Box pressures in mb, with shape ``(n_boxes,)``."""
        ...
    @property
    def temperature_k(self) -> np.ndarray:
        """Box temperatures in K, with shape ``(n_boxes,)``."""
        ...
    @property
    def air_density_cm3(self) -> np.ndarray:
        """Box air number densities in cm⁻³, with shape ``(n_boxes,)``."""
        ...
    @property
    def elapsed_seconds(self) -> np.ndarray:
        """Shared monotonic DIURN coordinate, with shape ``(n_timesteps,)``."""
        ...
    @property
    def time_hhmm(self) -> np.ndarray:
        """Shared cyclic local-time labels, with shape ``(n_timesteps,)``."""
        ...
    def species_profile(self, species_name: str) -> np.ndarray:
        """Return a daily-mean implicit-species profile with shape ``(n_boxes,)``."""
        ...
    def long_lived_profile(self, species_name: str) -> np.ndarray:
        """Return a daily-mean mixing-ratio profile with shape ``(n_boxes,)``."""
        ...
    def jvalue_profile(self, jvalue_name: str) -> np.ndarray:
        """Return a daily-mean J-value profile with shape ``(n_boxes,)``."""
        ...
    def species_grid(self, species_name: str) -> np.ndarray:
        """Return an implicit species as a 2-D array of shape ``(n_boxes, n_timesteps)``.

        Parameters
        ----------
        species_name:
            A name from :data:`pratmo.IMPLICIT_SPECIES_NAMES` (case-insensitive).
        """
        ...
    def jvalue_grid(self, jvalue_name: str) -> np.ndarray:
        """Return a J-value as a 2-D array of shape ``(n_boxes, n_timesteps)``.

        Parameters
        ----------
        jvalue_name:
            A name from :data:`pratmo.JVALUE_NAMES` (case-insensitive). Values
            are daily means repeated across the time dimension.
        """
        ...
    def __len__(self) -> int: ...

class CtmOutput:
    """Output from a CTM climatological run."""

    @property
    def boxes(self) -> list[BoxSnapshot]:
        """Final converged snapshot for each box."""
        ...
    @property
    def diagnostics(self) -> Diagnostics: ...
    @property
    def altitude_km(self) -> np.ndarray:
        """Box altitudes in km, with shape ``(n_boxes,)``."""
        ...
    @property
    def pressure_mb(self) -> np.ndarray:
        """Box pressures in mb, with shape ``(n_boxes,)``."""
        ...
    @property
    def temperature_k(self) -> np.ndarray:
        """Box temperatures in K, with shape ``(n_boxes,)``."""
        ...
    @property
    def air_density_cm3(self) -> np.ndarray:
        """Box air number densities in cm⁻³, with shape ``(n_boxes,)``."""
        ...
    def species_profile(self, species_name: str) -> np.ndarray:
        """Return an implicit species as a 1-D array of shape ``(n_boxes,)``.

        Names are case-insensitive and listed in
        :data:`pratmo.IMPLICIT_SPECIES_NAMES`.
        """
        ...
    def long_lived_profile(self, species_name: str) -> np.ndarray:
        """Return a mixing ratio as a 1-D array of shape ``(n_boxes,)``.

        Names are case-insensitive and listed in :data:`pratmo.LONG_LIVED_NAMES`.
        """
        ...
    def jvalue_profile(self, jvalue_name: str) -> np.ndarray:
        """Return a J-value as a 1-D array of shape ``(n_boxes,)``.

        Names are case-insensitive and listed in :data:`pratmo.JVALUE_NAMES`.
        """
        ...
    def __len__(self) -> int: ...

class No2ConstrainedDiurnConfig:
    """Iterative NOy scaling configuration for matching observed NO2."""

    diurn: DiurnConfig
    observed_no2_cm3: list[float]
    target_hhmm: int
    iterations: int
    def __init__(
        self,
        diurn: DiurnConfig,
        observed_no2_cm3: list[float],
        target_hhmm: int,
        iterations: int = 3,
    ) -> None: ...

class No2ConstrainedDiurnOutput:
    @property
    def output(self) -> DiurnOutput: ...
    @property
    def noy_scale(self) -> list[float]: ...
    @property
    def modeled_no2_cm3(self) -> list[float]: ...

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
    def run_diurn_no2_constrained(
        self,
        cfg: No2ConstrainedDiurnConfig,
    ) -> No2ConstrainedDiurnOutput:
        """Iteratively scale per-box NOy to match observed NO2 at a target HHMM."""
        ...
    def run_ctm(self, cfg: CtmConfig) -> CtmOutput:
        """Run the CTM climatological mode."""
        ...
    def __init__(self, data_dir: Optional[str] = None) -> None:
        """Create a model with embedded data or data loaded from *data_dir*."""
        ...
