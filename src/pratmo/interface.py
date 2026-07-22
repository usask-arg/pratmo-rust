"""User-oriented Python interface for configuring and running PRATMO."""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass
from datetime import date, datetime
from pathlib import Path
import warnings

import numpy as np

from pratmo._pratmo import (
    CtmBoxSpec,
    CtmConfig,
    CtmOutput,
    CustomAtmosphereProfile,
    DiurnBoxSpec,
    DiurnConfig,
    DiurnOutput,
    LongLivedMixingRatios,
    No2ConstrainedDiurnConfig,
    No2ConstrainedDiurnOutput,
    PratmoModel,
)
from pratmo.quantities import (
    altitude as convert_altitude,
    mixing_ratio,
    number_density,
    pressure as convert_pressure,
    surface_area_density,
    temperature as convert_temperature,
)


class PratmoWarning(UserWarning):
    """Base category for suspicious but technically valid PRATMO inputs."""


class ExperimentalFeatureWarning(PratmoWarning):
    """Warning emitted when an experimental mechanism is enabled."""


def _array(name: str, value: Sequence[float] | np.ndarray) -> np.ndarray:
    result = np.asarray(value, dtype=float)
    if result.ndim != 1:
        raise ValueError(f"{name} must be one-dimensional, got shape {result.shape}")
    if result.size == 0:
        raise ValueError(f"{name} must contain at least one value")
    if not np.all(np.isfinite(result)):
        raise ValueError(f"{name} must contain only finite values")
    return result


def _day_of_year(value: int | date | datetime | str) -> int:
    if isinstance(value, str):
        try:
            value = date.fromisoformat(value)
        except ValueError as exc:
            raise ValueError("day must be a day-of-year integer or an ISO date such as '2026-03-20'") from exc
    if isinstance(value, datetime):
        value = value.date()
    if isinstance(value, date):
        return value.timetuple().tm_yday
    if isinstance(value, bool):
        raise ValueError("day-of-year must be an integer from 1 through 366")
    result = int(value)
    if result != value or not 1 <= result <= 366:
        raise ValueError("day-of-year must be an integer from 1 through 366")
    return result


@dataclass(frozen=True)
class Box:
    """A chemistry box on a standard level or at an exact custom altitude.

    Parameters are in PRATMO's canonical units. Use :mod:`pratmo.quantities`
    helpers when source data uses other units.
    """

    level: int | None = None
    altitude_km: float | None = None
    aerosol_surface_area_um2_cm3: float | None = None
    sea_salt_surface_area_um2_cm3: float = 0.0
    temperature_offset_k: float = 0.0

    def __post_init__(self) -> None:
        if self.level is None and self.altitude_km is None:
            raise ValueError("Box requires either level= or altitude_km=")
        if self.level is not None:
            if not isinstance(self.level, int) or isinstance(self.level, bool):
                raise TypeError("Box level must be an integer")
            if not 1 <= self.level <= 81:
                raise ValueError("Box level must be between 1 and 81")
        if self.altitude_km is not None and (
            not np.isfinite(self.altitude_km) or self.altitude_km < 0.0
        ):
            raise ValueError("Box altitude_km must be finite and non-negative")
        for name, value in (
            ("aerosol_surface_area_um2_cm3", self.aerosol_surface_area_um2_cm3),
            ("sea_salt_surface_area_um2_cm3", self.sea_salt_surface_area_um2_cm3),
        ):
            if value is not None and (not np.isfinite(value) or value < 0.0):
                raise ValueError(f"Box {name} must be finite and non-negative")
        if not np.isfinite(self.temperature_offset_k):
            raise ValueError("Box temperature_offset_k must be finite")
        if abs(self.temperature_offset_k) > 30.0:
            warnings.warn(
                f"A temperature offset of {self.temperature_offset_k:g} K is unusually large; "
                "check that an absolute temperature was not supplied by mistake.",
                PratmoWarning,
                stacklevel=2,
            )

    @classmethod
    def at_level(cls, level: int, **kwargs: float) -> "Box":
        """Create a box at a 1-based PRATMO pressure level."""

        return cls(level=level, **kwargs)

    @classmethod
    def at_altitude(cls, value: float, *, unit: str = "km", **kwargs: float) -> "Box":
        """Create a box at an exact altitude in a custom atmosphere."""

        return cls(altitude_km=float(convert_altitude(value, unit)), **kwargs)


class Atmosphere:
    """Validated custom pressure, temperature, ozone, altitude, and aerosol profile.

    Unit names are explicit and converted once on construction. Ozone can be
    supplied as a mixing ratio (``fraction``, ``ppmv``, ``ppbv``, ``pptv``) or
    a number density (``cm-3`` or ``m-3``).
    """

    def __init__(
        self,
        *,
        pressure: Sequence[float] | np.ndarray,
        temperature: Sequence[float] | np.ndarray,
        ozone: Sequence[float] | np.ndarray,
        pressure_unit: str = "hPa",
        temperature_unit: str = "K",
        ozone_unit: str = "ppmv",
        altitude: Sequence[float] | np.ndarray | None = None,
        altitude_unit: str = "km",
        aerosol_surface_area: Sequence[float] | np.ndarray | None = None,
        aerosol_unit: str = "um2/cm3",
    ) -> None:
        self.pressure_mb = _array(
            "pressure", convert_pressure(pressure, pressure_unit)
        )
        self.temperature_k = _array(
            "temperature", convert_temperature(temperature, temperature_unit)
        )
        ozone_key = ozone_unit.strip().lower().replace(" ", "")
        if ozone_key in {
            "cm-3",
            "cm^-3",
            "cm⁻³",
            "1/cm3",
            "m-3",
            "m^-3",
            "m⁻³",
            "1/m3",
            "molecules/cm3",
            "molecules/m3",
        }:
            self.ozone = _array("ozone", number_density(ozone, ozone_unit))
            self.ozone_kind = "number_density"
        else:
            self.ozone = _array("ozone", mixing_ratio(ozone, ozone_unit))
            self.ozone_kind = "mixing_ratio"
        self.altitude_km = (
            None
            if altitude is None
            else _array("altitude", convert_altitude(altitude, altitude_unit))
        )
        self.aerosol_surface_area_um2_cm3 = (
            None
            if aerosol_surface_area is None
            else _array(
                "aerosol_surface_area",
                surface_area_density(aerosol_surface_area, aerosol_unit),
            )
        )
        self._validate()

    def _validate(self) -> None:
        n = self.pressure_mb.size
        arrays = {
            "temperature": self.temperature_k,
            "ozone": self.ozone,
            "altitude": self.altitude_km,
            "aerosol_surface_area": self.aerosol_surface_area_um2_cm3,
        }
        for name, values in arrays.items():
            if values is not None and values.size != n:
                raise ValueError(
                    f"{name} has {values.size} values but pressure has {n}; every profile must use the same grid"
                )
        if n > 81:
            raise ValueError("PRATMO supports at most 81 custom atmospheric levels")
        if np.any(self.pressure_mb <= 0.0):
            raise ValueError("pressure values must be positive")
        if n > 1 and np.any(np.diff(self.pressure_mb) >= 0.0):
            raise ValueError("pressure must decrease from the bottom to the top of the profile")
        if np.any(self.temperature_k <= 0.0):
            raise ValueError("temperature values must be above absolute zero")
        if np.any(self.ozone < 0.0):
            raise ValueError("ozone values must be non-negative")
        if self.altitude_km is not None:
            if np.any(self.altitude_km < 0.0):
                raise ValueError("altitude values must be non-negative")
            if n > 1 and np.any(np.diff(self.altitude_km) <= 0.0):
                raise ValueError("altitude must increase from the bottom to the top of the profile")
        if self.aerosol_surface_area_um2_cm3 is not None and np.any(
            self.aerosol_surface_area_um2_cm3 < 0.0
        ):
            raise ValueError("aerosol surface-area values must be non-negative")

        if np.any((self.temperature_k < 140.0) | (self.temperature_k > 330.0)):
            warnings.warn(
                "Some atmospheric temperatures lie outside 140–330 K. Check the temperature unit and profile.",
                PratmoWarning,
                stacklevel=3,
            )
        if self.ozone_kind == "mixing_ratio" and np.max(self.ozone) > 1.0e-2:
            warnings.warn(
                "Ozone exceeds a 1% mixing ratio. If these values are in ppmv, set ozone_unit='ppmv'.",
                PratmoWarning,
                stacklevel=3,
            )
        if self.ozone_kind == "number_density" and np.max(self.ozone) < 1.0e5:
            warnings.warn(
                "Ozone number densities are very small. Check whether mixing ratios were supplied as cm-3.",
                PratmoWarning,
                stacklevel=3,
            )
        if (
            self.aerosol_surface_area_um2_cm3 is not None
            and np.max(self.aerosol_surface_area_um2_cm3) > 100.0
        ):
            warnings.warn(
                "Aerosol surface area exceeds 100 µm² cm⁻³; check the aerosol unit.",
                PratmoWarning,
                stacklevel=3,
            )

    def to_native(self) -> CustomAtmosphereProfile:
        """Return the lower-level extension profile in canonical units."""

        return CustomAtmosphereProfile(
            pressure_mb=self.pressure_mb.tolist(),
            temperature_k=self.temperature_k.tolist(),
            altitude_km=None if self.altitude_km is None else self.altitude_km.tolist(),
            o3=self.ozone.tolist(),
            o3_kind=self.ozone_kind,
            aerosol_surface_area_um2_cm3=(
                None
                if self.aerosol_surface_area_um2_cm3 is None
                else self.aerosol_surface_area_um2_cm3.tolist()
            ),
        )

    def __len__(self) -> int:
        return int(self.pressure_mb.size)

    def __repr__(self) -> str:
        altitude = "provided altitude" if self.altitude_km is not None else "hydrostatic altitude"
        return f"Atmosphere({len(self)} levels, ozone={self.ozone_kind}, {altitude})"


@dataclass(frozen=True)
class ChemistryOptions:
    """Chemical mechanisms used by a high-level model run."""

    bromine: bool = True
    iodine: bool = False
    heterogeneous: bool = True

    def __post_init__(self) -> None:
        for name, value in (
            ("bromine", self.bromine),
            ("iodine", self.iodine),
            ("heterogeneous", self.heterogeneous),
        ):
            if not isinstance(value, bool):
                raise TypeError(f"{name} must be bool")


@dataclass(frozen=True)
class PhotolysisOptions:
    """Photolysis and radiative-transfer options."""

    solar_flux_scale: float = 1.0
    surface_albedo: float = 0.20
    aerosol_extinction: bool | None = None

    def __post_init__(self) -> None:
        if not np.isfinite(self.solar_flux_scale) or self.solar_flux_scale <= 0.0:
            raise ValueError("solar_flux_scale must be finite and positive")
        if not np.isfinite(self.surface_albedo) or not 0.0 <= self.surface_albedo <= 1.0:
            raise ValueError("surface_albedo must be between 0 and 1")
        if not 0.5 <= self.solar_flux_scale <= 1.5:
            warnings.warn(
                f"solar_flux_scale={self.solar_flux_scale:g} is far from the nominal value 1.0.",
                PratmoWarning,
                stacklevel=2,
            )
        if self.surface_albedo > 0.9:
            warnings.warn(
                "surface_albedo above 0.9 represents an unusually reflective lower boundary.",
                PratmoWarning,
                stacklevel=2,
            )
        if self.aerosol_extinction is not None and not isinstance(
            self.aerosol_extinction, bool
        ):
            raise TypeError("aerosol_extinction must be bool or None")


@dataclass(frozen=True)
class DiurnalOptions:
    """Numerical controls for a diurnal-cycle run."""

    integration_days: int = 20
    parallel_boxes: bool | None = None
    cpp_compatibility: bool = False
    elapsed_time_hours: Sequence[float] | None = None

    def __post_init__(self) -> None:
        if self.parallel_boxes is not None and not isinstance(self.parallel_boxes, bool):
            raise TypeError("parallel_boxes must be bool or None")
        if not isinstance(self.cpp_compatibility, bool):
            raise TypeError("cpp_compatibility must be bool")


@dataclass(frozen=True)
class CtmOptions:
    """Numerical controls for a climatological steady-state run."""

    integration_days: int = 40


def background_mixing_ratios(**overrides: float) -> LongLivedMixingRatios:
    """Return a complete, representative lower-stratospheric initial state.

    Values are dimensionless. Use :func:`pratmo.ppmv`, :func:`pratmo.ppbv`,
    and :func:`pratmo.pptv` for readable overrides. For quantitative work,
    prefer altitude- and season-specific climatological or measured inputs.
    """

    values = {
        "o3": 5.0e-6,
        "n2o": 300.0e-9,
        "noy": 10.0e-9,
        "ch4": 1.7e-6,
        "co": 50.0e-9,
        "clx": 3.0e-9,
        "cf2cl2": 500.0e-12,
        "cfcl3": 250.0e-12,
        "ccl4": 90.0e-12,
        "ch3cl": 550.0e-12,
        "ch3ccl3": 20.0e-12,
        "h2": 500.0e-9,
        "h2o": 5.0e-6,
        "nh3": 0.0,
        "c5h8": 0.0,
        "brx": 20.0e-12,
        "ch3br": 10.0e-12,
        "ocs": 500.0e-12,
        "iodx": 1.0e-12,
    }
    unknown = set(overrides) - set(values)
    if unknown:
        raise TypeError(f"Unknown long-lived field(s): {', '.join(sorted(unknown))}")
    values.update(overrides)
    return LongLivedMixingRatios(**values)


class Model:
    """Recommended Python entry point with runnable defaults and input checks."""

    def __init__(self, data_dir: str | Path | None = None) -> None:
        self._native_model = (
            PratmoModel.with_defaults()
            if data_dir is None
            else PratmoModel.from_data_dir(Path(data_dir))
        )

    @staticmethod
    def _check_days(days: int, mode: str) -> int:
        if not isinstance(days, int) or isinstance(days, bool) or days < 1:
            raise ValueError(f"{mode} integration_days must be a positive integer")
        if days < 3:
            warnings.warn(
                f"{mode} will integrate only {days} day(s); short runs often do not reach photochemical equilibrium. "
                "Inspect output.diagnostics before interpreting the result.",
                PratmoWarning,
                stacklevel=3,
            )
        return days

    @staticmethod
    def _check_latitude(latitude: float) -> float:
        value = float(latitude)
        if not np.isfinite(value) or not -90.0 <= value <= 90.0:
            raise ValueError("latitude must be finite and between -90 and 90 degrees")
        return value

    @staticmethod
    def _warn_chemistry(options: ChemistryOptions) -> None:
        if options.iodine:
            warnings.warn(
                "The iodine mechanism is experimental and has not been scientifically validated.",
                ExperimentalFeatureWarning,
                stacklevel=3,
            )

    @staticmethod
    def _resolved_boxes(
        boxes: Sequence[Box] | None,
        *,
        defaults: Sequence[int],
    ) -> list[Box]:
        return [Box.at_level(level) for level in defaults] if boxes is None else list(boxes)

    def ctm(
        self,
        *,
        latitude: float = 45.0,
        day: int | date | datetime | str = 80,
        boxes: Sequence[Box] | None = None,
        chemistry: ChemistryOptions = ChemistryOptions(),
        options: CtmOptions = CtmOptions(),
        solar_flux_scale: float = 1.0,
    ) -> CtmOutput:
        """Run the climatological steady-state workflow.

        By default, four standard levels spanning the lower and middle
        stratosphere are integrated for 40 days.
        """

        latitude_deg = self._check_latitude(latitude)
        ctm_latitude_index = int(np.floor((latitude_deg + 90.0) / 2.5))
        ctm_latitude_index = min(71, max(1, ctm_latitude_index))
        ctm_grid_latitude = -90.0 + 2.5 * ctm_latitude_index
        if not np.isclose(latitude_deg, ctm_grid_latitude):
            warnings.warn(
                f"CTM uses its {ctm_grid_latitude:g}° latitude grid point for "
                f"the requested latitude {latitude_deg:g}°. DIURN accepts the exact latitude.",
                PratmoWarning,
                stacklevel=2,
            )
        julian_day = _day_of_year(day)
        days = self._check_days(options.integration_days, "CTM")
        self._warn_chemistry(chemistry)
        if not np.isfinite(solar_flux_scale) or solar_flux_scale <= 0.0:
            raise ValueError("solar_flux_scale must be finite and positive")
        if not 0.5 <= solar_flux_scale <= 1.5:
            warnings.warn(
                f"solar_flux_scale={solar_flux_scale:g} is far from the nominal value 1.0.",
                PratmoWarning,
                stacklevel=2,
            )
        selected = self._resolved_boxes(boxes, defaults=(10, 15, 20, 25))
        if not selected:
            raise ValueError("CTM requires at least one Box")
        native_boxes = []
        for box in selected:
            if box.altitude_km is not None:
                raise ValueError("CTM boxes use standard level= indices; exact altitude requires a custom DIURN atmosphere")
            native_boxes.append(
                CtmBoxSpec(
                    altitude_level=box.level,
                    aerosol_surface_area_um2_cm3=(
                        0.0
                        if (
                            not chemistry.heterogeneous
                            or box.aerosol_surface_area_um2_cm3 is None
                        )
                        else box.aerosol_surface_area_um2_cm3
                    ),
                    sea_salt_surface_area_um2_cm3=(
                        box.sea_salt_surface_area_um2_cm3
                        if chemistry.heterogeneous
                        else 0.0
                    ),
                    temp_offset_k=box.temperature_offset_k,
                )
            )
        config = CtmConfig(
            latitude_deg=latitude_deg,
            julian_day=julian_day,
            integration_days=days,
            boxes=native_boxes,
            bromine=chemistry.bromine,
            iodine=chemistry.iodine,
            solar_flux_scale=solar_flux_scale,
        )
        return self._native_model.run_ctm(config)

    def _diurnal_config(
        self,
        *,
        latitude: float = 45.0,
        day: int | date | datetime | str = 80,
        boxes: Sequence[Box] | None = None,
        atmosphere: Atmosphere | None = None,
        initial_mixing_ratios: Sequence[LongLivedMixingRatios] | None = None,
        chemistry: ChemistryOptions = ChemistryOptions(),
        photolysis: PhotolysisOptions = PhotolysisOptions(),
        options: DiurnalOptions = DiurnalOptions(),
    ) -> DiurnConfig:
        """Build the validated native configuration for a DIURN workflow."""

        latitude_deg = self._check_latitude(latitude)
        julian_day = _day_of_year(day)
        days = self._check_days(options.integration_days, "DIURN")
        self._warn_chemistry(chemistry)
        if boxes is None:
            if atmosphere is None:
                selected = [Box.at_level(15)]
            else:
                selected = [
                    Box.at_level(
                        index + 1,
                        aerosol_surface_area_um2_cm3=(
                            None
                            if atmosphere.aerosol_surface_area_um2_cm3 is None
                            else float(atmosphere.aerosol_surface_area_um2_cm3[index])
                        ),
                    )
                    for index in range(len(atmosphere))
                ]
        else:
            selected = list(boxes)
        if atmosphere is None and not selected:
            raise ValueError("DIURN requires at least one Box when atmosphere is not provided")

        native_boxes = []
        for box in selected:
            level = box.level
            if box.altitude_km is not None:
                if atmosphere is None or atmosphere.altitude_km is None:
                    raise ValueError("Box.at_altitude requires an Atmosphere with an explicit altitude grid")
                if not atmosphere.altitude_km[0] <= box.altitude_km <= atmosphere.altitude_km[-1]:
                    raise ValueError("Box altitude lies outside the custom atmosphere")
                if level is None:
                    level = int(np.searchsorted(atmosphere.altitude_km, box.altitude_km, side="right"))
                    level = max(1, min(level, len(atmosphere)))
            if level is None:
                raise ValueError("Box requires a level")
            native_boxes.append(
                DiurnBoxSpec(
                    altitude_level=level,
                    altitude_km=box.altitude_km,
                    aerosol_surface_area_um2_cm3=(
                        0.25
                        if box.aerosol_surface_area_um2_cm3 is None
                        else box.aerosol_surface_area_um2_cm3
                    ),
                    sea_salt_surface_area_um2_cm3=box.sea_salt_surface_area_um2_cm3,
                    temp_offset_k=box.temperature_offset_k,
                )
            )

        if len(selected) > 25 or (not selected and atmosphere is not None and len(atmosphere) > 25):
            raise ValueError("DIURN supports at most 25 chemistry boxes; select a subset of custom levels")
        nbox = len(native_boxes)
        if initial_mixing_ratios is not None and len(initial_mixing_ratios) != nbox:
            raise ValueError(f"initial_mixing_ratios must contain one entry per box ({nbox})")
        if initial_mixing_ratios is not None:
            for index, ratios in enumerate(initial_mixing_ratios):
                values = ratios.to_dict()
                if sum(value == 0.0 for value in values.values()) > 10:
                    warnings.warn(
                        f"initial_mixing_ratios[{index}] contains many zeros. Use background_mixing_ratios() "
                        "or a climatological profile as a complete starting state.",
                        PratmoWarning,
                        stacklevel=2,
                    )
                if values["h2o"] > 1.0e-3:
                    warnings.warn(
                        f"initial_mixing_ratios[{index}].h2o exceeds 0.1%; check whether ppmv was supplied as a fraction.",
                        PratmoWarning,
                        stacklevel=2,
                    )

        has_profile_aerosol = (
            atmosphere is not None
            and atmosphere.aerosol_surface_area_um2_cm3 is not None
            and bool(np.any(atmosphere.aerosol_surface_area_um2_cm3 > 0.0))
        )
        radiative_aerosol = (
            has_profile_aerosol
            if photolysis.aerosol_extinction is None
            else photolysis.aerosol_extinction
        )
        if photolysis.aerosol_extinction is True and not has_profile_aerosol:
            warnings.warn(
                "aerosol_extinction=True has no profile-level aerosol surface area to use. "
                "Per-box aerosol affects heterogeneous chemistry but does not define a radiative profile.",
                PratmoWarning,
                stacklevel=2,
            )
        parallel = (
            nbox > 1 if options.parallel_boxes is None else options.parallel_boxes
        )
        config = DiurnConfig(
            latitude_deg=latitude_deg,
            julian_day=julian_day,
            integration_days=days,
            boxes=native_boxes,
            bromine=chemistry.bromine,
            iodine=chemistry.iodine,
            parallel_boxes=parallel,
            cpp_compatibility=options.cpp_compatibility,
            elapsed_time_hours=(
                None
                if options.elapsed_time_hours is None
                else list(options.elapsed_time_hours)
            ),
            solar_flux_scale=photolysis.solar_flux_scale,
            surface_albedo=photolysis.surface_albedo,
            heterogeneous_chemistry=chemistry.heterogeneous,
            radiative_aerosol=radiative_aerosol,
            atmosphere=None if atmosphere is None else atmosphere.to_native(),
            initial_mixing_ratios=(
                None if initial_mixing_ratios is None else list(initial_mixing_ratios)
            ),
        )
        return config

    def diurnal(
        self,
        *,
        latitude: float = 45.0,
        day: int | date | datetime | str = 80,
        boxes: Sequence[Box] | None = None,
        atmosphere: Atmosphere | None = None,
        initial_mixing_ratios: Sequence[LongLivedMixingRatios] | None = None,
        chemistry: ChemistryOptions = ChemistryOptions(),
        photolysis: PhotolysisOptions = PhotolysisOptions(),
        options: DiurnalOptions = DiurnalOptions(),
    ) -> DiurnOutput:
        """Run a resolved 24-hour photochemical cycle."""

        config = self._diurnal_config(
            latitude=latitude,
            day=day,
            boxes=boxes,
            atmosphere=atmosphere,
            initial_mixing_ratios=initial_mixing_ratios,
            chemistry=chemistry,
            photolysis=photolysis,
            options=options,
        )
        return self._native_model.run_diurn(config)

    def diurnal_no2_constrained(
        self,
        *,
        observed_no2_cm3: Sequence[float],
        target_hhmm: int,
        iterations: int = 3,
        latitude: float = 45.0,
        day: int | date | datetime | str = 80,
        boxes: Sequence[Box] | None = None,
        atmosphere: Atmosphere | None = None,
        initial_mixing_ratios: Sequence[LongLivedMixingRatios] | None = None,
        chemistry: ChemistryOptions = ChemistryOptions(),
        photolysis: PhotolysisOptions = PhotolysisOptions(),
        options: DiurnalOptions = DiurnalOptions(),
    ) -> No2ConstrainedDiurnOutput:
        """Scale total NOy so modeled NO2 matches observations at a local time.

        ``observed_no2_cm3`` must contain one non-negative number density per
        chemistry box. ``target_hhmm`` is local solar time in HHMM notation.
        """

        observed = _array("observed_no2_cm3", observed_no2_cm3)
        if np.any(observed < 0.0):
            raise ValueError("observed_no2_cm3 must be non-negative")
        if not isinstance(target_hhmm, int) or isinstance(target_hhmm, bool):
            raise TypeError("target_hhmm must be an integer in HHMM notation")
        hours, minutes = divmod(target_hhmm, 100)
        if not 0 <= hours <= 23 or not 0 <= minutes <= 59:
            raise ValueError("target_hhmm must be a valid local time from 0000 through 2359")
        if not isinstance(iterations, int) or isinstance(iterations, bool) or iterations < 1:
            raise ValueError("iterations must be a positive integer")

        config = self._diurnal_config(
            latitude=latitude,
            day=day,
            boxes=boxes,
            atmosphere=atmosphere,
            initial_mixing_ratios=initial_mixing_ratios,
            chemistry=chemistry,
            photolysis=photolysis,
            options=options,
        )
        if len(observed) != len(config.boxes):
            raise ValueError(
                f"observed_no2_cm3 must contain one value per chemistry box ({len(config.boxes)})"
            )
        constrained = No2ConstrainedDiurnConfig(
            diurn=config,
            observed_no2_cm3=observed.tolist(),
            target_hhmm=target_hhmm,
            iterations=iterations,
        )
        return self._native_model.run_diurn_no2_constrained(constrained)

    @property
    def native(self) -> PratmoModel:
        """Access the lower-level extension object for advanced workflows."""

        return self._native_model


__all__ = [
    "Atmosphere",
    "Box",
    "ChemistryOptions",
    "CtmOptions",
    "DiurnalOptions",
    "ExperimentalFeatureWarning",
    "Model",
    "PhotolysisOptions",
    "PratmoWarning",
    "background_mixing_ratios",
]
