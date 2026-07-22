"""Small, dependency-free helpers for PRATMO's canonical input units.

PRATMO intentionally does not require a general units package.  These helpers
make conversions visible at the call site and work with scalars, sequences,
and NumPy arrays.
"""

from __future__ import annotations

from collections.abc import Sequence
from typing import Union

import numpy as np


ArrayLike = Union[float, Sequence[float], np.ndarray]


def _return_like_input(value: ArrayLike, converted: np.ndarray) -> float | np.ndarray:
    if np.isscalar(value):
        return float(converted)
    return converted


def _normalized(unit: str) -> str:
    return (
        unit.strip()
        .lower()
        .replace(" ", "")
        .replace("²", "2")
        .replace("³", "3")
        .replace("⁻", "-")
        .replace("−", "-")
    )


def mixing_ratio(value: ArrayLike, unit: str = "fraction") -> float | np.ndarray:
    """Convert a volume/mole mixing ratio to PRATMO's dimensionless fraction.

    Accepted units are ``fraction`` (also ``vmr``), ``percent``, ``ppmv``,
    ``ppbv``, and ``pptv``.
    """

    scales = {
        "fraction": 1.0,
        "vmr": 1.0,
        "mol/mol": 1.0,
        "%": 1.0e-2,
        "percent": 1.0e-2,
        "ppmv": 1.0e-6,
        "ppm": 1.0e-6,
        "ppbv": 1.0e-9,
        "ppb": 1.0e-9,
        "pptv": 1.0e-12,
        "ppt": 1.0e-12,
    }
    key = _normalized(unit)
    if key not in scales:
        raise ValueError(
            f"Unknown mixing-ratio unit {unit!r}; use fraction, percent, ppmv, ppbv, or pptv"
        )
    array = np.asarray(value, dtype=float) * scales[key]
    return _return_like_input(value, array)


def ppmv(value: ArrayLike) -> float | np.ndarray:
    """Convert parts per million by volume to a dimensionless fraction."""

    return mixing_ratio(value, "ppmv")


def ppbv(value: ArrayLike) -> float | np.ndarray:
    """Convert parts per billion by volume to a dimensionless fraction."""

    return mixing_ratio(value, "ppbv")


def pptv(value: ArrayLike) -> float | np.ndarray:
    """Convert parts per trillion by volume to a dimensionless fraction."""

    return mixing_ratio(value, "pptv")


def pressure(value: ArrayLike, unit: str = "hPa") -> float | np.ndarray:
    """Convert pressure to PRATMO's canonical hPa (millibar)."""

    scales = {
        "hpa": 1.0,
        "mb": 1.0,
        "mbar": 1.0,
        "millibar": 1.0,
        "pa": 1.0e-2,
        "kpa": 10.0,
        "bar": 1.0e3,
        "atm": 1013.25,
    }
    key = _normalized(unit)
    if key not in scales:
        raise ValueError(f"Unknown pressure unit {unit!r}; use Pa, hPa/mb, kPa, bar, or atm")
    array = np.asarray(value, dtype=float) * scales[key]
    return _return_like_input(value, array)


def temperature(value: ArrayLike, unit: str = "K") -> float | np.ndarray:
    """Convert temperature to kelvin."""

    key = _normalized(unit)
    array = np.asarray(value, dtype=float)
    if key in {"k", "kelvin"}:
        converted = array
    elif key in {"c", "°c", "degc", "celsius"}:
        converted = array + 273.15
    else:
        raise ValueError(f"Unknown temperature unit {unit!r}; use K or degC")
    return _return_like_input(value, converted)


def altitude(value: ArrayLike, unit: str = "km") -> float | np.ndarray:
    """Convert geometric altitude to kilometres."""

    scales = {"km": 1.0, "m": 1.0e-3}
    key = _normalized(unit)
    if key not in scales:
        raise ValueError(f"Unknown altitude unit {unit!r}; use m or km")
    array = np.asarray(value, dtype=float) * scales[key]
    return _return_like_input(value, array)


def number_density(value: ArrayLike, unit: str = "cm-3") -> float | np.ndarray:
    """Convert number density to molecules cm⁻³.

    ``m-3`` values are divided by one million. ``cm-3`` and common spelling
    variants are passed through.
    """

    scales = {
        "cm-3": 1.0,
        "cm^-3": 1.0,
        "1/cm3": 1.0,
        "molecules/cm3": 1.0,
        "moleculescm-3": 1.0,
        "m-3": 1.0e-6,
        "m^-3": 1.0e-6,
        "1/m3": 1.0e-6,
        "molecules/m3": 1.0e-6,
        "moleculesm-3": 1.0e-6,
    }
    key = _normalized(unit)
    if key not in scales:
        raise ValueError(f"Unknown number-density unit {unit!r}; use cm-3 or m-3")
    array = np.asarray(value, dtype=float) * scales[key]
    return _return_like_input(value, array)


def surface_area_density(
    value: ArrayLike, unit: str = "um2/cm3"
) -> float | np.ndarray:
    """Convert aerosol surface-area density to µm² cm⁻³.

    A value in m² m⁻³ is multiplied by ``1e6``.
    """

    scales = {
        "um2/cm3": 1.0,
        "micron2/cm3": 1.0,
        "um2cm-3": 1.0,
        "um2cm^-3": 1.0,
        "m2/m3": 1.0e6,
        "m^-1": 1.0e6,
        "1/m": 1.0e6,
    }
    key = _normalized(unit).replace("µ", "u").replace("μ", "u")
    if key not in scales:
        raise ValueError(
            f"Unknown surface-area-density unit {unit!r}; use um2/cm3 or m2/m3"
        )
    array = np.asarray(value, dtype=float) * scales[key]
    return _return_like_input(value, array)


def mixing_ratio_as(value: ArrayLike, unit: str) -> float | np.ndarray:
    """Convert a dimensionless mixing ratio to a requested display unit."""

    one = float(mixing_ratio(1.0, unit))
    array = np.asarray(value, dtype=float) / one
    return _return_like_input(value, array)


def number_density_as(value: ArrayLike, unit: str) -> float | np.ndarray:
    """Convert molecules cm⁻³ to a requested display unit."""

    one = float(number_density(1.0, unit))
    array = np.asarray(value, dtype=float) / one
    return _return_like_input(value, array)


__all__ = [
    "altitude",
    "mixing_ratio",
    "mixing_ratio_as",
    "number_density",
    "number_density_as",
    "ppbv",
    "ppmv",
    "pptv",
    "pressure",
    "surface_area_density",
    "temperature",
]
