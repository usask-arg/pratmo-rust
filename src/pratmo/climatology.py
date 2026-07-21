"""PRATMO monthly climatologies for custom-atmosphere DIURN runs.

The tables are read from the original fixed-format ``fort03/04/05/51`` files.
They are numerically identical to the later ``PratmoClimatologies`` C++ tables,
but remain tracked as the authoritative model inputs.
"""

from __future__ import annotations

from dataclasses import dataclass
from datetime import date, timedelta
from pathlib import Path
import re

import numpy as np


_LATITUDES = np.arange(-85.0, 86.0, 10.0)
_NUMBER = re.compile(r"[-+]?(?:\d+\.?(?:\d*)?|\.\d+)(?:[Ee][-+]?\d+)?")


def _tagged_values(path: Path, tag: str, rows: str) -> list[float]:
    pattern = re.compile(rf"^{re.escape(tag)}[{rows}]\s", re.IGNORECASE)
    values: list[float] = []
    for line in path.read_text(errors="ignore").splitlines():
        stripped = line.lstrip()
        if pattern.match(stripped):
            values.extend(float(value) for value in _NUMBER.findall(stripped[2:]))
    return values


def _monthly_table(path: Path, tag: str, rows: str, heights: int) -> np.ndarray:
    expected = 12 * 18 * heights
    values = _tagged_values(path, tag, rows)[:expected]
    if len(values) != expected:
        raise ValueError(f"{path} contains {len(values)} {tag}-values; expected {expected}")
    return np.asarray(values).reshape(12, 18, heights)


def _aerosol_table(path: Path) -> np.ndarray:
    values: list[float] = []
    for line in path.read_text(errors="ignore").splitlines():
        fields = line.split()
        if fields and fields[0] in {"8", "9"}:
            values.extend(float(value) for value in fields[1:])
    expected = 4 * 18 * 17
    if len(values) != expected:
        raise ValueError(f"{path} contains {len(values)} aerosol values; expected {expected}")
    return np.asarray(values).reshape(4, 18, 17)


def _month_midpoints(day: date) -> tuple[date, date, float]:
    midpoint = day.replace(day=16)
    if day < midpoint:
        previous_day = midpoint.replace(day=1) - timedelta(days=1)
        lower = previous_day.replace(day=16)
        upper = midpoint
    else:
        next_month = (midpoint.replace(day=28) + timedelta(days=4)).replace(day=1)
        lower = midpoint
        upper = next_month.replace(day=16)
    weight = (day - lower).total_seconds() / (upper - lower).total_seconds()
    return lower, upper, weight


def _sample_month(table: np.ndarray, month: int, latitude: float, altitude_km: float) -> float:
    max_index = table.shape[2] - 1
    maximum_altitude = 2.0 * max_index
    latitude = float(np.clip(latitude, -85.0, 85.0))
    altitude = max(0.0, float(altitude_km))

    if latitude <= -85.0:
        lat0 = lat1 = 0
        lat_weight = 0.0
    elif latitude >= 85.0:
        lat0 = lat1 = 17
        lat_weight = 0.0
    else:
        lat0 = int(np.floor((latitude + 85.0) / 10.0))
        lat1 = lat0 + 1
        lat_weight = (latitude - _LATITUDES[lat0]) / 10.0

    scale = 1.0
    if altitude > maximum_altitude:
        ratio = table[month, lat0, max_index] / table[month, lat0, max_index - 1]
        scale = min(ratio, 1.0) ** (np.floor(altitude / 2.0) - max_index)
        altitude = maximum_altitude

    height0 = min(int(np.floor(altitude / 2.0)), max_index - 1)
    height1 = height0 + 1
    height_weight = (altitude - 2.0 * height0) / 2.0
    lower = (1.0 - lat_weight) * table[month, lat0, height0] + lat_weight * table[month, lat1, height0]
    upper = (1.0 - lat_weight) * table[month, lat0, height1] + lat_weight * table[month, lat1, height1]
    return float(((1.0 - height_weight) * lower + height_weight * upper) * scale)


def _pressure_height_grid(temperature_k: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return the C++ PRATMO pressure-level height and relative density grids."""
    temperature = np.asarray(temperature_k, dtype=float)
    pressure_ratio = 10.0 ** (-2.0 / 16.0)
    pressure_mb = 1000.0 * pressure_ratio ** np.arange(temperature.size)
    # This is the LOGP=1 branch of the C++ HYSTAT routine. Heights are in km.
    coefficient_km_per_k = -287.058 / (2.0 * 9.80665) * np.log(pressure_ratio) / 1000.0
    altitude_km = np.zeros(temperature.size)
    altitude_km[1:] = np.cumsum(
        (temperature[:-1] + temperature[1:]) * coefficient_km_per_k
    )
    relative_density = pressure_mb / temperature
    return altitude_km, relative_density


def _ch4_from_n2o(n2o_ppb: float) -> float:
    x = n2o_ppb
    x2 = x * x
    if 320.0 < x < 330.0:
        # Preserve the C++ compatibility kludge (including its pre-clamp x²).
        x = 320.0
    if 0.0 <= x < 45.0:
        value = 0.23370 + 0.028568 * x - 6.4358e-4 * x2 + 6.0186e-6 * x2 * x
    elif x < 100.0:
        value = 0.50499 + 6.9690e-3 * x - 3.0114e-5 * x2 + 7.2321e-8 * x2 * x
    elif x <= 205.0:
        value = 0.62567 + 4.0635e-3 * x - 5.6738e-6 * x2
    elif x < 280.0:
        value = 0.32720 + 4.3564e-3 * x
    elif x <= 320.0:
        value = 0.34841 + 4.2638e-3 * x
    else:
        value = 0.0
    growth = (1753.67 - 1730.975) / 4.0 / 1000.0
    return value * (1.666 + growth * 10.0) / 1.666


def _cfc11_from_n2o(n2o_ppb: float) -> float:
    if 130.0 <= n2o_ppb <= 310.0:
        return 252.92 - 4.5591 * n2o_ppb + 0.025644 * n2o_ppb**2 - 3.4639e-5 * n2o_ppb**3
    if n2o_ppb > 310.0:
        return 272.09
    return 17.52 * n2o_ppb / 130.0


def _bry_from_n2o(n2o_ppb: float, latitude: float) -> float:
    x = n2o_ppb * 310.37 / 319.02
    if abs(latitude) <= 30.0 and x >= 225.0:
        cfc11 = 39.92 - 0.33539 * x + 0.0033736 * x * x
    else:
        cfc11 = _cfc11_from_n2o(x)
    bry = 16.02 + 0.0033 * cfc11 - 5.305e-4 * cfc11**2 + 2.55e-6 * cfc11**3 - 5.37e-9 * cfc11**4
    bry = max(bry * 1.0296 ** (2003.5 - 1994.875), 0.0)
    if x < 130.0:
        n2o_grid = np.array([0.0, 20.0, 40.0, 60.0, 80.0, 100.0, 120.0, 140.0])
        age_grid = np.array([6.65, 6.10, 5.70, 5.50, 5.30, 5.10, 4.90, 4.65])
        bry /= 1.0296 ** (np.interp(x, n2o_grid, age_grid) - np.interp(130.0, n2o_grid, age_grid))
    return float(bry if abs(latitude) <= 30.0 else max(bry, 0.6))


@dataclass(frozen=True)
class PratmoClimatologyProfile:
    """Climatological values sampled onto an arbitrary altitude grid."""

    latitude_deg: float
    day: date
    altitude_km: np.ndarray
    temperature_k: np.ndarray
    o3: np.ndarray
    n2o: np.ndarray
    noy: np.ndarray
    aerosol_surface_area_um2_cm3: np.ndarray

    def initial_mixing_ratios(self, *, o3: np.ndarray | None = None) -> list[object]:
        """Build PRATMO long-lived inputs, optionally replacing climatological O3."""
        from pratmo._pratmo import LongLivedMixingRatios

        ozone = self.o3 if o3 is None else np.asarray(o3, dtype=float)
        if ozone.shape != self.n2o.shape:
            raise ValueError("o3 must have the same shape as the sampled altitude grid")
        result = []
        for ozone_value, n2o_value, noy_value in zip(ozone, self.n2o, self.noy):
            n2o_ppb = n2o_value * 1.0e9
            ch4 = _ch4_from_n2o(n2o_ppb) * 1.0e-6
            if 0.0 <= n2o_ppb <= 320.0:
                cly_ppb = max(
                    3.53876 - 2.67709e-3 * n2o_ppb - 1.91693e-5 * n2o_ppb**2 - 2.40584e-8 * n2o_ppb**3,
                    0.005,
                )
                ch3cl_ppb = 0.57898 - 0.00058507 * n2o_ppb - 7.219e-7 * n2o_ppb**2 - 9.2002e-9 * n2o_ppb**3
            else:
                cly_ppb = ch3cl_ppb = 0.0
            result.append(
                LongLivedMixingRatios(
                    o3=float(ozone_value), n2o=float(n2o_value), noy=float(noy_value),
                    ch4=ch4, co=3.0e-8, clx=cly_ppb * 1.0e-9,
                    cf2cl2=100.0e-12, cfcl3=100.0e-12, ccl4=100.0e-12,
                    ch3cl=ch3cl_ppb * 1.0e-9, ch3ccl3=100.0e-12, h2=0.5e-6,
                    h2o=max(7.0e-6 - 2.0 * ch4, 0.0), nh3=0.0, c5h8=0.0,
                    brx=max(_bry_from_n2o(n2o_ppb, self.latitude_deg) * 1.0e-12 + 2.0e-12, 0.5e-12),
                    ch3br=10.0e-12, ocs=0.0,
                    # The legacy climatology predates PRATMO's iodine extension.
                    # Retain its established background instead of interpreting
                    # an absent Iy field as a chemically singular zero family.
                    iodx=1.0e-12,
                )
            )
        return result


class PratmoClimatology:
    """Loader and interpolator for the PRATMO T/O3/N2O/NOy/ASA tables."""

    def __init__(self, data_dir: str | Path):
        directory = Path(data_dir)
        self._temperature = _monthly_table(directory / "fort03_LLM.x", "T", "1234", 41)
        self._o3 = _monthly_table(directory / "fort03_LLM.x", "Z", "123", 31) * 1.0e-6
        self._noy = _monthly_table(directory / "fort05.x", "n", "123", 31) * 1.0e-9
        self._n2o = _monthly_table(directory / "fort51.x", "n", "123", 31) * 1.0e-9
        self._aerosol = _aerosol_table(directory / "fort04.x")

    def sample(self, latitude_deg: float, day: date, altitude_km: np.ndarray) -> PratmoClimatologyProfile:
        """Sample all climatologies using the original C++ interpolation rules.

        The tabulated 0, 2, ..., 80 km coordinate is PRATMO's pseudo-height
        (pressure-level) coordinate. The C++ implementation first converts it
        to geometric height with ``HYSTAT`` and only then samples a requested
        atmosphere. Sampling the tables directly at geometric altitude gives a
        different NOy/N2O profile in the stratosphere.
        """
        altitudes = np.asarray(altitude_km, dtype=float)
        lower, upper, month_weight = _month_midpoints(day)

        pseudo_altitude = np.arange(41, dtype=float) * 2.0

        def monthly_on_pseudo_grid(table: np.ndarray) -> np.ndarray:
            return np.asarray([
                (1.0 - month_weight) * _sample_month(table, lower.month - 1, latitude_deg, altitude)
                + month_weight * _sample_month(table, upper.month - 1, latitude_deg, altitude)
                for altitude in pseudo_altitude
            ])

        temperature_grid = monthly_on_pseudo_grid(self._temperature)
        pressure_altitude, relative_density = _pressure_height_grid(temperature_grid)
        o3_grid = monthly_on_pseudo_grid(self._o3)
        n2o_grid = monthly_on_pseudo_grid(self._n2o)
        noy_grid = monthly_on_pseudo_grid(self._noy)

        # Above the 60 km top of the trace-gas tables, C++ continues number
        # density with the final tabulated level-to-level scale factor.
        for values in (o3_grid, n2o_grid, noy_grid):
            scale = min(
                values[30] * relative_density[30]
                / (values[29] * relative_density[29]),
                1.0,
            )
            for index in range(31, len(values)):
                values[index] = (
                    values[index - 1]
                    * relative_density[index - 1]
                    * scale
                    / relative_density[index]
                )

        season = 0 if day.month in {12, 1, 2} else 1 if day.month <= 5 else 2 if day.month <= 8 else 3
        latitude_index = int(np.clip(np.floor((latitude_deg + 90.0) / 10.0 - 1.0e-8), 0, 17))
        aerosol = []
        for altitude in pseudo_altitude:
            if altitude <= 8.0:
                value = self._aerosol[season, latitude_index, 0]
            elif altitude <= 40.0:
                value = self._aerosol[season, latitude_index, min(int(np.floor((altitude - 8.0) / 2.0)), 16)]
            else:
                value = self._aerosol[season, latitude_index, 16] * np.exp(-(altitude - 40.0) / 3.0)
            aerosol.append(value)
        aerosol_grid = np.asarray(aerosol)
        for index in range(21, len(aerosol_grid)):
            dz = pressure_altitude[index] - pressure_altitude[index - 1]
            aerosol_grid[index] = aerosol_grid[index - 1] * np.exp(-dz / 3.0)

        def at_geometric_height(values: np.ndarray) -> np.ndarray:
            return np.interp(altitudes, pressure_altitude, values)

        return PratmoClimatologyProfile(
            latitude_deg=float(latitude_deg), day=day, altitude_km=altitudes.copy(),
            temperature_k=at_geometric_height(temperature_grid),
            o3=at_geometric_height(o3_grid),
            n2o=at_geometric_height(n2o_grid),
            noy=at_geometric_height(noy_grid),
            aerosol_surface_area_um2_cm3=at_geometric_height(aerosol_grid),
        )
