---
jupytext:
  text_representation:
    format_name: myst
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# Atmospheric inputs and units

This guide shows the supported ways to specify atmospheric quantities. The
high-level interface converts everything to the exact units required by the
Rust core and validates the vertical grid before a run begins.

## Canonical model units

| Quantity | Internal unit | Accepted convenience units |
|---|---:|---|
| Pressure | hPa (mb) | Pa, hPa/mb, kPa, bar, atm |
| Temperature | K | K, °C |
| Altitude | km | m, km |
| Mixing ratio | fraction | fraction, %, ppmv, ppbv, pptv |
| Number density | molecules cm⁻³ | cm⁻³, m⁻³ |
| Aerosol surface area | µm² cm⁻³ | µm² cm⁻³, m² m⁻³ |

Helpers accept scalars, sequences, and NumPy arrays.

```{code-cell} ipython3
import numpy as np
from pratmo import (
    altitude,
    number_density,
    ppbv,
    ppmv,
    pptv,
    pressure,
    surface_area_density,
    temperature,
)

print("5 ppmv =", ppmv(5), "fraction")
print("300 ppbv =", ppbv(300), "fraction")
print("1.2 pptv =", pptv(1.2), "fraction")
print("5000 Pa =", pressure(5000, "Pa"), "hPa")
print("-55 degC =", temperature(-55, "degC"), "K")
print("25000 m =", altitude(25000, "m"), "km")
print("1e18 m-3 =", number_density(1e18, "m-3"), "cm-3")
```

## A complete custom atmosphere

`Atmosphere` requires pressure, temperature, and ozone. Altitude and aerosol
surface area are optional. Values must be ordered bottom to top: pressure
decreases and altitude increases.

Ozone defaults to **ppmv**, which is convenient for atmospheric profiles and
avoids the common mistake of passing values such as `5.0` as a dimensionless
fraction.

```{code-cell} ipython3
from pratmo import Atmosphere

atmosphere = Atmosphere(
    pressure=[100, 70, 45, 28, 17, 10],
    pressure_unit="hPa",
    temperature=[225, 222, 219, 216, 214, 216],
    temperature_unit="K",
    altitude=[16, 18, 21, 24, 27, 30],
    altitude_unit="km",
    ozone=[2.5, 3.8, 5.1, 6.4, 7.2, 7.0],
    ozone_unit="ppmv",
    aerosol_surface_area=[0.60, 0.45, 0.30, 0.18, 0.10, 0.05],
    aerosol_unit="um2/cm3",
)

print(atmosphere)
```

```{code-cell} ipython3
from pratmo.plotting import plot_atmosphere

plot_atmosphere(atmosphere)
```

Run all six profile levels as chemistry boxes by omitting `boxes`:

```python
from pratmo import Model

result = Model().diurnal(
    atmosphere=atmosphere,
    latitude=30.0,
    day="2026-06-21",
)
```

PRATMO supports up to 81 radiative levels but only 25 chemistry boxes in one
DIURN run. Radiative levels describe the column used by photolysis; chemistry
boxes are the selected locations where the mechanism is integrated. For a
denser radiative profile, explicitly select chemistry boxes.

## Pressure in pascals and temperature in Celsius

Units apply to an entire array. The following profile is identical in
structure but uses common retrieval-file units:

```{code-cell} ipython3
retrieval_atmosphere = Atmosphere(
    pressure=[10000, 7000, 4500, 2800, 1700, 1000],
    pressure_unit="Pa",
    temperature=[-48.15, -51.15, -54.15, -57.15, -59.15, -57.15],
    temperature_unit="degC",
    altitude=[16000, 18000, 21000, 24000, 27000, 30000],
    altitude_unit="m",
    ozone=[2.5, 3.8, 5.1, 6.4, 7.2, 7.0],
    ozone_unit="ppmv",
)

np.testing.assert_allclose(retrieval_atmosphere.pressure_mb, atmosphere.pressure_mb)
np.testing.assert_allclose(retrieval_atmosphere.temperature_k, atmosphere.temperature_k)
```

## Ozone as a dimensionless mixing ratio

Use `ozone_unit="fraction"` for data already expressed as mole or volume
fraction:

```python
fraction_atmosphere = Atmosphere(
    pressure=[80, 50, 30],
    temperature=[225, 220, 215],
    ozone=[3.5e-6, 5.0e-6, 6.0e-6],
    ozone_unit="fraction",
)
```

## Ozone as number density

Instrument products sometimes provide molecules cm⁻³ or molecules m⁻³. The
interface converts the latter and tells the core to skip mixing-ratio
conversion.

```python
ozone_density_atmosphere = Atmosphere(
    pressure=[80, 50, 30],
    temperature=[225, 220, 215],
    altitude=[18, 21, 24],
    ozone=[7.0e12, 5.5e12, 4.0e12],
    ozone_unit="cm-3",
)

ozone_density_si = Atmosphere(
    pressure=[8000, 5000, 3000],
    pressure_unit="Pa",
    temperature=[225, 220, 215],
    altitude=[18000, 21000, 24000],
    altitude_unit="m",
    ozone=[7.0e18, 5.5e18, 4.0e18],
    ozone_unit="m-3",
)
```

## Let PRATMO estimate altitude

Omit `altitude` to derive a hydrostatic grid from pressure and temperature:

```python
hydrostatic = Atmosphere(
    pressure=[100, 70, 45, 28],
    temperature=[225, 222, 219, 216],
    ozone=[2.5, 3.8, 5.1, 6.4],
    ozone_unit="ppmv",
)
```

The derived coordinate is relative: PRATMO anchors the first supplied level at
0 km and integrates upward. It preserves hydrostatic spacing but cannot infer
the absolute geometric altitude of a retrieval profile.

Exact `Box.at_altitude(...)` interpolation requires an explicit altitude grid;
hydrostatic profiles can still run every level or select `Box.at_level(...)`.

## Select levels or exact chemistry altitudes

Levels are 1-based indices into the custom radiative profile:

```python
from pratmo import Box, Model

selected = Model().diurnal(
    atmosphere=atmosphere,
    boxes=[Box.at_level(2), Box.at_level(4), Box.at_level(6)],
)
```

An exact altitude between two radiative shells interpolates pressure,
temperature, ozone, aerosol, and wavelength-resolved actinic flux:

```python
interpolated = Model().diurnal(
    atmosphere=atmosphere,
    boxes=[Box.at_altitude(22.5, unit="km")],
)
```

## Aerosol and sea-salt surface area

There are two related inputs:

- `Atmosphere(aerosol_surface_area=...)` defines a vertical sulfate aerosol
  profile for radiative extinction and for boxes created from profile levels.
- `Box(aerosol_surface_area_um2_cm3=...)` and
  `Box(sea_salt_surface_area_um2_cm3=...)` set the chemistry surface area for a
  particular box.

```python
from pratmo import Box, ChemistryOptions, PhotolysisOptions

volcanic = Model().diurnal(
    atmosphere=atmosphere,
    boxes=[
        Box.at_altitude(
            22.5,
            aerosol_surface_area_um2_cm3=8.0,
            sea_salt_surface_area_um2_cm3=0.02,
        )
    ],
    chemistry=ChemistryOptions(heterogeneous=True),
    photolysis=PhotolysisOptions(aerosol_extinction=True),
)
```

Heterogeneous chemistry and radiative aerosol are independent switches. This
allows chemistry-only, radiation-only, both, or neither sensitivity tests.
The complete control matrix is shown in {doc}`custom-atmosphere`.

## Long-lived initial mixing ratios

Supplying any custom long-lived state replaces the complete embedded state for
that box. Start from a complete profile rather than leaving most fields at
zero. `background_mixing_ratios()` supplies representative lower-stratospheric
values and accepts readable overrides:

```python
from pratmo import background_mixing_ratios, ppbv, ppmv, pptv

initial = background_mixing_ratios(
    o3=ppmv(5.0),
    n2o=ppbv(290),
    noy=ppbv(11),
    ch4=ppmv(1.65),
    h2o=ppmv(5.0),
    brx=pptv(20),
    iodx=pptv(1),
)

custom_initial = Model().diurnal(
    boxes=[Box.at_level(15)],
    initial_mixing_ratios=[initial],
)
```

`background_mixing_ratios()` is a complete numerical starting point, not an
altitude- or season-specific climatology. Do not use it as an observational
background without an explicit scientific justification.

For real profiles, {class}`pratmo.PratmoClimatology` samples the legacy
monthly temperature, O3, N2O, NOy, and aerosol tables and derives correlated
CH4, H2O, Bry, Cly, and CH3Cl inputs. See {doc}`custom-atmosphere` for the full
climatological and NO2-constrained workflows.

Atmospheric ozone and output ozone are easy to confuse: `Atmosphere.ozone`
sets the radiative/initial profile, `species_profile("o3")` returns modeled
number density, and `long_lived_profile("o3")` returns a mixing ratio. See
{doc}`species-and-units`.

## Convert output for display

Output remains in stable canonical units. Convert only at presentation or file
boundaries:

```python
from pratmo import mixing_ratio_as, number_density_as

ozone_per_m3 = number_density_as(result.species_profile("o3"), "m-3")
noy_ppbv = mixing_ratio_as(result.long_lived_profile("noy"), "ppbv")
water_ppmv = mixing_ratio_as(result.long_lived_profile("h2o"), "ppmv")
```
