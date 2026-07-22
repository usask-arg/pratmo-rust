---
jupytext:
  text_representation:
    format_name: myst
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# Custom atmospheres and observed NO2

Use {class}`pratmo.Atmosphere` when pressure, temperature, ozone, altitude, or
aerosol comes from a measurement or external model. DIURN can create one
chemistry box per supplied level or select a subset.

## Build and inspect a profile

```{code-cell} ipython3
from pratmo import Atmosphere, Model
from pratmo.plotting import plot_atmosphere

atmosphere = Atmosphere(
    pressure=[80.0, 60.0, 45.0, 33.0, 24.0],
    pressure_unit="hPa",
    temperature=[225.0, 222.0, 219.0, 216.0, 214.0],
    temperature_unit="K",
    ozone=[3.5, 4.4, 5.2, 5.8, 6.2],
    ozone_unit="ppmv",
    altitude=[18.0, 19.5, 21.0, 22.5, 24.0],
    altitude_unit="km",
    aerosol_surface_area=[0.35, 0.28, 0.20, 0.13, 0.08],
    aerosol_unit="um2/cm3",
)

plot_atmosphere(atmosphere)
```

Pressure must decrease with index, altitude must increase, and all arrays must
have the same length. Units are converted once when `Atmosphere` is created.

## Run every supplied level

Omitting `boxes` makes every custom radiative level a chemistry box, up to the
25-box limit.

```{code-cell} ipython3
:execution_timeout: 600

model = Model()
output = model.diurnal(
    latitude=0.0,
    day="2026-04-30",
    atmosphere=atmosphere,
)

print(output.altitude_km)
print(output.species_grid("no2").shape)
```

```{code-cell} ipython3
from pratmo.plotting import plot_profile

plot_profile(output, ["o3", "no2", "hno3"])
```

## Select levels or exact altitudes

On a custom atmosphere, level numbers are 1-based indices into the supplied
arrays—not the embedded standard grid.

```python
from pratmo import Box

selected = model.diurnal(
    atmosphere=atmosphere,
    boxes=[Box.at_level(1), Box.at_level(3), Box.at_level(5)],
)

interpolated = model.diurnal(
    atmosphere=atmosphere,
    boxes=[Box.at_altitude(20.3, unit="km")],
)
```

Exact-altitude boxes interpolate pressure, temperature, ozone, aerosol, and
the wavelength-dependent radiation field between surrounding radiative levels.

## Hydrostatic altitude

Altitude may be omitted. PRATMO then derives relative height hydrostatically
from pressure and temperature and anchors the bottom profile row at 0 km.

```python
relative_profile = Atmosphere(
    pressure=[80.0, 50.0, 30.0],
    temperature=[225.0, 220.0, 215.0],
    ozone=[3.5, 5.0, 6.0],
    ozone_unit="ppmv",
)
```

This preserves vertical spacing but does not reconstruct an instrument's
absolute geometric altitude. Supply `altitude` whenever absolute height or
`Box.at_altitude(...)` matters.

## Aerosol controls

Heterogeneous chemistry and radiative aerosol extinction are independent:

| Custom aerosol profile | `heterogeneous` | `aerosol_extinction` | Effect |
|---|---:|---:|---|
| Absent | `True` | `None` | Default per-box heterogeneous surface; no profile extinction |
| Present | `True` | `None` | Heterogeneous chemistry and automatic radiative extinction |
| Present | `False` | `True` | Radiation-only aerosol experiment |
| Present | `True` | `False` | Chemistry-only aerosol experiment |

```python
from pratmo import ChemistryOptions, PhotolysisOptions

radiation_only = model.diurnal(
    atmosphere=atmosphere,
    chemistry=ChemistryOptions(heterogeneous=False),
    photolysis=PhotolysisOptions(aerosol_extinction=True),
)
```

Per-box aerosol supplied through `Box` affects heterogeneous chemistry but is
not a complete radiative column. A radiative calculation requires the
profile-level aerosol field on `Atmosphere`.

## Match observed NO2 by scaling NOy

The high-level constrained workflow iteratively scales total NOy in each box so
modeled NO2 matches a non-negative observed number density at the nearest local
solar time.

```{code-cell} ipython3
:execution_timeout: 600

observed_no2 = [2.0e8, 1.4e8, 8.0e7, 4.0e7, 2.0e7]

matched = model.diurnal_no2_constrained(
    atmosphere=atmosphere,
    latitude=0.0,
    day="2026-04-30",
    observed_no2_cm3=observed_no2,
    target_hhmm=630,
    iterations=3,
)

print(matched.noy_scale)
print(matched.modeled_no2_cm3)
```

```{code-cell} ipython3
import plotly.graph_objects as go

comparison = go.Figure()
comparison.add_bar(x=output.altitude_km, y=observed_no2, name="Observed")
comparison.add_bar(x=output.altitude_km, y=matched.modeled_no2_cm3, name="Modeled after scaling")
comparison.update_layout(
    barmode="group",
    xaxis_title="Altitude (km)",
    yaxis_title="NO2 number density (cm⁻³)",
    template="plotly_white",
)
comparison
```

`target_hhmm` is local solar time from `0000` through `2359`. Report the scale
factors, final mismatch, input uncertainties, and chosen iteration count; a
successful numerical match does not validate the assumed partitioning of NOy.

## Initialize from the legacy climatology

`PratmoClimatology` reads `fort03_LLM.x`, `fort04.x`, `fort05.x`, and
`fort51.x` from an explicit legacy-data directory:

```python
from datetime import date
import numpy as np
from pratmo import Atmosphere, PratmoClimatology

altitude_km = np.array([18.0, 21.0, 24.0])
climatology = PratmoClimatology("path/to/legacy-fortran-data")
sample = climatology.sample(60.0, date(2026, 3, 16), altitude_km)

climatological_atmosphere = Atmosphere(
    pressure=[80.0, 45.0, 24.0],
    temperature=sample.temperature_k,
    ozone=sample.o3,
    ozone_unit="fraction",
    altitude=altitude_km,
    aerosol_surface_area=sample.aerosol_surface_area_um2_cm3,
)

result = model.diurnal(
    atmosphere=climatological_atmosphere,
    initial_mixing_ratios=sample.initial_mixing_ratios(),
)
```

Measured ozone can replace `sample.o3` while retaining climatological N2O,
NOy, and correlated tracer initialization.

## Lower-level compatibility API

`CustomAtmosphereProfile`, `DiurnConfig`, and
`No2ConstrainedDiurnConfig` remain available for code that already serializes
native configurations. They require canonical units. Prefer `Atmosphere` and
`Model` for new code.
