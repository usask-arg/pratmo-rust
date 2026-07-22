# Results, plotting, and export

CTM returns one snapshot per box. DIURN adds a shared noon-to-noon coordinate
and arrays shaped `(box, time)`.

## Array shapes and units

```python
profile.species_profile("o3")       # (box,), molecules cm-3
profile.long_lived_profile("noy")  # (box,), fraction
profile.jvalue_profile("no2")      # (box,), s-1

cycle.species_grid("oh")           # (box, time), molecules cm-3
cycle.elapsed_seconds               # (time,), monotonic 0..86400
cycle.time_hhmm                     # (time,), cyclic local-time labels
```

Use `elapsed_seconds` as the plotting and sorting coordinate. Both endpoints of
a DIURN run are labelled 12:00 local solar time, so `time_hhmm` is not
monotonic.

## Interpret diagnostics

| Field | What to check |
|---|---|
| `newraf_nonconvergence_count` | Zero is expected for an uncomplicated run |
| `rafday_nonconvergence_count` | Nonzero values require investigation |
| `rafday_max_final_relative_correction` | Smaller is better; compare it with the precision required by the study |
| `rafday_max_correction_iterations` | Large values indicate difficult daily relaxation |

Numerical convergence is necessary but does not validate the atmosphere or
chemistry. Also compare a normal run with a longer integration:

```python
import numpy as np
from pratmo import CtmOptions

normal = model.ctm(options=CtmOptions(integration_days=40))
longer = model.ctm(options=CtmOptions(integration_days=60))

relative_change = np.abs(
    longer.species_profile("o3") / normal.species_profile("o3") - 1.0
)
print(relative_change)
```

Choose acceptance thresholds from the scientific question; the interface does
not invent one.

## Plotting

Install `pratmo[plot]`, then use:

```python
from pratmo.plotting import plot_atmosphere, plot_diurnal, plot_profile

plot_profile(profile, ["o3", "no2", "oh"]).show()
plot_diurnal(cycle, ["oh", "ho2", "no2"], box=0).show()
plot_atmosphere(atmosphere).show()
```

Every interactive figure should be read together with its units and selected
box. Logarithmic axes make trace shapes comparable but can visually exaggerate
very small concentrations.

## Pandas and xarray

The core API returns NumPy arrays and does not require pandas or xarray.

```python
import pandas as pd

table = pd.DataFrame({
    "altitude_km": profile.altitude_km,
    "pressure_hpa": profile.pressure_mb,
    "o3_cm3": profile.species_profile("o3"),
    "noy_fraction": profile.long_lived_profile("noy"),
})
table.to_csv("pratmo-profile.csv", index=False)
```

```python
import xarray as xr

dataset = xr.Dataset(
    data_vars={
        "no2": (("box", "time"), cycle.species_grid("no2")),
        "oh": (("box", "time"), cycle.species_grid("oh")),
    },
    coords={
        "box": range(len(cycle)),
        "altitude_km": ("box", cycle.altitude_km),
        "elapsed_seconds": ("time", cycle.elapsed_seconds),
        "local_time_hhmm": ("time", cycle.time_hhmm),
    },
)
dataset.to_netcdf("pratmo-diurnal.nc")
```

Install `pratmo[io]` for the xarray and NetCDF dependencies.
