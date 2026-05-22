---
jupytext:
  text_representation:
    format_name: myst
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# DIURN mode

DIURN integrates the photochemistry through a full 24-hour diurnal cycle,
tracking 30 implicit species at 34 time steps (night → sunrise → noon → sunset).

:::{note}
`PratmoModel.with_defaults()` uses embedded science data whose `fort01.x` is
configured for CTM mode (`nd216 > 0`).  Running DIURN with embedded data
raises `ValueError` because the long-lived initial mixing ratios from
`fort02.x` are never loaded, leaving implicit species at zero and causing
Newton-Raphson divergence.

To run DIURN you must supply a data directory that contains a DIURN-mode
`fort01.x` (`nd216 = 0`) **and** a matching `fort02.x`:

```python
model = PratmoModel.from_data_dir("/path/to/fortran/data")
```
:::

## Configuration

```{code-cell} ipython3
from pratmo import DiurnConfig, DiurnBoxSpec, LongLivedMixingRatios

cfg = DiurnConfig(
    latitude_deg=0.0,
    julian_day=120,          # 30 April
    integration_days=20,
    boxes=[
        DiurnBoxSpec(altitude_level=20, albedo=0.05),
        DiurnBoxSpec(altitude_level=25, albedo=0.05),
    ],
    bromine=False,
    solar_flux_scale=1.0,
)
print(cfg)
```

## Supplying initial mixing ratios

Pass one `LongLivedMixingRatios` per box to override the defaults from file:

```{code-cell} ipython3
init = LongLivedMixingRatios(
    o3=5e-6,
    n2o=3.1e-7,
    ch4=1.8e-6,
    h2o=5e-3,
    noy=1e-8,
)
print(f"O₃={init.o3:.1e}  CH₄={init.ch4:.1e}  H₂O={init.h2o:.1e}")
```

To use custom initial conditions in a run, set `initial_mixing_ratios` to a
list with one entry per box:

```python
cfg_with_init = DiurnConfig(
    latitude_deg=0.0,
    julian_day=120,
    integration_days=20,
    boxes=[DiurnBoxSpec(altitude_level=20)],
    initial_mixing_ratios=[init],  # one per box
)
```

## Interpreting DIURN output

When called with a DIURN-mode data directory the output looks like this
(shown here as non-executed code; substitute your own data path):

```python
from pratmo import PratmoModel

model = PratmoModel.from_data_dir("/path/to/data")
out = model.run_diurn(cfg)

# Daily-mean box snapshots (one per box, equivalent to fort08.x)
for snap in out.boxes:
    print(f"{snap.altitude_km:.1f} km  O₃={snap.implicit.o3:.2e}  OH={snap.implicit.oh:.2e}")

# Full 24-hour time series (equivalent to fort07.x)
ts = out.time_series[0]   # first box
for step in ts.steps:
    hh = step.time_hhmm // 100
    mm = step.time_hhmm % 100
    print(f"  {hh:02d}:{mm:02d}  OH={step.implicit.oh:.2e} cm⁻³")

# Species as a numpy array: shape (n_boxes, n_timesteps)
import numpy as np
o3_grid  = out.species_grid("o3")   # O₃ time series for every box
oh_grid  = out.species_grid("oh")
j_no2    = out.jvalue_grid("no2")   # J(NO₂), constant within each day

times = [s.time_hhmm for s in out.time_series[0].steps]
print("Time steps (HHMM):", times)
print("O₃ grid shape:", o3_grid.shape)   # → (2, 34)
```

## Time-step layout

DIURN uses 34 time steps per 24-hour day:

| Phase       | Steps |
|-------------|-------|
| Night       | 6     |
| Sunrise     | 8     |
| Daytime     | 12    |
| Sunset      | 8     |

Time codes are integers in HHMM format (e.g. `1430` = 14:30 UTC).
J-values are recomputed at each solar zenith angle; nighttime steps carry
`J = 0` for photolysis channels.
