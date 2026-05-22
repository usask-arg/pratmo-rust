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
The embedded science data includes climatological initial conditions (`fort02.x`),
so no external files are needed.

## Basic run

```{code-cell} ipython3
from pratmo import PratmoModel, DiurnConfig, DiurnBoxSpec

model = PratmoModel.with_defaults()
cfg = DiurnConfig(
    latitude_deg=0.0,
    julian_day=120,          # 30 April
    integration_days=20,
    boxes=[DiurnBoxSpec(altitude_level=20)],
)
out = model.run_diurn(cfg)
snap = out.boxes[0]
print(f"{snap.altitude_km:.1f} km  O₃={snap.implicit.o3:.2e} cm⁻³  OH={snap.implicit.oh:.2e} cm⁻³")
```

## Diurnal time series

```{code-cell} ipython3
ts = out.time_series[0]
print(f"Time steps: {len(ts.steps)}")
print(f"\n{'HHMM':>6}  {'OH (cm⁻³)':>12}  {'O₃ (cm⁻³)':>12}  {'NO₂ (cm⁻³)':>12}")
print("-" * 48)
for step in ts.steps[::4]:   # every 4th step
    hh = step.time_hhmm // 100
    mm = step.time_hhmm % 100
    print(f"  {hh:02d}:{mm:02d}  {step.implicit.oh:12.2e}  {step.implicit.o3:12.2e}  {step.implicit.no2:12.2e}")
```

## Species grids as numpy arrays

`species_grid` returns shape `(n_boxes, n_timesteps)` — one row per box:

```{code-cell} ipython3
import numpy as np

cfg2 = DiurnConfig(
    latitude_deg=0.0,
    julian_day=120,
    integration_days=20,
    boxes=[DiurnBoxSpec(altitude_level=15), DiurnBoxSpec(altitude_level=20)],
)
out2 = model.run_diurn(cfg2)

o3  = out2.species_grid("o3")   # shape (2, 34)
oh  = out2.species_grid("oh")
times = [s.time_hhmm for s in out2.time_series[0].steps]

print(f"Grid shape: {o3.shape}  (boxes × timesteps)")
print(f"\nPeak daytime OH (cm⁻³):")
for i, snap in enumerate(out2.boxes):
    print(f"  {snap.altitude_km:.0f} km: {oh[i].max():.2e}")
```

## Supplying custom initial mixing ratios

Pass one `LongLivedMixingRatios` per box to override the embedded defaults.
The values must be physically consistent with the target altitude and season;
the model will iterate to convergence from there.

```python
from pratmo import LongLivedMixingRatios

init = LongLivedMixingRatios(
    o3=5e-6,
    n2o=3.1e-7,
    ch4=1.8e-6,
    h2o=5e-3,
    noy=1e-8,
)

cfg_custom = DiurnConfig(
    latitude_deg=0.0,
    julian_day=120,
    integration_days=20,
    boxes=[DiurnBoxSpec(altitude_level=20)],
    initial_mixing_ratios=[init],  # one entry per box
)
out_custom = model.run_diurn(cfg_custom)
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
