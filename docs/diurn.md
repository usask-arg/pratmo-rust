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

DIURN integrates a full noon-to-noon photochemical cycle and retains the 34
resolved time steps. Use it when local solar time matters or when the
atmosphere comes from an instrument or retrieval.

## Basic run

```{code-cell} ipython3
:execution_timeout: 600

from pratmo import Box, Model

model = Model()
cycle = model.diurnal(
    latitude=0.0,
    day="2026-04-30",
    boxes=[Box.at_level(15)],
)

print(cycle)
print(cycle.altitude_km)
print(cycle.species_grid("oh").shape)
```

The default atmosphere and initial conditions are embedded. `Box.at_level(15)`
is a standard atmospheric level; see {doc}`standard-levels`.

## Time coordinates

`elapsed_seconds` increases monotonically from zero to 86,400 seconds.
`time_hhmm` is a cyclic local-solar-time label, so both endpoints are 1200.

```{code-cell} ipython3
for elapsed, hhmm in zip(cycle.elapsed_seconds[::4], cycle.time_hhmm[::4]):
    print(f"{elapsed / 3600:5.1f} elapsed hours  {hhmm:04d} local solar time")
```

Plot and sort with elapsed time. Use `time_hhmm` only as a clock label.

## Interactive time series

```{code-cell} ipython3
from pratmo.plotting import plot_diurnal

plot_diurnal(cycle, ["oh", "ho2", "no", "no2"])
```

## Multiple boxes

```{code-cell} ipython3
:execution_timeout: 600

from pratmo import DiurnalOptions

profile_cycle = model.diurnal(
    latitude=0.0,
    day="2026-04-30",
    boxes=[Box.at_level(level) for level in (8, 12, 18, 22)],
    options=DiurnalOptions(integration_days=20),
)

print(profile_cycle.species_grid("o3").shape)  # (box, time)
```

```{code-cell} ipython3
from pratmo.plotting import plot_profile

plot_profile(profile_cycle, ["o3", "no2", "hno3"])
```

With `parallel_boxes=None`, multiple boxes run in parallel automatically. Their
output order remains the requested box order.

## Exact time coordinates

Supply elapsed hours when matching an external sampling grid:

```python
import numpy as np
from pratmo import DiurnalOptions

half_hourly = DiurnalOptions(
    elapsed_time_hours=np.arange(0.0, 24.5, 0.5),
)
cycle = model.diurnal(options=half_hourly)
```

The coordinate must contain 2–64 strictly increasing values, start at 0, and
end at 24. Every value is both an integration and a photolysis-evaluation time.

## Custom initial mixing ratios

Pass one complete long-lived state per chemistry box. The helper supplies
representative lower-stratospheric values so unspecified gases are not
silently set to zero.

```python
from pratmo import background_mixing_ratios, ppbv, ppmv

initial = background_mixing_ratios(
    o3=ppmv(5.0),
    n2o=ppbv(300),
    noy=ppbv(10),
)

cycle = model.diurnal(
    boxes=[Box.at_level(15)],
    initial_mixing_ratios=[initial],
)
```

These values are examples, not a climatology. Quantitative work should use
altitude-, latitude-, and season-appropriate measured or climatological inputs.

## J-value limitation

`jvalue_profile` exposes daily-mean photolysis frequencies. `jvalue_grid`
repeats those means across the time dimension for array-shape compatibility;
it is not a resolved actinic-flux time series.

## Lower-level compatibility API

The extension-level `PratmoModel` and `DiurnConfig` classes remain available
for serialized legacy configurations and exact compatibility work. They use
canonical units and fewer user-facing checks. New workflows should start with
`Model.diurnal()`.
