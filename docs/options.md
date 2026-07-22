---
jupytext:
  text_representation:
    format_name: myst
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# Model options by example

The high-level API groups switches by what they affect. Defaults represent an
ordinary, runnable stratospheric calculation: established bromine chemistry is
on, experimental iodine is off, heterogeneous chemistry is on, solar flux is
unscaled, lower-boundary albedo is 0.20, and integrations are long enough for a
useful first equilibrium attempt.

{doc}`defaults` provides a compact table of every high-level default and its
rationale.

## Chemistry mechanisms

```python
from pratmo import ChemistryOptions, Model

model = Model()

# Recommended default
standard = ChemistryOptions()

# Gas-phase only: no sulfate or sea-salt heterogeneous reactions
gas_phase_only = ChemistryOptions(heterogeneous=False)

# Experimental iodine extension (emits ExperimentalFeatureWarning)
with_iodine = ChemistryOptions(iodine=True)

# Mechanism isolation experiment
no_halogen_families = ChemistryOptions(bromine=False, iodine=False)

output = model.diurnal(chemistry=gas_phase_only)
```

Iodine is opt-in because the extension is not part of the Fortran reference
and has not been scientifically validated. Disabling a family removes its
active chemistry; it does not turn this implementation into an independently
validated model.

Compare heterogeneous and gas-phase-only chemistry at the same box:

```{code-cell} ipython3
:execution_timeout: 600

import plotly.graph_objects as go
from pratmo import Box, ChemistryOptions, Model

comparison_model = Model()
standard_cycle = comparison_model.diurnal(boxes=[Box.at_level(15)])
gas_only_cycle = comparison_model.diurnal(
    boxes=[Box.at_level(15)],
    chemistry=ChemistryOptions(heterogeneous=False),
)

chemistry_figure = go.Figure()
hours = standard_cycle.elapsed_seconds / 3600.0
chemistry_figure.add_scatter(
    x=hours, y=standard_cycle.species_grid("hno3")[0], name="Default heterogeneous"
)
chemistry_figure.add_scatter(
    x=hours, y=gas_only_cycle.species_grid("hno3")[0], name="Gas phase only"
)
chemistry_figure.update_layout(
    xaxis_title="Elapsed hours from local noon",
    yaxis_title="HNO3 number density (cm⁻³)",
    yaxis_type="log",
    template="plotly_white",
)
chemistry_figure
```

## Solar flux, albedo, and aerosol extinction

```python
from pratmo import PhotolysisOptions

# Defaults
nominal = PhotolysisOptions()

# Bright lower boundary
snow = PhotolysisOptions(surface_albedo=0.80)

# Scale the top-of-atmosphere solar spectrum by 2%
solar_sensitivity = PhotolysisOptions(solar_flux_scale=1.02)

# Explicitly include a custom Atmosphere's aerosol profile in radiation
with_aerosol_radiation = PhotolysisOptions(aerosol_extinction=True)

output = model.diurnal(photolysis=snow)
```

`aerosol_extinction=None` is the automatic default: it activates when a custom
`Atmosphere` contains a nonzero aerosol profile. `True` without such a profile
emits a warning. Per-box aerosol values still control heterogeneous chemistry
but cannot define the full radiative column.

The effect of an albedo choice is best inspected as a paired experiment:

```{code-cell} ipython3
:execution_timeout: 600

from pratmo import PhotolysisOptions

bright_cycle = comparison_model.diurnal(
    boxes=[Box.at_level(15)],
    photolysis=PhotolysisOptions(surface_albedo=0.80),
)

albedo_figure = go.Figure()
albedo_figure.add_scatter(
    x=hours, y=standard_cycle.species_grid("oh")[0], name="Albedo 0.20"
)
albedo_figure.add_scatter(
    x=hours, y=bright_cycle.species_grid("oh")[0], name="Albedo 0.80"
)
albedo_figure.update_layout(
    xaxis_title="Elapsed hours from local noon",
    yaxis_title="OH number density (cm⁻³)",
    yaxis_type="log",
    template="plotly_white",
)
albedo_figure
```

## Integration length and parallel boxes

```python
from pratmo import CtmOptions, DiurnalOptions

ctm_long = CtmOptions(integration_days=60)

diurn_parallel = DiurnalOptions(
    integration_days=30,
    parallel_boxes=True,
)
```

`Model.ctm()` defaults to 40 integration days. `Model.diurnal()` defaults to
20. Runs shorter than three days emit a warning because their endpoints are
unlikely to represent photochemical equilibrium. The final authority is the
`output.diagnostics` object, not the requested day count.

`parallel_boxes=None` automatically enables box-level parallelism when more
than one box is requested. Set it explicitly when reproducible scheduling or a
single-threaded embedding matters. Parallel evaluation does not change output
ordering: rows remain in the requested box order.

## Dates and latitude

All of these describe the same seasonal input:

```python
from datetime import date

model.ctm(latitude=52.1, day=172)
model.ctm(latitude=52.1, day="2026-06-21")
model.ctm(latitude=52.1, day=date(2026, 6, 21))
```

Latitude is geodetic degrees from -90 to +90. The calendar year is used only
to derive day of year; the model does not include year-specific emissions or
meteorology. CTM bins 52.1° to its lower 50° climatology grid point and warns;
DIURN evaluates the solar geometry at 52.1° exactly.

## Standard boxes

```python
from pratmo import Box

# A standard pressure level
box = Box.at_level(15)

# Add sulfate and sea-salt surface area for heterogeneous chemistry
aerosol_box = Box.at_level(
    15,
    aerosol_surface_area_um2_cm3=0.8,
    sea_salt_surface_area_um2_cm3=0.03,
)

# Temperature sensitivity around the standard atmosphere
warm_box = Box.at_level(15, temperature_offset_k=5.0)
```

The standard grid is 1-based, from level 1 near the surface through level 41
at the top. CTM uses only standard levels. DIURN can additionally use exact
altitudes in a custom atmosphere; see {doc}`atmospheric-inputs`.

## Explicit time coordinates

DIURN normally generates the 34-point legacy noon-to-noon grid. Supply elapsed
hours to integrate and evaluate photolysis at exact archived times:

```python
import numpy as np

half_hourly = DiurnalOptions(
    elapsed_time_hours=np.arange(0.0, 24.5, 0.5),
)
output = model.diurnal(options=half_hourly)
```

The coordinate must contain 2–64 strictly increasing values, begin at 0, and
end at 24. Use `output.elapsed_seconds` for plots. `output.time_hhmm` is a
cyclic clock label, so sorting it would splice the trajectory incorrectly.

## C++ compatibility mode

```python
cpp_comparison = DiurnalOptions(cpp_compatibility=True)
output = model.diurnal(options=cpp_comparison)
```

This selects the later C++ box model's fixed-cos(SZA) time grid and whole-day
endpoint convergence test. It is a comparison policy, not a generally more
accurate setting. An explicit `elapsed_time_hours` grid takes precedence.

Compatibility controls are for reproducing or comparing archived workflows;
they should not be enabled merely because a run is difficult to converge.

## Lower-level configuration objects

The original typed extension API remains public for workflows that serialize
or mutate configuration objects directly:

```python
from pratmo import CtmBoxSpec, CtmConfig, PratmoModel

native_model = PratmoModel()  # embedded data by default
config = CtmConfig()          # runnable four-level profile by default
config.latitude_deg = 60.0
config.boxes = [CtmBoxSpec(altitude_level=20)]
output = native_model.run_ctm(config)
```

Use this layer when exact compatibility with existing code matters. Prefer
`Model` for explicit units, dates, automatic parallel defaults, and warnings.
