---
jupytext:
  text_representation:
    format_name: myst
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# Quick start tutorial

This tutorial introduces the model, runs both main workflows, explains the
output, and creates interactive plots. It uses the high-level Python interface,
which is the best starting point for new code.

## 1. What the model does

PRATMO is a one-dimensional collection of independent atmospheric **boxes**.
Each box has pressure, temperature, altitude, long-lived chemical families,
and 40 short-lived species. It computes photolysis and gas-phase and
heterogeneous chemistry through a full solar day.

There are two workflows:

- **CTM** integrates selected standard atmospheric levels toward a
  climatological photochemical steady state. Use it for vertical profiles and
  parameter sweeps.
- **DIURN** retains the resolved noon-to-noon trajectory. Use it for local-time
  behavior, observations made at a particular time, or a custom atmosphere.

Both workflows treat boxes independently; this is a photochemical box model,
not a transport model. Values are representative only after you have checked
the atmosphere, initialization, convergence diagnostics, and scientific
validity for your application. {doc}`concepts` explains the model structure in
more detail.

## 2. Install and import

Install the package with interactive plotting support:

```bash
python -m pip install "pratmo[plot]"
```

Create a model. Embedded science data is used automatically, so no external
input directory is required.

```{code-cell} ipython3
from pratmo import Model

model = Model()
```

## 3. Run a vertical profile

`ctm()` has runnable defaults: 45°N near the March equinox, four standard
levels, established bromine chemistry, experimental iodine disabled, and a
40-day integration.

```{code-cell} ipython3
:execution_timeout: 600

profile = model.ctm()
print(profile)
print("altitude (km):", profile.altitude_km)
```

The output provides NumPy arrays. Short-lived species are number densities in
molecules cm⁻³; long-lived species and families are dimensionless mixing
ratios; J-values are photolysis rates in s⁻¹.

```{code-cell} ipython3
o3_density = profile.species_profile("o3")
noy_fraction = profile.long_lived_profile("noy")
j_no2 = profile.jvalue_profile("no2")

for altitude_km, ozone, noy, jvalue in zip(
    profile.altitude_km, o3_density, noy_fraction, j_no2
):
    print(
        f"{altitude_km:5.1f} km  O3={ozone:.3e} cm-3  "
        f"NOy={noy:.3e}  J(NO2)={jvalue:.3e} s-1"
    )
```

For this illuminated stratospheric profile, ozone and J(NO2) should be finite
and positive. A zero photolysis profile, negative concentration, or non-finite
value is a reason to stop and investigate rather than continue plotting.

Plotly figures remain interactive in the built documentation: hover for exact
values, zoom, pan, and toggle traces in the legend.

```{code-cell} ipython3
from pratmo.plotting import plot_profile

plot_profile(profile, ["o3", "no2", "oh"])
```

Mixing-ratio and photolysis profiles use the same helper.

```{code-cell} ipython3
plot_profile(
    profile,
    ["o3", "ch4", "h2o"],
    kind="mixing_ratio",
    unit="ppmv",
)
```

NOy is normally easier to read in ppbv:

```{code-cell} ipython3
plot_profile(profile, "noy", kind="mixing_ratio", unit="ppbv")
```

## 4. Choose location, date, and boxes

Dates can be ISO strings, `datetime.date` objects, or day-of-year integers.
Standard levels are 1-based. Exact geometric altitudes are available with a
custom atmosphere. The orientation table in {doc}`standard-levels` maps common
level numbers to approximate heights and pressures.

```{code-cell} ipython3
from datetime import date
from pratmo import Box, CtmOptions

southern_summer = model.ctm(
    latitude=-35.0,
    day=date(2026, 1, 15),
    boxes=[Box.at_level(level) for level in (8, 12, 16, 20, 24)],
    options=CtmOptions(integration_days=30),
)
print(southern_summer.altitude_km)
```

## 5. Run a resolved diurnal cycle

The default DIURN run uses one mid-stratospheric standard box. Plot against
`elapsed_seconds`, the monotonic noon-to-noon coordinate; `time_hhmm` is a
cyclic local-solar-time label and both endpoints are noon.

```{code-cell} ipython3
:execution_timeout: 600

cycle = model.diurnal(latitude=0.0, day="2026-04-30")
print(cycle)
print("steps:", len(cycle.time_series[0]))
print("clock endpoints:", cycle.time_hhmm[0], cycle.time_hhmm[-1])
```

```{code-cell} ipython3
from pratmo.plotting import plot_diurnal

plot_diurnal(cycle, ["oh", "ho2", "no2"])
```

## 6. Use readable quantities

PRATMO stores mixing ratios as dimensionless fractions. The unit helpers make
input intent explicit without adding a units-package dependency.

```{code-cell} ipython3
from pratmo import ppbv, ppmv, pptv, background_mixing_ratios

initial = background_mixing_ratios(
    o3=ppmv(5.2),
    n2o=ppbv(285),
    noy=ppbv(9.5),
    iodx=pptv(1.0),
)

print(initial.o3, initial.n2o, initial.noy, initial.iodx)
```

See {doc}`atmospheric-inputs` for pressure, temperature, altitude, ozone number
density, aerosol, and complete custom-profile examples.

## 7. Understand warnings and diagnostics

Hard inconsistencies—such as increasing pressure with altitude, mismatched
array lengths, or negative concentrations—raise `ValueError`. Inputs that are
technically possible but look suspicious emit `PratmoWarning`. Examples include
a one-day equilibrium run, an unusually large temperature offset, ozone that
looks like unconverted ppmv, or a very high surface albedo.

```python
from pratmo import CtmOptions, PratmoWarning
import warnings

with warnings.catch_warnings():
    warnings.simplefilter("error", PratmoWarning)
    model.ctm(options=CtmOptions(integration_days=1))
```

Every output also carries numerical diagnostics:

```{code-cell} ipython3
diagnostics = profile.diagnostics
print("Newton non-convergence count:", diagnostics.newraf_nonconvergence_count)
print("Daily non-convergence count:", diagnostics.rafday_nonconvergence_count)
print("Maximum final relative correction:", diagnostics.rafday_max_final_relative_correction)
```

Warnings and successful solver convergence do not establish scientific
validity. They are guardrails for common interface mistakes.

## 8. Where to go next

- {doc}`atmospheric-inputs` covers every atmospheric input and unit.
- {doc}`options` covers chemistry, radiation, integration, and compatibility
  switches with examples.
- {doc}`diurn` goes deeper into time grids and array layouts.
- {doc}`explorer` interactively compares every implicit species and chemical
  family.
- {doc}`results` covers convergence checks and export to pandas or xarray.
- {doc}`science-scope` states what the box model does and does not represent.
