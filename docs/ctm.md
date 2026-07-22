---
jupytext:
  text_representation:
    format_name: myst
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# CTM vertical profiles

CTM combines the embedded seasonal climatology with selected standard levels
and integrates each box toward a repeatable daily chemical state.

```{code-cell} ipython3
:execution_timeout: 600

from pratmo import Box, CtmOptions, Model, mixing_ratio_as

model = Model()
profile = model.ctm(
    latitude=60.0,
    day="2026-03-16",
    boxes=[Box.at_level(level) for level in (10, 15, 20, 25)],
    options=CtmOptions(integration_days=40),
)

print(profile.altitude_km)
print(mixing_ratio_as(profile.long_lived_profile("noy"), "ppbv"))
```

Use the coordinates returned by the run because standard-level altitudes vary
with climatological temperature. See {doc}`standard-levels` for orientation.
The CTM climatology is tabulated on a 2.5° grid from 87.5°S to 87.5°N. A
latitude between grid points uses the lower grid point, matching the original
workflow, and emits `PratmoWarning`; DIURN uses the requested latitude exactly.

## Plot chemistry and photolysis

```{code-cell} ipython3
from pratmo.plotting import plot_profile

plot_profile(profile, ["o3", "no2", "oh"])
```

```{code-cell} ipython3
plot_profile(profile, ["no2", "o3_o1d"], kind="jvalue")
```

## Check convergence

```{code-cell} ipython3
d = profile.diagnostics
print("Newton failures:", d.newraf_nonconvergence_count)
print("Daily failures:", d.rafday_nonconvergence_count)
print("Largest final relative correction:", d.rafday_max_final_relative_correction)
```

Forty days is a starting default, not a universal guarantee. Repeat important
calculations with a longer integration and compare the target quantities as
described in {doc}`results`.
