# Custom atmospheres and observed NO2

DIURN can run on a caller-provided pressure, temperature, altitude, and ozone
profile. This is the preferred interface for instrument grids and retrieval
studies. Pressure must decrease with index, altitude must increase, and every
array must describe the same number of levels.

## Run every custom level

Leave `boxes` empty to create one chemistry box per atmosphere level. PRATMO
supports at most 25 boxes in one run.

```python
from pratmo import CustomAtmosphereProfile, DiurnConfig, PratmoModel

atmosphere = CustomAtmosphereProfile(
    pressure_mb=[80.0, 50.0, 30.0],
    temperature_k=[225.0, 220.0, 215.0],
    altitude_km=[18.0, 21.0, 24.0],
    o3=[3.5e-6, 5.0e-6, 6.0e-6],
    o3_kind="mixing_ratio",
)

model = PratmoModel.with_defaults()
output = model.run_diurn(
    DiurnConfig(
        latitude_deg=0.0,
        julian_day=120,
        integration_days=5,
        atmosphere=atmosphere,
        boxes=[],
        bromine=True,
        iodine=True,
        parallel_boxes=True,
    )
)

print(output.altitude_km)
print(output.species_grid("no2").shape)
```

Set `o3_kind="number_density"` when `o3` is already in cm⁻³. If
`altitude_km` is omitted, PRATMO estimates altitude hydrostatically from the
pressure and temperature profile.

## Select custom levels

When an atmosphere contains radiative-transfer levels that should not all be
chemistry boxes, pass explicit `DiurnBoxSpec` objects. Their `altitude_level`
values are 1-based indices into the custom profile.

```python
from pratmo import DiurnBoxSpec

config = DiurnConfig(
    atmosphere=atmosphere,
    boxes=[DiurnBoxSpec(altitude_level=1), DiurnBoxSpec(altitude_level=3)],
)
```

## Constrain total NOy with observed NO2

`run_diurn_no2_constrained` iteratively rescales each box's total NOy so that
modeled NO2 matches an observation at the nearest local-solar-time step.

```python
from pratmo import No2ConstrainedDiurnConfig

result = model.run_diurn_no2_constrained(
    No2ConstrainedDiurnConfig(
        diurn=DiurnConfig(
            latitude_deg=0.0,
            julian_day=120,
            integration_days=3,
            atmosphere=atmosphere,
            boxes=[],
        ),
        observed_no2_cm3=[2.0e8, 5.0e7, 1.0e7],
        target_hhmm=630,
        iterations=3,
    )
)

print(result.noy_scale)
print(result.modeled_no2_cm3)
print(result.output.species_grid("no2"))
```

`target_hhmm` is a local solar time from `0000` through `2359`. Nearest-step
selection is cyclic across midnight. The observation array must contain one
finite, non-negative number density per chemistry box.

See [the compact custom-atmosphere example](../examples/no2_constrained_custom_atmosphere.py)
and [the OMPS NetCDF batch workflow](../examples/omps_no2_constrained_batch.py)
for complete scripts.
