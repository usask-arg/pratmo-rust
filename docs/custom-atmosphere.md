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
    aerosol_surface_area_um2_cm3=[0.35, 0.20, 0.08],
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

`DiurnConfig.surface_albedo` sets the Lambertian lower-boundary reflectivity
used by the photolysis calculation (default `0.20`). Set
`heterogeneous_chemistry=False` to disable all sulfate-aerosol and sea-salt
heterogeneous reaction channels. This switch is independent of the per-box
aerosol surface-area values.

Set `radiative_aerosol=True` to include the custom sulfate aerosol profile in
photolysis optical depth and scattering. The radiative conversion uses the
original PRATMO 300-nm surface-area-to-extinction factor and wavelength law.
It is independent of `heterogeneous_chemistry`, so chemistry and radiative
aerosol effects can be tested separately.

Set `o3_kind="number_density"` when `o3` is already in cm⁻³. If
`altitude_km` is omitted, PRATMO estimates altitude hydrostatically from the
pressure and temperature profile.

Custom DIURN atmospheres may contain up to 81 radiative levels, matching the
later C++ model's default 0–80 km one-kilometre shell grid. The embedded legacy
Fortran atmosphere and CTM inputs remain 41-level grids. Set
`cpp_compatibility=True` to also use the C++ fixed-cos(SZA) time grid and its
whole-day endpoint convergence test.

For comparisons against an externally archived time coordinate, pass
`elapsed_time_hours=[0.0, ..., 24.0]`. The values must be strictly increasing,
start at local noon (0 h), end at the following local noon (24 h), and contain
at most 64 points. Each supplied point is used as both an integration and a
photolysis-evaluation time, overriding the generated legacy or C++ grid.

Structured API runs follow the C++ default of zero rainout. Rainout is distinct
from heterogeneous aerosol chemistry: sulfate reactions remain enabled by
default whenever a nonzero aerosol surface area is supplied.

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

For a chemistry box between radiative shells, also set `altitude_km`. PRATMO
then interpolates pressure, temperature, O3, and aerosol to the exact chemistry
altitude and linearly interpolates each wavelength's actinic flux between its
two surrounding radiative shells, matching the later C++ box model's layout.

```python
config = DiurnConfig(
    atmosphere=atmosphere,
    boxes=[DiurnBoxSpec(altitude_level=2, altitude_km=22.5)],
    radiative_aerosol=True,
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

## PRATMO climatological initialization

`PratmoClimatology(data_dir)` reads the original `fort03_LLM.x`, `fort04.x`,
`fort05.x`, and `fort51.x` monthly tables. Its `sample()` method interpolates
temperature, O3, N2O, and NOy in date, latitude, and altitude using the later
PRATMO C++ climatology rules, and samples the seasonal Thomason aerosol surface
area. `profile.initial_mixing_ratios()` derives CH4, H2O, Bry, Cly, and CH3Cl
from N2O and supplies the constant background gases expected by PRATMO.

Measured profiles can take precedence selectively. For example, pass measured
O3 to `initial_mixing_ratios(o3=measured_o3)` while retaining climatological N2O,
NOy, and correlated tracers. Aerosol values can be attached to `DiurnBoxSpec`
even when heterogeneous chemistry is disabled, making that chemistry switch
explicit and reversible.

The legacy climatology predates the iodine mechanism and therefore has no Iy
field. The profile builder retains PRATMO-Rust's 1 ppt background Iy value;
using zero for this absent field makes the iodine-family Newton system singular.
