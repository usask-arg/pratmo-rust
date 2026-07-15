---
jupytext:
  text_representation:
    format_name: myst
kernelspec:
  display_name: Python 3
  language: python
  name: python3
---

# Quickstart: CTM mode

The CTM (climatological transport model) mode integrates a set of altitude
boxes to a photochemical steady state. It runs entirely on compiled-in
science data — no external files needed.

## Basic run

```{code-cell} ipython3
from pratmo import IMPLICIT_SPECIES_NAMES, PratmoModel, CtmConfig, CtmBoxSpec

model = PratmoModel.with_defaults()

cfg = CtmConfig(
    latitude_deg=60.0,
    julian_day=75,        # 16 March
    integration_days=40,
    boxes=[CtmBoxSpec(altitude_level=lvl) for lvl in [10, 15, 20, 25]],
)

out = model.run_ctm(cfg)
print(f"Ran {len(out.boxes)} boxes")
for snap in out.boxes:
    print(f"  {snap.altitude_km:.1f} km  O₃={snap.implicit.o3:.2e} cm⁻³  OH={snap.implicit.oh:.2e} cm⁻³")
```

## Species profiles as numpy arrays

`CtmOutput.species_profile` returns a 1-D numpy array with one value per box.

```{code-cell} ipython3
import numpy as np

o3  = out.species_profile("o3")
oh  = out.species_profile("oh")
no2 = out.species_profile("no2")
alts = out.altitude_km

print("altitude (km) | O3 (cm⁻³)   | OH (cm⁻³)   | NO2 (cm⁻³)")
print("-" * 60)
for z, o, h, n in zip(alts, o3, oh, no2):
    print(f"  {z:5.1f}       | {o:.3e}  | {h:.3e}  | {n:.3e}")
```

## J-value profiles

Photolysis rates are returned the same way via `jvalue_profile`.

```{code-cell} ipython3
j_no2 = out.jvalue_profile("no2")
j_o3d = out.jvalue_profile("o3_o1d")

print("altitude (km) | J(NO₂) s⁻¹  | J(O³→O¹D) s⁻¹")
print("-" * 50)
for z, jn, jo in zip(alts, j_no2, j_o3d):
    print(f"  {z:5.1f}       | {jn:.3e}   | {jo:.3e}")
```

## Accessing long-lived species and box metadata

Each `BoxSnapshot` carries atmospheric state alongside the chemistry. For
array workflows, `long_lived_profile` extracts any name in
`pratmo.LONG_LIVED_NAMES` without a Python loop.

```{code-cell} ipython3
snap = out.boxes[-1]   # highest box
print(f"Box {snap.box_index}:  {snap.altitude_km:.1f} km,  {snap.pressure_mb:.2f} mb,  {snap.temperature_k:.1f} K")
print(f"  Air density : {snap.air_density_cm3:.3e} cm⁻³")
print(f"  O₃ (long-lived MR): {snap.long_lived.o3:.3e}")
print(f"  CH₄          :      {snap.long_lived.ch4:.3e}")
print(f"  H₂O          :      {snap.long_lived.h2o:.3e}")

noy = out.long_lived_profile("noy")
ch4 = out.long_lived_profile("ch4")
print(f"Profile shapes: NOy={noy.shape}, CH4={ch4.shape}")
```

## Discovering available fields

The three exported name tuples are the source of truth for string-based array
extraction. Names passed to profile and grid methods are case-insensitive.

```{code-cell} ipython3
from pratmo import JVALUE_NAMES, LONG_LIVED_NAMES

print(f"Implicit species ({len(IMPLICIT_SPECIES_NAMES)}):", IMPLICIT_SPECIES_NAMES)
print(f"Long-lived fields ({len(LONG_LIVED_NAMES)}):", LONG_LIVED_NAMES)
print(f"J-values ({len(JVALUE_NAMES)}):", JVALUE_NAMES)
assert np.array_equal(out.species_profile("O3"), out.species_profile("o3"))
```

## Sweeping latitude

Run the same set of boxes at multiple latitudes to build a latitude–altitude grid.

```{code-cell} ipython3
latitudes = [-60, -30, 0, 30, 60]
oh_grid = []   # shape: (n_lat, n_box)

for lat in latitudes:
    cfg_lat = CtmConfig(
        latitude_deg=float(lat),
        julian_day=180,       # June solstice
        integration_days=20,
        boxes=[CtmBoxSpec(altitude_level=lvl) for lvl in [10, 15, 20, 25]],
    )
    result = model.run_ctm(cfg_lat)
    oh_grid.append(result.species_profile("oh"))

oh_grid = np.array(oh_grid)
print("OH (cm⁻³) — rows=latitude, cols=altitude")
print("Lat\\Alt  ", "  ".join(f"{a:6.1f}km" for a in alts))
for lat, row in zip(latitudes, oh_grid):
    print(f"  {lat:+4d}°   ", "  ".join(f"{v:.2e}" for v in row))
```
