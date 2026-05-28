"""Run DIURN on a custom atmosphere and scale NOy to match observed NO2.

This example shows the workflow for an instrument profile:

1. Provide pressure, temperature, and O3 on your own vertical grid.
2. Provide observed NO2 at a chosen local solar time.
3. Run an iterative loop that scales total NOy in each box by
   observed_NO2 / modeled_NO2 at that local solar time.
4. Read the final species as arrays shaped (box, LST step).

The arrays below are intentionally tiny. Replace them with your instrument or
retrieval grid.
"""

from __future__ import annotations

import numpy as np

from pratmo import (
    CustomAtmosphereProfile,
    DiurnConfig,
    No2ConstrainedDiurnConfig,
    PratmoModel,
)


def main() -> None:
    # Pressure must decrease with altitude. O3 is supplied here as a dimensionless
    # volume mixing ratio; use o3_kind="number_density" for cm^-3 instead.
    pressure_mb = [80.0, 50.0, 30.0]
    temperature_k = [225.0, 220.0, 215.0]
    altitude_km = [18.0, 21.0, 24.0]
    o3_vmr = [3.5e-6, 5.0e-6, 6.0e-6]

    observed_no2_cm3 = [2.0e8, 5.0e7, 1.0e7]
    target_hhmm = 630

    atmosphere = CustomAtmosphereProfile(
        pressure_mb=pressure_mb,
        temperature_k=temperature_k,
        o3=o3_vmr,
        o3_kind="mixing_ratio",
        altitude_km=altitude_km,
    )

    diurn = DiurnConfig(
        latitude_deg=0.0,
        julian_day=120,
        integration_days=3,
        boxes=[],  # empty => one DIURN box per custom atmosphere level
        bromine=True,
        iodine=False,
        parallel_boxes=True,
        atmosphere=atmosphere,
    )

    constrained = No2ConstrainedDiurnConfig(
        diurn=diurn,
        observed_no2_cm3=observed_no2_cm3,
        target_hhmm=target_hhmm,
        iterations=3,
    )

    model = PratmoModel.with_defaults()
    result = model.run_diurn_no2_constrained(constrained)

    out = result.output
    times_hhmm = np.array([step.time_hhmm for step in out.time_series[0].steps])
    no2 = out.species_grid("no2")
    bro = out.species_grid("bro")
    o3 = out.species_grid("o3")

    print("Final NOy scale factors:")
    for i, scale in enumerate(result.noy_scale):
        print(
            f"  box={i:02d} z={altitude_km[i]:5.1f} km "
            f"pressure={pressure_mb[i]:6.2f} mb scale={scale:.6g}"
        )

    print(f"\nModeled NO2 nearest {target_hhmm:04d} LST:")
    for i, value in enumerate(result.modeled_no2_cm3):
        print(f"  box={i:02d} modeled={value:.6e} observed={observed_no2_cm3[i]:.6e}")

    print("\nFinal species arrays:")
    print(f"  times_hhmm shape: {times_hhmm.shape}")
    print(f"  NO2 shape:        {no2.shape}")
    print(f"  BrO shape:        {bro.shape}")
    print(f"  O3 shape:         {o3.shape}")
    print("\nFirst-box NO2 time series:")
    for hhmm, value in zip(times_hhmm, no2[0]):
        print(f"  {hhmm:04d}  {value:.6e}")


if __name__ == "__main__":
    main()
