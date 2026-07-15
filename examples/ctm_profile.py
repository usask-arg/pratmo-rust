"""Run a CTM altitude profile and inspect NumPy outputs."""

from __future__ import annotations

from pratmo import (
    IMPLICIT_SPECIES_NAMES,
    CtmBoxSpec,
    CtmConfig,
    PratmoModel,
)


def main() -> None:
    model = PratmoModel.with_defaults()
    output = model.run_ctm(
        CtmConfig(
            latitude_deg=60.0,
            julian_day=75,
            integration_days=10,
            boxes=[CtmBoxSpec(altitude_level=level) for level in (10, 15, 20, 25)],
            bromine=True,
            iodine=True,
        )
    )

    ozone = output.species_profile("o3")
    hydroxyl = output.species_profile("oh")
    noy = output.long_lived_profile("noy")

    print(" altitude      pressure             O3             OH        NOy VMR")
    print("      (km)          (mb)         (cm-3)         (cm-3)               ")
    for altitude, pressure, o3, oh, noy_vmr in zip(
        output.altitude_km,
        output.pressure_mb,
        ozone,
        hydroxyl,
        noy,
    ):
        print(
            f"{altitude:10.2f}  {pressure:12.4g}  {o3:13.5e}  "
            f"{oh:13.5e}  {noy_vmr:13.5e}"
        )

    print(f"\n{len(IMPLICIT_SPECIES_NAMES)} implicit species are available.")
    print(f"Solver diagnostics: {output.diagnostics!r}")


if __name__ == "__main__":
    main()
