"""Run a CTM altitude profile and inspect NumPy outputs."""

from __future__ import annotations

from pratmo import (
    IMPLICIT_SPECIES_NAMES,
    Box,
    ChemistryOptions,
    CtmOptions,
    Model,
    mixing_ratio_as,
)


def main() -> None:
    model = Model()
    output = model.ctm(
        latitude=60.0,
        day="2026-03-16",
        boxes=[Box.at_level(level) for level in (10, 15, 20, 25)],
        chemistry=ChemistryOptions(),
        options=CtmOptions(integration_days=40),
    )

    ozone = output.species_profile("o3")
    hydroxyl = output.species_profile("oh")
    noy = mixing_ratio_as(output.long_lived_profile("noy"), "ppbv")

    print(" altitude      pressure             O3             OH      NOy (ppbv)")
    print("      (km)         (hPa)         (cm-3)         (cm-3)                ")
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
