"""Plot the inorganic iodine diurnal cycle near 24 km."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from pratmo import DiurnBoxSpec, DiurnConfig, PratmoModel


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("iodine_24km.png"),
        help="Output PNG path (default: %(default)s)",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    model = PratmoModel.with_defaults()
    output = model.run_diurn(
        DiurnConfig(
            latitude_deg=0.0,
            julian_day=120,
            integration_days=20,
            boxes=[DiurnBoxSpec(altitude_level=13)],
            iodine=True,
        )
    )

    air_density = output.air_density_cm3[0]
    elapsed_hours = output.elapsed_seconds / 3600.0
    ppt_per_cm3 = 1.0e12 / air_density
    species = {
        "I": output.species_grid("i")[0] * ppt_per_cm3,
        "IO": output.species_grid("io")[0] * ppt_per_cm3,
        "HOI": output.species_grid("hoi")[0] * ppt_per_cm3,
        "IONO2": output.species_grid("iono2")[0] * ppt_per_cm3,
        "HI": output.species_grid("hi")[0] * ppt_per_cm3,
    }

    figure, axis = plt.subplots(figsize=(10, 5))
    for name, values in species.items():
        axis.plot(elapsed_hours, values, linewidth=1.8, label=name)

    ticks = list(range(0, 25, 2))
    axis.set_xticks(ticks, [f"{(12 + hour) % 24:02d}:00" for hour in ticks])
    axis.set_xlim(0, 24)
    axis.set_xlabel("Local solar time (noon to noon)")
    axis.set_ylabel("Mixing ratio (ppt)")
    axis.set_title(
        f"Iodine diurnal cycle at {output.altitude_km[0]:.0f} km "
        "(0°N, Julian day 120)"
    )
    axis.legend()
    axis.grid(True, alpha=0.3)
    figure.tight_layout()

    args.output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.output, dpi=150)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
