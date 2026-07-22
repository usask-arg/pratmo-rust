"""Run a multi-box DIURN cycle and inspect its shared time coordinate."""

from __future__ import annotations

import numpy as np

from pratmo import Box, DiurnalOptions, Model


def format_hhmm(value: int) -> str:
    """Format an integer local-time label without treating it as elapsed time."""
    return f"{value // 100:02d}:{value % 100:02d}"


def main() -> None:
    model = Model()
    output = model.diurnal(
        latitude=0.0,
        day="2026-04-30",
        boxes=[Box.at_level(level) for level in (10, 15, 20)],
        options=DiurnalOptions(integration_days=5, parallel_boxes=True),
    )

    no2 = output.species_grid("no2")
    oh = output.species_grid("oh")
    elapsed_hours = output.elapsed_seconds / 3600.0

    print(f"NO2 grid shape: {no2.shape} (box, time)")
    print(" altitude   local time at peak NO2        peak NO2         peak OH")
    print("      (km)                                     (cm-3)          (cm-3)")
    for box_index, altitude in enumerate(output.altitude_km):
        peak_index = int(np.argmax(no2[box_index]))
        print(
            f"{altitude:10.2f}  "
            f"{format_hhmm(int(output.time_hhmm[peak_index])):>23}  "
            f"{no2[box_index, peak_index]:14.5e}  {oh[box_index].max():14.5e}"
        )

    assert elapsed_hours[0] == 0.0
    assert elapsed_hours[-1] == 24.0
    print("\nThe orbit is ordered continuously from noon to the following noon.")


if __name__ == "__main__":
    main()
