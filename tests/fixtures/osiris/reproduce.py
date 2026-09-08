"""Replay captured OSIRIS inputs using only PRATMO and its dependencies.

Examples:
    python reproduce.py                         # failing 17.5 km box
    python reproduce.py --all-boxes             # failing full baseline
    python reproduce.py --scan 37375022         # successful neighboring box
    python reproduce.py --scan 37375022 --all-boxes

Exceptions intentionally propagate and give a nonzero exit status.
"""

import argparse
import json
from pathlib import Path

import pratmo


def load_arguments(scan: int) -> dict:
    source = Path(__file__).parent / f"o{scan}-pratmo-inputs.json"
    saved = json.loads(source.read_text())
    return {
        "latitude": saved["latitude"],
        "day": saved["day"],
        "atmosphere": pratmo.Atmosphere(**saved["atmosphere"]),
        "boxes": [pratmo.Box(**box) for box in saved["boxes"]],
        "initial_mixing_ratios": [
            pratmo.LongLivedMixingRatios(**mix)
            for mix in saved["initial_mixing_ratios"]
        ],
        "chemistry": pratmo.ChemistryOptions(**saved["chemistry"]),
        "photolysis": pratmo.PhotolysisOptions(**saved["photolysis"]),
        "options": pratmo.DiurnalOptions(**saved["options"]),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scan", type=int, choices=(37375021, 37375022), default=37375021)
    parser.add_argument("--all-boxes", action="store_true")
    parser.add_argument("--box-index", type=int, default=5)
    args = parser.parse_args()
    kwargs = load_arguments(args.scan)
    if not args.all_boxes:
        index = args.box_index
        if not 0 <= index < len(kwargs["boxes"]):
            parser.error("box-index must be in [0, 22]")
        kwargs["boxes"] = [kwargs["boxes"][index]]
        kwargs["initial_mixing_ratios"] = [kwargs["initial_mixing_ratios"][index]]
    print(f"PRATMO module: {pratmo.__file__}", flush=True)
    print(f"Scan: {args.scan}; boxes: {len(kwargs['boxes'])}", flush=True)
    result = pratmo.Model().diurnal(**kwargs)
    print("Integration returned successfully.")
    diagnostics = result.diagnostics
    for name in (
        "newraf_nonconvergence_count",
        "rafday_nonconvergence_count",
        "rafday_max_final_relative_correction",
        "rafday_max_correction_iterations",
        "rafday_warnings",
    ):
        print(f"{name}: {getattr(diagnostics, name, [])}")


if __name__ == "__main__":
    main()
