#!/usr/bin/env python3
"""Compare the final RAXLOOP/RADCOUNT diagnostics from two runs."""

from __future__ import annotations

import re
import sys
from pathlib import Path


def fortran_metrics(text: str) -> tuple[int, int]:
    matches = re.findall(r"(?m)^\s*(\d+(?:\.\d+)?)\s+(\d+(?:\.\d+)?)\s*$", text)
    if not matches:
        raise ValueError("Fortran diagnostics line not found")
    raxloop, radcount = matches[-1]
    return round(float(raxloop)), round(float(radcount))


def rust_metrics(text: str) -> tuple[int, int]:
    match = re.search(r"RAXLOOP=([0-9]+)\s+RADCOUNT=([0-9]+)", text)
    if match is None:
        raise ValueError("Rust diagnostics line not found")
    return int(match.group(1)), int(match.group(2))


def main() -> int:
    reference, candidate = map(Path, sys.argv[1:3])
    expected = fortran_metrics(reference.read_text(errors="replace"))
    got = rust_metrics(candidate.read_text(errors="replace"))
    if expected[1] != got[1]:
        print(f"RADCOUNT mismatch: Rust={got[1]} Fortran={expected[1]}")
        return 1
    print(f"diagnostics: RADCOUNT={got[1]}, RAXLOOP Rust={got[0]} Fortran={expected[0]}")
    return 0


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit("usage: compare_diagnostics.py fortran.stdout rust.stdout")
    raise SystemExit(main())
