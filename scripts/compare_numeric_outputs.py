#!/usr/bin/env python3
"""Compare fixed-format PRATMO output files at numeric-token precision.

The Fortran and Rust writers intentionally differ in some headers and spacing,
so this compares numeric fields line-by-line rather than requiring byte
identity.  It reports the largest absolute and relative differences.
"""

from __future__ import annotations

import argparse
import math
import re
import sys
from pathlib import Path
from typing import Optional, Tuple


NUMBER = re.compile(
    r"(?<![A-Za-z])[-+]?(?:(?:\d+\.?(?:\d*)?)|(?:\.\d+))(?:[EeDd][-+]?\d+)?"
)


def numbers(line: str) -> list[float]:
    return [float(token.replace("D", "E").replace("d", "e")) for token in NUMBER.findall(line)]


def compare(reference: Path, candidate: Path, atol: float, rtol: float, min_reference: float) -> tuple[bool, str]:
    ref_lines = reference.read_text(errors="replace").splitlines()
    got_lines = candidate.read_text(errors="replace").splitlines()

    # boxout.dat contains implementation-specific headers before the fixed
    # separator.  The rows after that separator are the comparable numerical
    # record stream.
    if reference.name == "boxout.dat" and candidate.name == "boxout.dat":
        def data_rows(lines: list[str]) -> list[str]:
            for index, line in enumerate(lines):
                if line.count("-") > 50:
                    return lines[index + 1 :]
            raise ValueError("boxout separator not found")

        ref_lines = data_rows(ref_lines)
        got_lines = data_rows(got_lines)
    if len(ref_lines) != len(got_lines):
        return False, f"line count differs: {len(reference)}={len(ref_lines)} {len(candidate)}={len(got_lines)}"

    max_abs = 0.0
    max_rel = 0.0
    max_where = ""
    numeric_fields = 0
    first_mismatch: Optional[Tuple[int, int, float, float, float, float, str, str]] = None
    for line_no, (ref_line, got_line) in enumerate(zip(ref_lines, got_lines), 1):
        ref_values = numbers(ref_line)
        got_values = numbers(got_line)
        if len(ref_values) != len(got_values):
            return False, (
                f"numeric field count differs at line {line_no}: "
                f"{len(ref_values)} != {len(got_values)}\n"
                f"  ref: {ref_line}\n  got: {got_line}"
            )
        for field_no, (ref_value, got_value) in enumerate(zip(ref_values, got_values), 1):
            numeric_fields += 1
            diff = abs(got_value - ref_value)
            scale = max(abs(ref_value), abs(got_value), 1.0e-300)
            rel = diff / scale
            if diff > max_abs or rel > max_rel:
                max_abs = max(max_abs, diff)
                max_rel = max(max_rel, rel)
                max_where = f"line {line_no}, field {field_no}"
            if (
                (abs(ref_value) >= min_reference or abs(got_value) >= min_reference)
                and not math.isclose(got_value, ref_value, rel_tol=rtol, abs_tol=atol)
                and first_mismatch is None
            ):
                first_mismatch = (
                    line_no,
                    field_no,
                    got_value,
                    ref_value,
                    diff,
                    rel,
                    ref_line,
                    got_line,
                )
    summary = f"{numeric_fields} numeric fields; max abs={max_abs:.3e}, max rel={max_rel:.3e} ({max_where})"
    if first_mismatch is not None:
        line_no, field_no, got_value, ref_value, diff, rel, ref_line, got_line = first_mismatch
        return False, (
            f"numeric mismatch at {line_no}:{field_no}: got={got_value:.17e} "
            f"ref={ref_value:.17e} abs={diff:.3e} rel={rel:.3e}\n"
            f"  ref: {ref_line[:240]}\n  got: {got_line[:240]}\n{summary}"
        )
    return True, summary


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("reference", type=Path)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("--atol", type=float, default=5.0e-13)
    parser.add_argument("--rtol", type=float, default=5.0e-4)
    parser.add_argument("--min-reference", type=float, default=1.0e-12)
    args = parser.parse_args()
    ok, message = compare(args.reference, args.candidate, args.atol, args.rtol, args.min_reference)
    print(f"{args.reference.name}: {message}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
