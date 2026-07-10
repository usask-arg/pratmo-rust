#!/usr/bin/env bash
# Reproducible clean-room differential against the original Fortran.
#
# The default `all` mode compares both the standard CTM case and a DIURN/TPATH
# case. Inputs are normalized only in temporary directories; repository files
# and generated outputs are never modified.

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REQUEST="${1:-all}"
RUST_MODE="parity"
if [[ "$REQUEST" == normal ]]; then
    REQUEST="ctm"
    RUST_MODE="normal"
fi
case "$REQUEST" in
    all|ctm|diurn) ;;
    *) echo "usage: $0 [all|ctm|diurn|normal]" >&2; exit 2 ;;
esac
if ! command -v gfortran >/dev/null 2>&1; then
    echo "gfortran is required for the Fortran differential" >&2
    exit 2
fi

WORK="$(mktemp -d "${TMPDIR:-/tmp}/pratmo-differential.XXXXXX")"
if [[ "${KEEP_WORK:-0}" != 1 ]]; then
    trap 'rm -rf "$WORK"' EXIT
fi

INPUT_FILES=(
    fort01.x fort03_LLM.x fort04.x fort05.x fort10_cam06.x fort11_jpl09.x
    fort13.x fort14.x fort51.x J_H2O_SZA0.dat boxin_gui.dat
)
SOURCE_FILES=(
    bchem.f bctin.f bctmx.f bdiel.f bdiff.f bhetp.f bjval.f bpath.f bread.f
    bry_vs_tracers.gen.vers03.f butil.f HUNT.FOR
)

for file in "${SOURCE_FILES[@]}"; do
    tr -d '\r' < "$ROOT/fortran/$file" > "$WORK/$file"
done
tr -d '\r' < "$ROOT/fortran/batmo.f" > "$WORK/batmo.f"
tr -d '\r' < "$ROOT/fortran/bcomm.h" > "$WORK/bcomm.h"

# gfortran rejects the original Windows DATA statement for CINPDIR. Inject a
# Unix working-directory assignment into the temporary source only.
awk '{ if ($0 == "      RAXLOOP = 0") print "      cinpdir = '\''./'\''"; print }' \
    "$WORK/batmo.f" > "$WORK/batmo_unix.f"

gfortran -O2 -fdefault-real-8 -fdefault-double-8 \
    -ffixed-form -ffixed-line-length-none -fdec -fno-backslash -w \
    -I"$WORK" \
    -isysroot "$(xcrun --show-sdk-path 2>/dev/null || printf '/')" \
    -o "$WORK/pratmo_gf" \
    "$WORK/batmo_unix.f" \
    "${SOURCE_FILES[@]/#/$WORK/}"

prepare_run_tree() {
    local case_name="$1"
    local side="$2"
    local run_dir="$WORK/$case_name/$side"
    mkdir -p "$run_dir"
    for file in "${INPUT_FILES[@]}"; do
        tr -d '\r' < "$ROOT/fortran/$file" > "$run_dir/$file"
    done

    # The fixture contains the original Windows output path. Redirect the
    # temporary reference and Rust runs to a local boxout.dat.
    awk 'NR == 13 { print "./boxout.dat"; next } { print }' \
        "$run_dir/boxin_gui.dat" > "$run_dir/boxin_gui.local.dat"
    mv "$run_dir/boxin_gui.local.dat" "$run_dir/boxin_gui.dat"

    if [[ "$case_name" == diurn ]]; then
        # Fortran reads unconnected unit 2 from fort.2; Rust's reader uses the
        # explicit fort02.x name. Both receive the same checked-in fixture.
        cp "$ROOT/pratmo-core/data/fort02.x" "$run_dir/fort.2"
        cp "$ROOT/pratmo-core/data/fort02.x" "$run_dir/fort02.x"
        python3 - "$run_dir/fort01.x" <<'PY'
from pathlib import Path
import sys

path = Path(sys.argv[1])
lines = path.read_text().splitlines(True)
for index, line in enumerate(lines):
    if line.startswith("XRF/RLX/LB"):
        lines[index] = line.replace("   44   44", "    0    0")
        break
else:
    raise SystemExit("XRF/RLX/LB control line not found")
path.write_text("".join(lines))
PY
    fi
}

run_case() {
    local case_name="$1"
    local fortran_run="$WORK/$case_name/fortran"
    local rust_run="$WORK/$case_name/rust"
    prepare_run_tree "$case_name" fortran
    prepare_run_tree "$case_name" rust

    (cd "$fortran_run" && "$WORK/pratmo_gf" > fortran.stdout 2> fortran.stderr)
    if [[ "$RUST_MODE" == parity ]]; then
        cargo run --quiet --release -p pratmo-cli --features fortran-parity -- \
            --input-dir "$rust_run" > "$rust_run/rust.stdout" 2> "$rust_run/rust.stderr"
    else
        cargo run --quiet --release -p pratmo-cli -- \
            --input-dir "$rust_run" > "$rust_run/rust.stdout" 2> "$rust_run/rust.stderr"
    fi
    python3 "$ROOT/scripts/compare_diagnostics.py" \
        "$fortran_run/fortran.stdout" "$rust_run/rust.stdout"

    if [[ "$case_name" == ctm ]]; then
        python3 "$ROOT/scripts/compare_numeric_outputs.py" \
            "$fortran_run/boxout.dat" "$rust_run/boxout.dat"
    else
        python3 "$ROOT/scripts/compare_numeric_outputs.py" \
            "$fortran_run/fort.7" "$rust_run/fort07.x"
        python3 "$ROOT/scripts/compare_numeric_outputs.py" \
            "$fortran_run/fort.8" "$rust_run/fort08.x"
        python3 "$ROOT/scripts/compare_numeric_outputs.py" \
            "$fortran_run/fort.9" "$rust_run/fort09.x"
    fi
}

if [[ "$REQUEST" == all || "$REQUEST" == ctm ]]; then
    run_case ctm
fi
if [[ "$REQUEST" == all || "$REQUEST" == diurn ]]; then
    run_case diurn
fi

echo "Fortran differential passed for $REQUEST in $RUST_MODE mode."
if [[ "${KEEP_WORK:-0}" == 1 ]]; then
    echo "Kept differential workspace: $WORK"
fi
