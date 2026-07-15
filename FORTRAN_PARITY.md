# Fortran parity mode

The `fortran-parity` Cargo feature reproduces known quirks in the original
PRATMO v6.0 executable without making them the default behavior of the Rust
model.

```bash
cargo run --release -p pratmo-cli --features fortran-parity -- --input-dir fortran
```

The feature currently enables these legacy behaviors:

- replace the `GMU0` value read from `fort01.x` with `-0.14`;
- replace the atmospheric temperature profile with 280 K in the READIN
  calibration block (CTM later loads its climatological temperature);
- leave `XJDO` at zero during each `NEWRAF` solve, updating H2O photolysis only
  for the post-solve diagnostic `CHEMS` call;
- use the original seasonal `fort04.x` aerosol climatology rather than the
  corrected/default aerosol handling;
- reproduce DIURN's zero-initialized `SSF` spectral multipliers (all DIURN
  photolysis is therefore zero in the reference executable);
- reproduce DIURN's zero heterogeneous rate constants (170–177);
- omit R177 from HNO3 production in the DIURN chemistry balance, as the
  reference `CHEMPL` does;
- preserve the pre-loop `PUNCH(0,0)` metadata altitude and the final `LEND`
  mixing-ratio dump.

The feature is forwarded by `pratmo-cli` and `pratmo-py`. Normal builds keep
the corrected behavior.

## Repairs enabled in every build

These are translation or input-compatibility fixes, not intentional legacy
quirks, so they are always active:

- accept both the original 40-species `fort10_cam06.x` and the extended
  iodine dataset;
- honor `aero_sf`, `ivarO3`, `irunclim`, and `iTfull` controls;
- keep CTM output open across multiple climatology cases and emit the same
  39-level altitude range as Fortran;
- retain `N2BOX` semantics when a negative `NBOX` requests propagation of the
  first `fort02.x` column;
- preserve and apply the result of `FIXRAT` in RAFDAY;
- correct two DSTEP one-based index translations; and
- repair NEWATM multi-line profile reads, old-density rescaling, and per-box
  `FIXMIX` context;
- parse the fixed-format ODF ISR field before its trailing band annotation; and
- preserve the Fortran `SETDAY` evening-time mirror and average-mode WTIME
  indices;
- select the pressure-level `HYSTAT(1)` grid for `P:*` atmospheres; and
- accept embedded blanks in malformed fixed-width `fort02.x` E-fields, as
  gfortran's formatted reader does.

## Source-level policy map

Every conditional is documented next to its implementation.  This map gives
the quick navigation point and the normal-mode counterpart:

| Source | Fortran-parity behavior | Normal Rust behavior |
|---|---|---|
| `reader.rs` | hard-code `GMU0=-0.14`; flatten the calibration temperature to 280 K | honor `fort01.x` GMU0; retain the supplied temperature profile |
| `jvalue.rs` | pass `PIRAY`/`PIAER` in the legacy atmospheric-level order | pass them in optical-depth order |
| `ctm.rs` | reload seasonal aerosol optical depths from `fort04.x` | retain the corrected/current aerosol profile |
| `chemistry.rs` | zero DIURN rates 170–177 and omit R177 from HNO3 production | keep heterogeneous rates and the complete R177 balance |
| `diurnal.rs` | zero uninitialized `SSF`; emit PUNCH metadata before the box loop; repeat the legacy RAFDAY relaxation after every slow-species correction | keep `SSF=1`; emit metadata after selecting the first box; relax once so HNO3/BrONO2 RAFDAY corrections are retained |
| `path.rs` | solve NEWRAF with stale-zero `XJDO`; emit the final LEND dump | include H2O photolysis in the solve; omit the legacy dump |

The fixed `HYSTAT(1)` and embedded-blank input handling are unconditional
Fortran-format compatibility repairs, not physical-policy switches.

## Automated parity-policy coverage

`pratmo-core/tests/fortran_parity_feature.rs` is compiled in both configurations
so the feature cannot silently drift from normal Rust behavior. It checks the
original 40-species spectral input and ODF annotations, the GMU0/280-K setup
policy, DIURN's `SSF` initialization, the zero heterogeneous-rate block, the
R177 HNO3 balance, a public one-box DIURN photolysis smoke run, and a full
25-box CTM smoke run in parity mode. It also exercises representative tropical,
mid-latitude, and polar-night DIURN configurations for finite-state safety. The
normal configuration asserts the corrected/input-driven counterpart for each
policy where applicable. These are policy and execution guards; the clean-room
gfortran differential below is the authoritative numerical regression.

Run the focused guards with:

```bash
cargo test -p pratmo-core --test fortran_parity_feature -- --test-threads=1
cargo test -p pratmo-core --features fortran-parity \
  --test fortran_parity_feature -- --test-threads=1
```

For a clean-room differential covering both CTM and DIURN/TPATH, run:

```bash
scripts/fortran_differential.sh
```

The harness copies and normalizes the legacy inputs into a temporary workspace,
compiles the reference with gfortran, redirects the original Windows output
path locally, runs Rust with `fortran-parity`, and compares CTM `boxout.dat` plus
DIURN/TPATH `fort07`/`fort08`/`fort09` numeric records. Repository inputs and
generated outputs are not modified. Values below `1e-12` use the absolute
floating-point floor; all larger fields use a `5e-4` relative tolerance. Pass
`ctm` or `diurn` to run one mode; `normal` is a CTM diagnostic for the
intentionally corrected Rust behavior, while the default parity invocation is
the release gate.

## Differential result (60°N, Julian day 75)

The reference was freshly compiled with gfortran 14.2 using 8-byte default
reals and run for the standard 25-box, 40-day case. The original CRLF inputs
were normalized in a temporary copy; the repository inputs were not changed.

| Quantity | Rust `fortran-parity` vs gfortran |
|---|---:|
| altitude, T, pressure, air | exact at output precision |
| O3, N2O, CH4, H2O, NOy, Cly, Bry | exact at output precision |
| J(BrNO3), J(NO2) | exact to printed precision (maximum relative difference <0.01%) |
| NO, NO2, NO3, HNO3, N2O5, OH, HO2 and other radicals | exact to printed precision for fields ≥1e-12 (maximum relative difference <0.01%) |
| solver RADCOUNT | 33,066 (same as Fortran) |
| solver RAXLOOP | 121,141 vs Fortran 121,162 (same converged solution) |

### DIURN / TPATH mode (ND216 = 0)

The standard 25-box DIURN case now matches the compiled reference through
`fort07.x`, `fort08.x`, `fort09.x`, and the final `LEND` snapshot when built
with `--features fortran-parity`.  The only observed difference is a single
last-digit formatter roundoff in one twilight radical field.  Without the
feature, Rust keeps `SSF=1` and heterogeneous chemistry active, which is the
physically useful behavior and intentionally differs from the legacy
executable.

The only conspicuous relative errors are twilight halogen fields at the
floating-point floor (below `5e-15` in the output volume-mixing-ratio units);
the largest absolute difference is below `2e-18`. All chemically meaningful
fields match the reference output at the formatter's precision. A one-day
differential run also matches the full 25-box diurnal state to printed
precision.
