# PRATMO Rust Port — Status

Experimental, AI-assisted Rust rewrite of PRATMO v6.0 (Prather stratospheric
photochemical box model, Fortran 77), with additional experimental iodine
chemistry. The rewrite and iodine extension have not been scientifically
validated. The comparisons below measure implementation agreement with selected
Fortran fixtures; they do not establish scientific validity.

## Where we are

**CTM and short-lived species now match the compiled Fortran reference in the
opt-in parity mode.**
The original 40-species Fortran input now runs directly, and an opt-in
`fortran-parity` feature reproduces known undesirable executable behavior.
In the standard 60°N, day-75 case, atmospheric fields, long-lived families,
photolysis rates, and all chemically meaningful radical fields match at output
precision. See [FORTRAN_PARITY.md](FORTRAN_PARITY.md) for the differential
command and measurements.

---

## Fortran comparison summary (60°N day 75, 25 boxes, 40 days)

| Quantity | Status | Notes |
|---|---|---|
| O3, N2O, NOy, CH4, H2O | **✓ 4 sig figs** | Exact altitude grid match |
| Air density M, T, pressure | **✓ 4 sig figs** | |
| output J(BrNO3), J(NO2) | **✓ <0.01%** | Exact to printed precision in `fortran-parity` build |
| OH, HO2 final noon | **✓ <0.01%** | Exact to printed precision |
| HNO3 | **✓ <0.01%** | Exact to printed precision |
| N2O5 | **✓ <0.01%** | Exact to printed precision |
| NTIMDO | **✓ 34** | Matches Fortran after jj_fortran fix |
| RADCOUNT | **✓ 33066** | Matches the compiled Fortran run |

---

## Release sign-off gate

The reproducible clean-room differential now passes in `fortran-parity` mode:

```bash
scripts/fortran_differential.sh
```

This compiles gfortran in a temporary normalized input tree, compares CTM
`boxout.dat`, DIURN/TPATH `fort07`–`fort09`, and checks exact `RADCOUNT`. The
default and parity Rust test matrices are also green. The remaining open gaps
below are unsupported or broader-than-fixture Fortran entry paths, not failures
of the CTM/DIURN implementation-comparison gate.

---

## Module map

| Fortran | Rust | Status |
|---|---|---|
| `bcomm.h` (24 COMMON blocks) | `state.rs` (`ModelState`) | ✓ Complete |
| `bchem.f` | `chemistry.rs` | ✓ Complete |
| `bjval.f` | `jvalue.rs` | ✓ Complete |
| `butil.f` | `solver.rs` + `output.rs` | ✓ Complete |
| `bread.f` + SETDAY/HYSTAT | `reader.rs` | ✓ Complete |
| `bdiel.f` | `diurnal.rs` | ✓ Complete |
| `bpath.f` | `path.rs` | ✓ Complete |
| `bhetp.f` | `heterogeneous.rs` | ✓ Complete |
| `bry_vs_tracers.gen.vers03.f` | `tracers.rs` | ✓ Complete |
| `bctin.f` | `init.rs` | ✓ Complete |
| `ClNO3code.f` | `clno3.rs` | ✓ Reference impl |
| `bctmx.f` | `ctm.rs` | ✓ CTM + boxout.dat output |
| `batmo.f` | `main.rs` | ✓ Single-case (no multi-case loop) |
| `bdiff.f` | — | T1DIFF is a no-op in Fortran; skipped |

---

## Bugs found and fixed

| # | Bug | Fix |
|---|---|---|
| 1 | `ssf[NQ] = 0.0` — silenced ALL J-values | Initialize to 1.0 in `ModelState::new()` |
| 2 | `gmu0 = -0.14` hardcoded | Normal builds read fort01.x; `fortran-parity` reproduces -0.14 |
| 3 | NTIMDO off by 2: `jj` didn't count noon step | `jj_fortran = jj + 1` in setday() |
| 4 | N2O read from fort51.x: F7.4 format used | fort51.x uses F8.4 (8-char fields) |
| 5 | CTM unit conversion missing | ppmv×M×1e-6 for O3, ppbv×M×1e-9 for NOy/N2O |
| 6 | nd216/nd216s from fort01.x, not boxin_gui | Compute from jdaydo/xlatdo via HUNT + JDDO array |
| 7 | CTM altitude grid: geometric 2 km spacing | CTM always uses `HYSTAT(1)` (log-pressure, T-dependent Z) |
| 8 | `jcomp = ntim + 1 - jj` — wrong evening/morning mirror indices for utime/ztime/jtim | `jcomp = ntim - 1 - jj` (0-based conversion of Fortran `NTIM+1-JJ` where `JJ=jj+1`) |
| 9 | `read_initial_densities`: `fixmix` called without setting `s.ibox`/`s.ialt` | Set `s.ibox = ib; s.ialt = ialt` before `fixmix(s)` — mirrors Fortran `IBOX=I; IALT=IABS(NBOXDO(IBOX))` |
| 10 | DIURN test-scaling (do3ref/dm × 1.049/1.066, T=280 K) applied before fort02.x was read | Moved scaling block to `read_all` after `read_initial_densities`, matching Fortran's bread.f order |
| 11 | `diurn_unit7_header` wrote `s.n[j] + 1` for NT values; `s.n[j]` is already 1-based | Remove `+ 1` from NT computation in `diurn_unit7_header` |
| 12 | `SETDAY` mirrored `DTIME` one slot late, changing the first evening integration weights | Use `NTIM-1-j` for the zero-based destination |
| 13 | ODF ISR values were parsed together with their trailing wavelength-band annotation | Parse the first fixed-format token only |
| 14 | Parity-mode CHEMPL omitted R177 from HNO3 in every mode, although Fortran includes it in CTM | Restrict the legacy omission to `nd216 == 0` (DIURN) |

### `bugs.txt` follow-up

- The current `boxin_gui.dat` has all 24 IPFR flags, and the Rust CTM reader
  now rejects short or non-integer IPF/IPFR/IPJV lines instead of silently
  padding them with zero.
- The tracer boundary typo from item 12 is already corrected in Rust; tests
  cover both Cly and CH3Cl out-of-range sentinel behavior.
- Normal Rust Newton solves now use a conservative f64 cancellation floor for
  `RAFPML`, covering the dark-BrCl precision case from item 19. The
  `fortran-parity` feature retains the original strict Fortran criterion.

---

## Test coverage

90 Rust tests pass in the default configuration; 5 policy-specific short-lived
assertions are ignored in the normal build. The `fortran-parity` all-features
configuration has 64 active Rust tests (55 unit tests plus 9 feature-policy
tests). The Python binding suite adds 36 passing integration tests.

| Suite | Count | What it covers |
|---|---|---|
| `tests/ctm_integration.rs` | 11 | 6 active CTM regressions + 5 policy-specific ignored comparisons |
| `tests/custom_diurn.rs` | 4 | Public custom-atmosphere DIURN API |
| `tests/fortran_parity_feature.rs` | 8 normal / 9 parity | Original 40-species input, ODF parsing, policy switches, public DIURN smoke, representative latitude/season matrix, and end-to-end parity CTM smoke coverage |
| `tests/iodine_chemistry.rs` | 13 | Iodine chemistry, photolysis, convergence, and heterogeneous recycling |
| `chemistry::iodine_tests` | 7 | Reaction stoichiometry, atom conservation, analytic Jacobian, and evaluated rate points |
| `reader::tests` | 23 | `parse_e_field`, `parse_fixed_i32s`, `parse_fixed_f64s_fw`, `hystat` (logp + geometric), `setday` (ntim, weights, time-grid mirror, polar night, nday modes) |
| `ctm::tests` | 8 | `interp2` (corners, midpoint, linear, uniform), strict boxin integer counts, `read_f8_4`, `read_f7_1` |
| `output::tests` | 9 | `hunt` (interior, boundary, below/above range, hint, JDDO array), `fmt_e10p3` |
| `solver::tests` | 4 | machine-precision convergence floor, cancellation guard, and iodine-family scaling |
| `tracers::tests` | 2 | Cly/CH3Cl fit boundaries and out-of-range sentinels |
| `api::tests` | 6 | Configuration validation, cyclic HHMM lookup, J-value mapping, and rejection of DERIVS/PZSTD |
| `tests/test_pratmo.py` | 36 | Python configuration/output objects, discoverable names, NumPy profiles/grids, and CTM/DIURN runs |

---

## Open gaps

### 1. DIURN mode comparison (nd216 = 0) — matched in parity mode

DIURN mode has been compared against the compiled Fortran reference (equatorial
30°N May, 25 boxes at levels 1–30). The parity build now matches the complete
DIURN/TPATH output sequence; the normal build intentionally retains physical
photolysis and heterogeneous chemistry.

The setup and output repairs include:

- **Bug #9**: `read_initial_densities` did not set `s.ibox`/`s.ialt` before calling `fixmix` — fixmix always operated on stale box 0 instead of the current box, corrupting initial densities for all boxes.
- **Bug #10**: Test-scaling block (do3ref/do3int × 1.066, dm × 1.049, T=280 K) was applied inside `read_ozone_profile` *before* fort02.x was read. Fortran applies it *after*. This caused initial species densities to be 1.049× too large in DIURN mode.
- **Bug #11 (output)**: `diurn_unit7_header` wrote `s.n[j] + 1` for NT values but `s.n[j]` is already 1-based — producing NT values 1 too large.
- **DIURN parity quirk**: the Fortran executable leaves `SSF` and
  heterogeneous rates 170–177 at zero. The Rust feature gate reproduces this
  observable behavior; default Rust runs keep both active.
- **Path/output parity**: restored the pre-loop `PUNCH(0,0)` metadata record,
  fixed fixed-width E-field blank handling in `fort02.x`, and emitted the
  `LEND` final mixing-ratio snapshot.

**Comparison results** (equatorial May, 25 boxes levels 1–30, parity feature):
- **fort07.x**: DIURN time series and metadata match to printed precision,
  with one last-digit twilight roundoff.
- **fort08.x / fort09.x**: all segment/day/box species and rate rows match.
- **LEND snapshot**: final mixing-ratio records match the Fortran output.

### 2. CTM grid coverage
The clean-room differential now exercises the standard 60°N/day-75 CTM case,
and the public DIURN matrix covers tropical, mid-latitude, and polar-night
atmospheres. The other 70 latitudes × 24 half-months in the 71×24 CTM
climatology grid remain untested as a full sweep. A future campaign can run
with `nd216s=1` and compare the complete `boxout.dat` grid. The invariant
matrix deliberately uses polar-night cases; polar-day edge cases that exceed
the legacy 16-angle `STORJV` storage remain outside the supported fixture.

### 3. Known latent bug in `setday` (unreachable in standard run)
When any `ATIME[j] * SSET > 0.5 * DAYSEC` (requires ATIME > ~2.4 for a 5-hour sunset), Fortran uses `DTIME[j]` in `SNIT` but Rust uses `DTIME[jj]` (the step before). This path is never triggered by the standard `fort01.x` (all ATIME ≤ 1.05), but would produce wrong night-step intervals for unusual configurations.

### 4. DERIVS mode (nd216 < 0)
Sensitivity Jacobians (d(P-L)/d(O3) etc.). Fortran `SUBROUTINE DERIVS` perturbs species and re-runs DSTEP. Not translated. The CLI and public API now reject this mode explicitly instead of running a false-success no-op.

### 5. PZSTD
Converts pressure to standard Z grid. Called from NEWATM when NPSTD > 0. Currently unsupported; the CLI, public API, and TPATH now reject `NPSTD > 0` explicitly. The standard fort01.x has NPSTD=0.

### 6. Multi-case loop
`batmo.f` can loop over multiple cases (READIN → run → READIN with LEND=TRUE → T1DIFF → repeat). T1DIFF is a no-op. `main.rs` handles only the first case.

---

## Key design decisions

- **COMMON blocks → `ModelState`**: all 24 COMMON blocks unified into one `Box<ModelState>` (~3 MB heap).
- **EQUIVALENCE arrays**: `DDDDDD`, `FFFFFF`, `VVVVVV`, `RCOLUM` handled via `den_get/den_set`, `fff_get/fff_set`, `jval_get/jval_set` rather than unsafe aliasing.
- **N1..N30 indices**: `s.n[j]` stores 1-based XN slot values. Use `s.n[j] - 1` for 0-based array access.
- **SSF must be 1.0**: `s.ssf[k] = 1.0` (solar flux scale factor) is set by CTMLFQ in the Fortran; in Rust it is initialized to 1.0 in `ModelState::new()`.
- **CTM always uses HYSTAT(1)**: log-pressure altitude grid (Z from T), not geometric 2 km spacing.
- **Input files unchanged**: same fixed-format Fortran files. Input files need CRLF→LF cleaning on macOS (see `agents.md`).
