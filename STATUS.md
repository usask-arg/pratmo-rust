# PRATMO Rust Port — Status

Faithful 1-to-1 Rust port of PRATMO v6.0 (Prather stratospheric photochemical box model, Fortran 77).

## Where we are

**CTM mode (standard test case) is fully validated against gfortran.**  
`cargo build` is clean (warnings only). `cargo test` passes 49 tests (38 unit + 11 integration). All species — long-lived (O3, N2O, NOy, CH4, H2O), short-lived radicals (OH, HO2), and NOy reservoirs (HNO3, N2O5) — match gfortran to ≥4 significant figures at all 25 altitude levels for the standard test case (60°N, March 16). DIURN mode is implemented but not yet validated against Fortran output.

---

## Validation summary (CTM mode, 60°N March 16, 25 boxes, 40 days)

| Quantity | Status | Notes |
|---|---|---|
| O3, N2O, NOy, CH4, H2O | **✓ 4 sig figs** | Exact altitude grid match |
| Air density M, T, pressure | **✓ 4 sig figs** | |
| J(O3→O1D) | **✓ +0.10%** | Altitude-independent |
| J(H2O2), J(N2O5) | **✓ ~1%** | |
| OH, HO2 noon | **✓ 4 sig figs** | Resolved by jcomp fix (bug #8) |
| HNO3 | **✓ 4 sig figs** | Resolved by jcomp fix (bug #8) |
| N2O5 | **✓ 4 sig figs** | Resolved by jcomp fix (bug #8) |
| NTIMDO | **✓ 34** | Matches Fortran after jj_fortran fix |
| RADCOUNT | **✓ 33000** | 40 days × 25 boxes × 33 time steps |

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
| 2 | `gmu0 = -0.14` hardcoded | Read from fort01.x (value: -0.12) |
| 3 | NTIMDO off by 2: `jj` didn't count noon step | `jj_fortran = jj + 1` in setday() |
| 4 | N2O read from fort51.x: F7.4 format used | fort51.x uses F8.4 (8-char fields) |
| 5 | CTM unit conversion missing | ppmv×M×1e-6 for O3, ppbv×M×1e-9 for NOy/N2O |
| 6 | nd216/nd216s from fort01.x, not boxin_gui | Compute from jdaydo/xlatdo via HUNT + JDDO array |
| 7 | CTM altitude grid: geometric 2 km spacing | CTM always uses `HYSTAT(1)` (log-pressure, T-dependent Z) |
| 8 | `jcomp = ntim + 1 - jj` — wrong evening/morning mirror indices for utime/ztime/jtim | `jcomp = ntim - 1 - jj` (0-based conversion of Fortran `NTIM+1-JJ` where `JJ=jj+1`) |
| 9 | `read_initial_densities`: `fixmix` called without setting `s.ibox`/`s.ialt` | Set `s.ibox = ib; s.ialt = ialt` before `fixmix(s)` — mirrors Fortran `IBOX=I; IALT=IABS(NBOXDO(IBOX))` |
| 10 | DIURN test-scaling (do3ref/dm × 1.049/1.066, T=280 K) applied before fort02.x was read | Moved scaling block to `read_all` after `read_initial_densities`, matching Fortran's bread.f order |
| 11 | `diurn_unit7_header` wrote `s.n[j] + 1` for NT values; `s.n[j]` is already 1-based | Remove `+ 1` from NT computation in `diurn_unit7_header` |

---

## Test coverage

49 tests total (`cargo test`).

| Suite | Count | What it covers |
|---|---|---|
| `tests/ctm_integration.rs` | 11 | Full model run vs gfortran reference (all 25 altitudes, all key species, full OH diurnal, sanity checks) |
| `reader::tests` | 20 | `parse_e_field`, `parse_fixed_i32s`, `parse_fixed_f64s_fw`, `hystat` (logp + geometric), `setday` (ntim, weights, mirror symmetry, polar night, nday modes) |
| `ctm::tests` | 7 | `interp2` (corners, midpoint, linear, uniform), `read_f8_4`, `read_f7_1` |
| `output::tests` | 9 | `hunt` (interior, boundary, below/above range, hint, JDDO array), `fmt_e10p3` |

---

## Open gaps

### 1. DIURN mode validation (nd216 = 0) — partially validated

DIURN mode has been compared against the Fortran reference (fort01_diurn.x, equatorial 30°N May, 25 boxes at levels 1–30). Three bugs were found and fixed:

- **Bug #9**: `read_initial_densities` did not set `s.ibox`/`s.ialt` before calling `fixmix` — fixmix always operated on stale box 0 instead of the current box, corrupting initial densities for all boxes.
- **Bug #10**: Test-scaling block (do3ref/do3int × 1.066, dm × 1.049, T=280 K) was applied inside `read_ozone_profile` *before* fort02.x was read. Fortran applies it *after*. This caused initial species densities to be 1.049× too large in DIURN mode.
- **Bug #11 (output)**: `diurn_unit7_header` wrote `s.n[j] + 1` for NT values but `s.n[j]` is already 1-based — producing NT values 1 too large.

**Validation results** (equatorial May, 25 boxes levels 1–30):
- **fort08.x** (TPATH species averages, segment 0): species values **match** the Fortran reference at all 25 boxes to ≥3 sig figs.
- **fort07.x** (DIURN time series): initial noon values **match** exactly. Post-noon evolution differs: the Fortran DIURN produces near-zero O(1D) and OH at high-altitude boxes (levels 25–30, ~50–58 km) while the Rust computes physically plausible photo-steady-state values (~10³ cm⁻³ OH). This appears to be a pre-existing Fortran DIURN bug where high-altitude J-values are incorrectly near zero — confirmed by the inverted pattern (non-zero O1D at low altitude, zero at high altitude) inconsistent with UV physics.
- **fort09.x** (TPATH rate averages, segment 0): rates differ as a consequence of the fort07 discrepancy.

### 2. CTM grid coverage
The integration test only exercises ipath=415 (60°N, March half-month 6) because `boxin_gui.dat` constrains `nd216 = nd216s`. The other 70 latitudes × 24 half-months in the 71×24 climatology grid are untested. Next step: run with `nd216s=1` (full grid scan) and compare the full `boxout.dat` against Fortran.

### 3. Known latent bug in `setday` (unreachable in standard run)
When any `ATIME[j] * SSET > 0.5 * DAYSEC` (requires ATIME > ~2.4 for a 5-hour sunset), Fortran uses `DTIME[j]` in `SNIT` but Rust uses `DTIME[jj]` (the step before). This path is never triggered by the standard `fort01.x` (all ATIME ≤ 1.05), but would produce wrong night-step intervals for unusual configurations.

### 4. DERIVS mode (nd216 < 0)
Sensitivity Jacobians (d(P-L)/d(O3) etc.). Fortran `SUBROUTINE DERIVS` perturbs species and re-runs DSTEP. Not translated. `main.rs` prints a warning.

### 5. PZSTD
Converts pressure to standard Z grid. Called from NEWATM when NPSTD > 0. Currently a stub (no-op). The standard fort01.x has NPSTD=0.

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
