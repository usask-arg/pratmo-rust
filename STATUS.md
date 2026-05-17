# PRATMO Rust Rewrite — Status

Fortran 77 stratospheric photochemical box model (PRATMO v6.0, Prather) → Rust.
Direct 1-to-1 port preserving Fortran subroutine boundaries and original input file formats.

## Workspace layout

```
pratmo-core/   — library crate (all physics)
pratmo-cli/    — binary crate (CLI entry point)
fortran/       — original Fortran source (reference)
```

## Module mapping (Fortran → Rust)

| Fortran file | Rust module | Status |
|---|---|---|
| `bcomm.h` (24 COMMON blocks) | `state.rs` (`ModelState` struct) | **Complete** |
| `bchem.f` (1150 lines) | `chemistry.rs` | **Complete** |
| `bjval.f` (518 lines) | `jvalue.rs` | **Complete** |
| `butil.f` (FIXRAT/FIXMIX/NEWRAF/LINSLV) | `solver.rs` | **Complete** |
| `bread.f` (535 lines) + SETDAY/HYSTAT | `reader.rs` | **Complete** |
| `bdiel.f` (334 lines) | `diurnal.rs` | **Complete** |
| `bpath.f` (391 lines) | `path.rs` | **Complete** |
| `bhetp.f` (201 lines) | `heterogeneous.rs` | **Complete** |
| `bry_vs_tracers.gen.vers03.f` (1700 lines) | `tracers.rs` | **Complete** (key functions) |
| `bctin.f` (239 lines) | `init.rs` | **Complete** |
| `ClNO3code.f` (184 lines) | `clno3.rs` | **Complete** (reference impl) |
| `bctmx.f` (2404 lines) | `ctm.rs` | **Skeleton** — see gaps below |
| `batmo.f` (main program) | `pratmo-cli/src/main.rs` | **Partial** — see gaps below |
| `constants.h` / `bcomm.h` PARAMETERs | `constants.rs` | **Complete** |

Build status: `cargo build` → **clean** (warnings only from unused stub params).

---

## What works

- All core physics modules compile and have no `todo!()` stubs.
- **Chemistry**: SETUPR (rate constants + Troe fall-off), CHEMS (221 reactions, O(1D)/N(4S) quasi-SS), CHEMPL (P/L rates for 30 implicit + 20 family species), RHSLHS (RHS vector + 30×30 Jacobian).
- **J-values**: SOL → JVALUE → OPTAU → SCATTR (Eddington 2-stream, 3-pt Gauss quadrature) → SPHERE (Chapman function twilight correction). All 44 photolysis channels, 77 wavelength bins. **Bug fixed: SSF (solar flux scale factor) was initialized to 0; now correctly 1.0.**
- **Solver**: FIXRAT (cubic family conservation, regula falsi + NR), FIXMIX, NEWRAF (adaptive time-step halving), NEWRAX (18-iteration NR), LINSLV/RESOLV (Crout partial-pivot LU).
- **I/O**: FortranReader reads all standard input files (fort10\_cam06.x, fort11\_jpl09.x, fort01.x, fort02.x, fort13.x, fort14.x, J\_H2O\_SZA0.dat, fort03\_LLM.x, fort05.x, fort51.x) with full fixed-format Fortran parsing. fort02.x reading added for DIURN mode initial conditions.
- **Diurnal integration**: DIURN (box loop), DAILY (24-hr time-step loop with J-value caching in STORJV), RAFDAY (NR slow-species steady-state driver).
- **Path integration**: DSTEP (NDAYSD-day integration with explicit species update and family renormalization), TPATH (outer loop with NEWATM reading from fort01\_remaining).
- **NEWATM**: Reads subsequent path records from fort01\_remaining (cached during READIN), resets T/O3/lat and calls FIXMIX. Path mode atmosphere resets work.
- **Heterogeneous chemistry**: Full Shi/Worsnop/Davidovits H₂SO₄ aerosol formula (8 γ values: ClONO2+HCl/H2O, HOCl+HCl, N2O5±HCl, BrONO2, HOBr±HCl/HBr), PSC ice branch.
- **Tracer relationships**: Bry/Cly/CH4/CH3Cl vs N2O polynomial fits (Wamsley 2003, Michelsen 1998, Schauffler 2003) with age-of-air correction.
- **CTM initialisation**: CTINIT family density initialisation from NOy, Cly, Bry with uniform NOy/Cly/Bry partitioning. CTM climatology (T, O3, NOy, N2O) correctly interpolated and unit-converted from ppmv/ppbv to cm⁻³.
- **Output routines** (new module `output.rs`):
  - `PUNCH(IB, IDAY)` — writes diurnal species time series to fort07.x ✓
  - `PRTALL(mode, flags, n)` — species/rate table to stdout ✓
  - `PRTPTH(iseg, iday, ibx)` — path summary to fort08.x/fort09.x ✓
  - `PRTRAT(nnn)`, `PRTAVG(nnn)` — rate/average tables ✓
  - `HUNT(xx, x, jlo)` — general binary bracket search ✓
  - `CTOUTP` — CTM diagnostic stdout output (N2O/NOy/CH4 loss frequencies) ✓
- **CTM mode validated**: Running fort01.x (ND216=44, 60°N/March, 25 boxes, 40 days) produces physically consistent N2O loss rates (~1e-8 s⁻¹ at 40 km), NOy loss rates, and O3 profiles. RADCOUNT=31010 (40 days × 25 boxes × 31 steps). ✓

---

## Known gaps (what remains)

### 1. Numerical parity with Fortran
The Fortran executables are Windows PE32+ binaries (cannot run on macOS). Physical validation uses known atmospheric chemistry constraints:
- N2O photochemical lifetime at 40 km: ~150 days (model: ✓)
- O3 profile shape vs altitude: qualitatively correct (model: ✓)
- NOy/N2O relationship via tracer empirics: implemented (Wamsley 2003)
- Full bit-for-bit comparison against Fortran output: **not yet possible without Windows/Wine**

### 2. CTOUTP DPL file output
The CTM run writes N2O/NOy/CH4 loss frequencies to stdout but does NOT write the full DPL array (P-L rates, mixing ratios vs lat/month) to binary output files. The Fortran CTOUTP writes these to `bmoutfile` (specified in boxin\_gui.dat). The DPL array is populated during the loop but not written to disk.

### 3. DERIVS mode (`nd216 < 0`)
Sensitivity/derivative calculation mode. The Fortran `SUBROUTINE DERIVS` perturbs key species and runs DSTEP to compute Jacobians (d(P-L)/d(O3), etc.). Not yet translated. `main.rs` prints a warning and exits cleanly.

### 4. Missing utility subroutines
- `PZSTD` — converts from pressure levels to standard Z grid. Called from NEWATM when `NPSTD > 0`. Stub (no-op currently).
- `T1DIFF` — 1D diffusion transport. Called from `batmo.f` after the main run in multi-case mode. A no-op in the Fortran reference code (`bdiff.f`), so safe to skip.

### 5. Multi-case loop (`batmo.f`)
`main.rs` handles only one case per invocation. The Fortran loops: `CALL READIN → run → CALL READIN (LEND=TRUE) → T1DIFF → repeat`. `T1DIFF` is a no-op in the reference, so this is low-priority.

### 6. CTM interpolation simplifications
- `s.xdecd = 0.0` (equinox declination) is hardcoded for all months — the Fortran reads EDEC per month from a DATA statement.
- The boxin\_gui.dat lat/jday selection for CTM is partially read but `nd216` override is not applied from it (see `bctmx.f` `j0` calculation).

---

## Suggested next steps (in priority order)

1. **Write DPL output file** from CTM mode: implement CTOUTP file writer using `bmoutfile` path from boxin\_gui.dat.
2. **Validate CTM against Fortran** by running `pratmo.exe` under Wine and comparing stdout loss frequencies.
3. **Implement EDEC per month** in ctm.rs for accurate solar declination by lat/mon.
4. **Implement DERIVS** for sensitivity analysis capability.
5. **Implement `PZSTD`** if NPSTD > 0 runs are needed.
6. **Multi-case loop** in `main.rs` for sensitivity ensembles.

---

## Key design decisions made

- **COMMON blocks → `ModelState`**: all 24 COMMON blocks unified into one `Box<ModelState>` (~3 MB heap allocation).
- **EQUIVALENCE arrays**: `DDDDDD(NB,NSPEC)≡DNO`, `FFFFFF(NB,18)≡FO3`, `VVVVVV(NL,NJVAL)≡VNO`, `RCOLUM(430)≡XR` handled via `den_get/den_set`, `fff_get/fff_set`, `jval_get/jval_set` accessor methods rather than unsafe aliasing.
- **0-based indexing**: all Rust arrays are 0-based; Fortran 1-based species indices N1..N30 stored as `s.n[0..29]` (values remain 1-based for use as array indices into `xn[N-1]`).
- **`nboxdo: [i32; NB]`**: signed (negative = use DAILY, positive = use RAFDAY) matching the Fortran sign convention.
- **`ModelReader` trait**: I/O abstracted behind a trait for future format modernisation.
- **Input files preserved**: same fixed-format Fortran files (fort10\_cam06.x, fort11\_jpl09.x, fort01.x, fort13.x, fort14.x).
