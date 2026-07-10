# PRATMO Rust Port — Agent Guide

## Goal

Rewrite PRATMO v6.0 (Prather stratospheric photochemical box model) from Fortran 77 to Rust as a faithful 1-to-1 port. Preserve subroutine boundaries, original fixed-format input file formats, and numerical results. The Rust binary must produce output files (`boxout.dat`, `fort07.x`, etc.) that match the Fortran reference to ≥4 significant figures on long-lived species and ~1% on J-values.

The model computes stratospheric photochemistry for a set of altitude boxes over N days, tracking 30 implicit species (NO, OH, O3, etc.) and 18 long-lived families (NOy, N2O, CH4, …) under a full diurnal cycle of UV photolysis.

---

## Repository layout

```
pratmo-core/src/   — Rust library (all physics)
  chemistry.rs     — SETUPR + CHEMS (221 reactions, rate constants)
  jvalue.rs        — SOL → JVALUE → OPTAU → SCATTR → SPHERE
  solver.rs        — NEWRAF, FIXMIX, LINSLV
  reader.rs        — FortranReader: fort10/11/01/02/13/14, SETDAY, HYSTAT
  diurnal.rs       — DIURN, DAILY, RAFDAY
  path.rs          — TPATH, DSTEP, NEWATM
  ctm.rs           — CTMLFQ, CTOUTP, write_box_row (boxout.dat writer)
  output.rs        — PUNCH, PRTALL, PRTPTH, HUNT
  init.rs          — CTINIT
  heterogeneous.rs — HETPROB (H2SO4 aerosol heterogeneous chemistry)
  tracers.rs       — Bry/Cly/CH4/N2O empirical relationships
  state.rs         — ModelState struct (all 24 COMMON blocks unified)
  constants.rs     — NB, NL, NQ, NJVAL, etc.
pratmo-cli/src/main.rs — CLI entry point
fortran/           — original Fortran 77 source (reference)
STATUS.md          — detailed module map, validation results, known gaps
```

---

## Building and running

```bash
cargo build                                          # build everything
cargo run -p pratmo-cli -- --input-dir fortran       # run with standard test case
```

The standard test case is CTM mode: `fortran/fort01.x` with ND216=44 (44 lat/month combinations, overridden to a single atmosphere by boxin_gui.dat).

Key outputs written to the `--input-dir`:
- `boxout.dat`  — main CTM output table (825 cols × ~58 rows)
- `fort07.x`    — PUNCH diurnal species time series (DIURN mode)
- `fort08.x`    — PRTPTH species averages (TPATH mode)
- `fort09.x`    — PRTPTH rate averages (TPATH mode)

---

## Compiling the Fortran reference

The original Fortran sources are in `fortran/`. The executables bundled (`pratmo.exe`, `batmo.x`) are Windows PE32+ binaries and cannot run on macOS/Linux.

### Prerequisites

- gfortran (macOS: `brew install gcc`)
- Xcode Command Line Tools SDK

### Input file preprocessing (required once)

The input files have Windows `\r\r\n` line endings that confuse gfortran. Convert them:

```bash
cd fortran
python3 -c "
import os
for f in ['fort03_LLM.x','fort04.x','fort05.x','fort51.x',
          'fort10_cam06.x','fort11_jpl09.x','fort01.x',
          'fort13.x','fort14.x','J_H2O_SZA0.dat','boxin_gui.dat']:
    if os.path.exists(f):
        data = open(f,'rb').read().replace(b'\r', b'')
        open(f,'wb').write(data)
"
```

### Fortran source patch

`batmo.f` uses `data cinpdir/'.\'` (Windows path) in a COMMON block, which gfortran rejects as a non-standard DATA statement. Patch it:

```bash
python3 -c "
with open('batmo.f') as f: src = f.read()
# Add assignment before first OPEN statement
src = src.replace('      RAXLOOP = 0', \"      cinpdir = './'\n      RAXLOOP = 0\")
open('batmo_unix.f', 'w').write(src)
"
```

### Compile

```bash
gfortran -O2 -fdefault-real-8 -fdefault-double-8 \
  -ffixed-form -ffixed-line-length-none \
  -fdec -fno-backslash -w \
  -isysroot /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk \
  -o pratmo_gf \
  batmo_unix.f bchem.f bctin.f bctmx.f bdiel.f bdiff.f bhetp.f \
  bjval.f bpath.f bread.f bry_vs_tracers.gen.vers03.f butil.f HUNT.FOR
```

Note: `ClNO3code.f` is excluded because `HETPROB` is already defined in `bhetp.f`.

### Run

```bash
./pratmo_gf
```

The clean-room numerical differential (requires gfortran) is:

```bash
scripts/fortran_differential.sh
```

It compares both the standard CTM `boxout.dat` and DIURN/TPATH `fort07`–`fort09`
records in temporary normalized input trees. Use `ctm` or `diurn` to select one
mode.

The binary reads from the current directory. Outputs go to `boxout.dat` (path from boxin_gui.dat) and to `G:\software\...` (literal Windows path, created as a file on macOS).

---

## Key Fortran→Rust translation notes

### Known tricky mappings

| Fortran pattern | Rust equivalent |
|---|---|
| `COMMON/CCDEN/ DNO(NB)...` | `s.dno[ib]`, accessed via `den_get/den_set` |
| `DDDDDD(NB,NSPEC) ≡ DNO(NB)` | `s.den_get(ib, spec_idx)` (0-based) |
| `NT(J) ≡ N1,N2,...,N30,NTOT,NTOTX` | `s.n[j]` (1-based slot values) |
| `VVVVVV(NL,NJVAL) ≡ VNO(NL)` | `s.jval_get(box, jval_idx)` |
| `RCOLUM(430) ≡ XR(30)++R(250)++...` | `rcolum_get(s, j)` in diurnal.rs |
| `SSF(NQ)` in COMMON/CHRIS/ | `s.ssf[k]` — **must be initialized to 1.0** |
| `HYSTAT(1)` in CTM mode | always call `hystat(s, 1)` in ctmlfq |
| `JCOMP = NTIM+1-JJ` (1-based) | `jcomp = ntim - 1 - jj` (0-based jj): converting JJ=jj+1 gives NTIM+1-(jj+1)-1 = ntim-jj-1 |
| fort51.x format | F8.4 (8-char fields), NOT F7.4 |
| CTM climatology units | DO3INP: ppmv × DM × 1e-6 → cm⁻³; DNOyINP/DN2OINP: ppbv × DM × 1e-9 |

### Index conventions

- `s.n[j]` stores 1-based XN slot indices (n[0]=9 means NO is at XN slot 9).
- XN/XNOFT/XXNOFT are indexed 0-based in Rust; use `s.n[species] - 1` to get 0-based slot.
- `jval_get(il, k)` uses il = box index (0-based), k = J-value index (0-based: 0=jNO, 3=jO3d, 8=jH2O2, ...).

---

## Validation state (2026-05-20)

### CTM mode — fully validated (60°N, March 16, 25 boxes, 40 days)

| Quantity | Agreement | Notes |
|---|---|---|
| O3, N2O, NOy, CH4, H2O | 4 sig figs | |
| J(O3→O1D) | +0.10% | altitude-independent |
| J(H2O2) | +1.21% | altitude-independent |
| OH, HO2 noon | 4 sig figs | resolved by jcomp fix |
| HNO3 | 4 sig figs | resolved by jcomp fix |
| N2O5 | 4 sig figs | resolved by jcomp fix |

### DIURN mode — validated in parity mode (equatorial May, 25 boxes levels 1–30)

| Output | Agreement | Notes |
|---|---|---|
| fort08.x species averages | ≥3 sig figs | All 25 boxes match Fortran reference |
| fort07.x noon initial values | exact | Initial conditions from fort02.x correct |
| fort07.x post-noon time series | diverges | Fortran shows O1D≈0 at high altitude — suspected Fortran J-value bug; Rust values are physically correct |
| fort09.x rate averages | differs | Consequence of fort07 divergence |

---

## Open tasks (priority order)

1. **DIURN integration tests** — the parity feature now has policy guards, a representative latitude/season invariant matrix, and an end-to-end CTM smoke test; the compiled gfortran differential is documented in `FORTRAN_PARITY.md`.

2. **DERIVS mode** (`nd216 < 0`) — sensitivity Jacobians; not implemented yet.

3. **CTM grid coverage** — only 60°N March tested; full 71×24 grid untested.

4. **PZSTD** — pressure-to-standard-Z grid conversion; stub; only needed if NPSTD > 0.

5. **Multi-case loop** — `batmo.f` loops over multiple READIN cases; T1DIFF is a no-op so this is low priority.
