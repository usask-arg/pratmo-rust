# OSIRIS orbit 37375 convergence regression

The two `*-pratmo-inputs.json` files are verbatim copies of baseline
nitrogen/bromine inputs exported by the OSIRIS investigation in
`arg-processing/.pytest_cache/pratmo-orbit-37375-handoff`. They contain the
81-level atmosphere and all 23 boxes (5–60 km), with explicit units and options.
No retrieval packages, NetCDF files, mounted data, or OSIRIS dependencies are
needed. The observation date is 2007-12-31; the failing scan is 37375021 and
37375022 is the neighboring control.
The supplied `reproduce.py` is retained with optional printing of the new
warning diagnostics; it also runs with 0.3.0.

`o37375022-v0.3.0-species.json` contains all 40 species profiles from the
successful full-column run before the fix, using the 0.3.0 source at main
commit `e768b55`. Units are number densities in cm^-3. The regression tolerance
is 5e-5 relative, with a 1e-20 absolute floor for inactive species.

## Reproduce

After building/installing the local PRATMO package:

```sh
python tests/fixtures/osiris/reproduce.py
python tests/fixtures/osiris/reproduce.py --all-boxes
python tests/fixtures/osiris/reproduce.py --scan 37375022
python tests/fixtures/osiris/reproduce.py --scan 37375022 --all-boxes
python -m pytest tests/test_osiris_convergence.py -q
```

## Failure diagnosis

Version 0.3.0 fails only at 17.5 km (box index 5). Its seven Newton coordinates
are HBr, CH3OOH, N2O5, HNO3, HNO4, H2O2 and HCl. At zero-based iteration 7,
HCl has density 1.075948267e8 cm^-3, mean P-L +107.898435 cm^-3 s^-1,
and Newton correction +3.978392961e8 cm^-3. Subtracting the correction gives
-2.902444694e8 cm^-3. The early component bounds have just expired. The first
unbounded HCl step was already negative at iteration 1; clipping subsequent
components did not bring the coupled system toward its root.

The finite-difference Jacobian changed XNOLD but not the named densities used
by SETUPR/HETPROB. Consequently it omitted aerosol uptake's dependence on HCl,
HBr and ClONO2. A proposed correction changed those rates only on the following
iteration, so the Jacobian and accepted chemistry described different systems.

## Fix and checks

Version 0.3.1 evaluates every perturbed or trial state with both density
representations synchronized. Daily evaluations start with the same
deterministic guesses: warm orbit guesses can select different adaptive
integration paths, making small residual comparisons inconsistent. A single
step fraction preserves the Newton direction; FIXRAT then enforces the
coupled NOy/Cly/Bry/Iy targets. Backtracking accepts only finite, nonnegative,
family-conserving states that reduce the scaled daily P-L norm. Convergence
requires both an undamped relative correction and fractional daily residual
below DAYERR (3e-5), so damping alone cannot trigger success.

The bromine-off HBr equation has identically zero P-L. Its coordinate is held
fixed in the slow solve, rather than allowing a zero pivot to generate NaNs.
Active singular, inaccurate or nonfinite solves are rejected explicitly.
Two iodine shape regressions were updated from their old, falsely converged
bromine-off orbit: HOI peaks one sample before the OH maximum, and the peak
IxOy/Iy fraction is about 2.007e-4. The common iodine fixture now also asserts
both solvers converged. The realistic-Bry abundance and conservation checks
retain their original tolerances and reference values.

Measured on macOS/arm64 with the release build:

| Run | Newton iterations (maximum) | Final relative correction (maximum) |
| --- | ---: | ---: |
| 37375021, 17.5 km | 5 | 2.20e-7 |
| 37375021, all 23 boxes | 5 | 2.87e-5 |
| 37375022, all 23 boxes | 6 | 2.81e-5 |

All runs have zero NEWRAF and RAFDAY nonconvergence counts. Serial/parallel
and 3/20/40/100-day regressions check positivity, independent noon-to-noon
closure and NOy/Cly/Bry conservation at every time sample and in snapshots.
The maximum change from the neighboring scan's released species profiles is
1.82e-5 relative (0.00182%). Core tests deliberately overshoot a chemistry
Newton direction by 100x to verify damping, residual reduction, conservation,
and that rejected trials leave the saved state unchanged.

## Recoverable failures

If a subsequent numerical solve cannot progress, RAFDAY restores the last
successfully integrated, validated cycle. It increments the nonconvergence
count and records the box, altitude and failure reason in `rafday_warnings`.
The Python API emits one aggregated `RuntimeWarning`; the CLI prints warnings.
Iteration exhaustion follows the same warning/output contract. A retained
cycle is finite usable output, but is explicitly **not established equilibrium**.
If no valid cycle has been computed, the error still propagates.

Regression tests force both iteration exhaustion and an invalid Jacobian
perturbation in serial and parallel full-column runs. All 23 boxes remain
available, family totals remain conserved and warnings identify each affected
box. The exact captured scan must still converge without exercising recovery.

The opt-in `fortran-parity` build retains the legacy correction algorithm;
the compiled Fortran differential passes for CTM and DIURN/TPATH outputs.
