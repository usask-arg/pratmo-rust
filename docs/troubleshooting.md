# Troubleshooting

## The import fails

Confirm that the wheel supports the active Python and platform, then run:

```bash
python -c "import pratmo; from importlib.metadata import version; print(version('pratmo'))"
```

Install `pratmo[plot]` when a plotting helper reports that Plotly is missing,
and `pratmo[io]` for xarray/NetCDF examples.

## A standard level is at an unexpected altitude

Standard-level pressure is fixed, but altitude is derived from the seasonal
temperature profile. Always use `output.altitude_km`; the table in
{doc}`standard-levels` is only an orientation.

## The custom atmosphere is rejected

Check that all arrays are finite and the same length, pressure decreases with
index, altitude increases, and units are named explicitly. Ozone values around
`1–10` are probably ppmv, whereas a dimensionless fraction is around `1e-6`.

## Exact altitude selection fails

`Box.at_altitude(...)` requires a custom `Atmosphere` with explicit altitude.
The requested height must lie inside that grid.

## Aerosol has no radiative effect

Per-box aerosol controls heterogeneous chemistry. Radiative extinction needs a
profile-level `Atmosphere.aerosol_surface_area` and
`PhotolysisOptions(aerosol_extinction=True)` or the automatic `None` default.

## A short run warns or results change with integration length

Increase `integration_days`, inspect `output.diagnostics`, and compare the
scientifically important fields between successively longer runs. A solver
success alone is not an equilibrium test.

## DIURN local times appear out of order

The integration starts at noon and ends at noon the following day. Plot and
sort with `elapsed_seconds`, not the cyclic `time_hhmm` labels.

## The J-value grid does not vary with time

The current structured DIURN output exposes daily-mean J-values repeated across
the time dimension. It is shape-compatible with species grids but is not a
time-resolved actinic-flux product.

## The iodine run warns

This is intentional. Iodine is experimental and opt-in. Read
{doc}`IODINE_CHEMISTRY` and {doc}`validation` before using it.

## Performance

Runtime grows approximately with integration days, number of boxes, and
mechanism complexity. Start with one box while developing a workflow, then use
the automatic multi-box parallel default. Avoid shortening production runs
only to reduce runtime without rechecking convergence.
