# Validation, provenance, and citation

```{warning}
Implementation agreement is not scientific validation. The Rust rewrite and
the iodine extension require independent scientific review before their output
is used for conclusions or decisions.
```

## What has been compared

The repository includes automated Rust and Python tests plus a clean-room
comparison with a freshly compiled Fortran reference. In opt-in parity mode,
the documented 60°N/day-75 CTM fixture and the DIURN/TPATH output sequence
match the reference at their printed precision for chemically meaningful
fields. Normal mode deliberately repairs several observable legacy quirks.

These checks establish implementation consistency for selected fixtures. They
do not validate the underlying reaction mechanism, climatologies, cross
sections, or behavior over the complete latitude/season grid.

The detailed evidence and unsupported legacy paths are maintained in
[STATUS.md](https://github.com/usask-arg/pratmo-rust/blob/main/STATUS.md) and
[FORTRAN_PARITY.md](https://github.com/usask-arg/pratmo-rust/blob/main/FORTRAN_PARITY.md).

## Known boundaries

- Full 71-latitude by 24-season CTM coverage has not been compared.
- `DERIVS`, `PZSTD`, and the original multi-case driver are unsupported.
- The iodine mechanism has mechanism-level tests but no Fortran reference.
- The legacy climatologies and photochemical data require scientific provenance
  review for any proposed application.

## Reproducible reporting

Record at least:

- the `pratmo` version and source revision;
- normal or `fortran-parity` policy;
- latitude, date/day of year, selected levels, and integration length;
- every atmospheric profile and unit;
- chemistry, photolysis, aerosol, and compatibility options;
- solver diagnostics; and
- the exact species accessor and output units used.

## Citation

The project does not yet advertise a versioned archival DOI. Until one is
available, cite the exact software release or Git revision and identify it as
the experimental Rust rewrite of PRATMO v6.0. Cite the scientific source for
any added mechanism separately; the iodine references are listed in
{doc}`IODINE_CHEMISTRY`.
