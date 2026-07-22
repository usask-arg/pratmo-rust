# Python examples

Install the package first:

```bash
python -m pip install "pratmo[plot,io]"
```

## Starting points

| Example | Purpose | External data |
|---|---|---|
| `ctm_profile.py` | Four-level climatological profile and diagnostics | None |
| `diurnal_cycle.py` | Multi-box noon-to-noon chemistry | None |
| `iodine_diurnal_cycle.py` | Experimental iodine time series and PNG | None |
| `no2_constrained_custom_atmosphere.py` | Small custom profile with observed-NO2 scaling | None |

Run one from the repository root, for example:

```bash
python examples/ctm_profile.py
```

## Data-driven workflows

`omps_no2_constrained_batch.py` and `osiris_pratmo_comparison.py` are advanced
instrument workflows. They require caller-supplied NetCDF files; inspect each
script's `--help` output for the exact paths and variables. The OSIRIS workflow
also samples the original PRATMO climatology tables from a legacy data
directory.

```bash
python examples/omps_no2_constrained_batch.py --help
python examples/osiris_pratmo_comparison.py --help
```

Every example uses the high-level `Model`, `Atmosphere`, and `Box` interface
unless a native compatibility object is specifically required by an archived
workflow. Experimental iodine examples emit a warning by design.
