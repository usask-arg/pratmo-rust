# Default settings

Defaults are chosen to make the first run executable and conservative about
experimental mechanisms. They are not universal recommendations for every
scientific problem.

## Workflow defaults

| Setting | CTM | DIURN | Reason |
|---|---:|---:|---|
| Latitude | 45°N | 45°N | Mid-latitude example |
| Day of year | 80 | 80 | Near March equinox |
| Chemistry boxes | Levels 10, 15, 20, 25 | Level 15 | Stratospheric starting selection |
| Integration length | 40 days | 20 days | Useful first equilibrium attempt |
| Atmosphere | Embedded climatology | Embedded atmosphere | Runnable without external files |

## Chemistry defaults

| Option | Default | Meaning |
|---|---:|---|
| `bromine` | `True` | Established bromine mechanism enabled |
| `iodine` | `False` | Experimental iodine mechanism requires opt-in |
| `heterogeneous` | `True` | Sulfate/sea-salt heterogeneous reactions enabled |

## Photolysis defaults

| Option | Default | Meaning |
|---|---:|---|
| `solar_flux_scale` | `1.0` | No user multiplier |
| `surface_albedo` | `0.20` | Lambertian lower boundary |
| `aerosol_extinction` | `None` | Automatically use a nonzero custom aerosol profile |

## DIURN numerical defaults

| Option | Default | Meaning |
|---|---:|---|
| `parallel_boxes` | `None` | Parallel when more than one box is selected |
| `cpp_compatibility` | `False` | Use the normal Rust/legacy time policy |
| `elapsed_time_hours` | `None` | Generate the 34-point noon-to-noon grid |

## Guardrails

- Fewer than three integration days emits `PratmoWarning`.
- A CTM latitude that is not on the internal 2.5° grid emits `PratmoWarning`
  and reports the grid point used.
- Iodine emits `ExperimentalFeatureWarning`.
- Solar scaling outside 0.5–1.5 and albedo above 0.9 emit warnings.
- Impossible units, negative concentrations, inconsistent vertical ordering,
  and mismatched array lengths raise exceptions.

See {doc}`options` for runnable alternatives and {doc}`science-scope` for what
these checks cannot establish.
