# Python API reference

Start with the high-level interface. Results and plotting are separated from
the native compatibility layer so ordinary workflows do not need to navigate
the extension's full configuration surface.

```{toctree}
:maxdepth: 1

api-high-level
api-results
api-low-level
```

Field names accepted by result accessors are listed in
{doc}`species-and-units` and exposed programmatically as
`IMPLICIT_SPECIES_NAMES`, `LONG_LIVED_NAMES`, and `JVALUE_NAMES`.
