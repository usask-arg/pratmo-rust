# PRATMO Python guide

```{warning}
`pratmo` is an experimental, AI-assisted Rust rewrite of the PRATMO v6.0
stratospheric photochemical box model. It has not been scientifically
validated. The iodine extension is experimental as well. Treat the outputs as
research software results that require independent scientific review.
```

PRATMO calculates gas-phase, photolytic, and heterogeneous chemistry in a set
of independent atmospheric boxes. The recommended Python interface is
{class}`pratmo.Model`, with {class}`pratmo.Atmosphere` for custom profiles and
{class}`pratmo.Box` for selecting chemistry levels.

## Install

Install the model with the interactive plotting helpers used in this guide:

```bash
python -m pip install "pratmo[plot]"
```

Then run a climatological profile:

```python
from pratmo import Model

model = Model()
profile = model.ctm()
print(profile.altitude_km)
print(profile.species_profile("o3"))
```

No external input directory is needed for the standard workflows. Developers
building from source should follow the repository README; those commands are
intentionally kept out of the user installation path.

## Choose a path

- Start with {doc}`quickstart` for a guided CTM profile and DIURN cycle.
- Read {doc}`concepts` before interpreting model output.
- Use {doc}`atmospheric-inputs` for pressure, temperature, ozone, altitude,
  aerosol, and long-lived-family inputs.
- Use {doc}`options` for chemistry, radiation, and numerical experiments.
- Check {doc}`science-scope` and {doc}`validation` before quantitative work.

```{toctree}
:maxdepth: 1
:caption: Getting started

quickstart
concepts
standard-levels
species-and-units
```

```{toctree}
:maxdepth: 1
:caption: Workflows

ctm
diurn
atmospheric-inputs
custom-atmosphere
options
results
explorer
```

```{toctree}
:maxdepth: 1
:caption: Science and validation

science-scope
validation
IODINE_CHEMISTRY
```

```{toctree}
:maxdepth: 1
:caption: Reference

defaults
troubleshooting
api
```
