# pratmo documentation

```{warning}
`pratmo` is an experimental, AI-assisted Rust rewrite of the PRATMO v6.0
stratospheric photochemical box model. It has not been scientifically
validated. The iodine chemistry is an additional experimental mechanism and
has not been validated either. Do not treat model output as suitable for
scientific, operational, or safety-critical use without independent review and
validation.
```

The `pratmo` Python package provides PyO3 bindings to the Rust core, giving
direct access to CTM steady-state and DIURN diurnal-cycle workflows. Selected
comparisons with the Fortran executable are implementation checks only; they
do not validate the model, its inputs, or its scientific predictions.

## Feature overview

- Rust core, legacy-compatible command-line runner, and typed Python/NumPy API.
- 40 implicit gas-phase species and 19 long-lived species or families.
- Standard and custom atmospheric profiles with altitude, pressure,
  temperature, and density coordinates.
- CTM profiles, full diurnal time series, J-values, diagnostics, and
  observation-constrained NO2 runs.
- Configurable aerosol and sea-salt surface areas.
- Experimental inorganic iodine chemistry with ten gas-phase iodine species,
  photolysis, higher oxides, and heterogeneous recycling.
- Runnable profile, diurnal-cycle, iodine, custom-atmosphere, and batch examples.

```{toctree}
:maxdepth: 2
:caption: Contents

quickstart
diurn
custom-atmosphere
explorer
IODINE_CHEMISTRY
iodine_saiz_lopez_2014_upgrade
api
releasing
```

## Installation

Create the development environment and build the extension from source:

```bash
uv sync --dev
uv run maturin develop --release
```

Then import as usual:

```python
import pratmo
```

The published distribution and import are both named `pratmo`. The core package
requires NumPy. Install `pratmo[io]` for the xarray/NetCDF
batch example and `pratmo[plot]` for plotting examples.

* {ref}`genindex`
* {ref}`modindex`
