# pratmo

> [!WARNING]
> `pratmo` is an experimental, AI-assisted rewrite of the PRATMO v6.0
> stratospheric photochemical box model. It has **not been scientifically
> validated**, including the additional iodine chemistry. Do not use its
> output for scientific conclusions, operational decisions, or safety-critical
> work without independent review and validation.

`pratmo` rewrites the original Fortran 77 model in Rust and exposes it through
a command-line runner and a typed Python/NumPy interface. The project aims to
preserve the original model structure while making the code easier to test,
inspect, and embed. It also adds an experimental inorganic iodine mechanism
that is not present in PRATMO v6.0.

Implementation tests and selected numerical comparisons with the compiled
Fortran program are included, but these are software consistency checks—not
scientific validation of the Rust rewrite, its inputs, or its predictions.

## Features

- Rust model core with structured CTM steady-state and full DIURN-cycle APIs.
- Python package named `pratmo`, with typed configuration objects and NumPy
  arrays for profiles, grids, coordinates, time series, and diagnostics.
- 40 implicit gas-phase species and 19 long-lived species or chemical families.
- Experimental inorganic iodine chemistry, including ten iodine species,
  iodine photolysis, higher iodine oxides, and heterogeneous recycling.
- Standard atmospheres, custom atmospheric profiles, configurable aerosol and
  sea-salt surface areas, and observation-constrained NO2 workflows.
- Compatibility paths for the original fixed-format inputs and legacy CLI
  outputs, plus opt-in tooling for numerical comparisons with Fortran.
- Runnable CTM, DIURN, iodine, custom-atmosphere, and OMPS batch examples.

The iodine mechanism's assumptions and known omissions are documented in
[Iodine chemistry](docs/IODINE_CHEMISTRY.md).

## Python quickstart

The published distribution and import name are both `pratmo`. To install a
released wheel:

```bash
python -m pip install pratmo
```

For development, build the extension into the project environment:

```bash
uv sync --dev
uv run maturin develop --release
```

Run a small CTM altitude profile:

```python
from pratmo import CtmBoxSpec, CtmConfig, PratmoModel

model = PratmoModel.with_defaults()
output = model.run_ctm(
    CtmConfig(
        latitude_deg=60.0,
        julian_day=75,
        integration_days=40,
        boxes=[CtmBoxSpec(altitude_level=level) for level in (10, 15, 20, 25)],
    )
)

print(output.altitude_km)
print(output.species_profile("o3"))       # cm^-3
print(output.long_lived_profile("noy"))  # dimensionless mixing ratio
```

The package exposes `IMPLICIT_SPECIES_NAMES`, `LONG_LIVED_NAMES`, and
`JVALUE_NAMES` for programmatic discovery. Field lookup is case-insensitive.
See the [Python documentation](docs/index.md) and [examples](examples/) for CTM,
DIURN, custom-atmosphere, and observation-constrained workflows.

## Rust and CLI

```bash
cargo build --workspace
cargo run -p pratmo-cli -- --input-dir fortran
```

`pratmo-core` provides the structured Rust API. `pratmo-cli` preserves the
original file-based workflow and writes Fortran-compatible output files into
the selected input directory.

## Development verification

These checks test software behavior and selected Fortran agreement. Passing
them does not mean that the model is scientifically validated.

Run the normal development checks with:

```bash
cargo fmt --all -- --check
cargo clippy --workspace --all-targets --all-features -- -D warnings
cargo test --workspace
cargo test --workspace --all-features
uv run pytest
uv run sphinx-build -W --keep-going -b html docs docs/_build/html
```

GitHub Actions runs these Rust checks, the Python tests, a strict documentation
build, and the complete cross-platform wheel build on every push and pull
request. Built distributions are retained as workflow artifacts. The publish
job runs only for release tags matching `v*` and uses PyPI Trusted Publishing.

The clean-room compiled-Fortran differential is available separately:

```bash
scripts/fortran_differential.sh
```

The detailed numerical comparisons, parity policy, and unsupported legacy modes
are tracked in [STATUS.md](STATUS.md) and [FORTRAN_PARITY.md](FORTRAN_PARITY.md).

## Current boundaries

- `DERIVS` sensitivity mode and `PZSTD` conversion are rejected explicitly.
- The original multi-case Fortran driver loop is not implemented.
- Full CTM climatology-grid coverage remains broader than the comparison
  fixture.
- The iodine extension has mechanism-level tests but no Fortran reference and
  has not been scientifically validated.

## License

`pratmo` is available under the [MIT License](LICENSE).
