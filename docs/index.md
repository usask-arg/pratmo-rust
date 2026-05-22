# PRATMO Python API

PRATMO is a stratospheric photochemical box model. The `pratmo` Python package
provides PyO3 bindings to the Rust core, giving direct access to the CTM
(climatological transport) and DIURN (diurnal cycle) modes.

```{toctree}
:maxdepth: 2
:caption: Contents

quickstart
diurn
explorer
api
```

## Installation

Build the extension from source using [maturin](https://www.maturin.rs/):

```bash
uv run maturin develop --release
```

Then import as usual:

```python
import pratmo
```

* {ref}`genindex`
* {ref}`modindex`
