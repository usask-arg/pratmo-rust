# Releasing `pratmo`

The Python distribution and import name are both `pratmo`. The current package
version is `0.1.0`.

Distributions are built by `.github/workflows/release.yml` with maturin on every
push and pull request. The workflow builds wheels for Linux, musllinux, Windows,
and macOS, plus a source distribution, and retains them as workflow artifacts.
For tags beginning with `v`, the same successfully built files are uploaded to
PyPI using Trusted Publishing.

The extension targets PyO3's `abi3-py39` stable ABI. Each platform therefore
produces one `cp39-abi3` wheel that can be installed on Python 3.9 and newer,
instead of rebuilding the extension for every Python minor version.

Every pushed commit and pull request also runs `.github/workflows/ci.yml`. CI
checks Rust formatting and Clippy, runs both the default and all-feature Rust
test suites, builds and installs a representative Linux wheel, runs the Python
tests, and builds the documentation with warnings treated as errors. The CI
workflow has read-only repository permissions and contains no publishing step.

## One-time PyPI setup

Configure a PyPI Trusted Publisher with these values:

- Owner: `usask-arg`
- Repository: `pratmo-rust`
- Workflow: `release.yml`
- Environment: `pypi`

Create a protected GitHub environment named `pypi` as well. No long-lived PyPI
API token is needed.

## Publishing a release

Before tagging, keep the version synchronized in `pyproject.toml` and the root
Cargo workspace metadata, then run the complete verification commands from the
project README. Tag names use the Python version with a leading `v`:

```bash
git tag -a v0.1.0 -m "v0.1.0"
git push origin v0.1.0
```

The publish job only runs for matching tags and only after all wheel and source
distribution jobs finish successfully. Publishing packages does not change the
project's experimental, unvalidated status.
