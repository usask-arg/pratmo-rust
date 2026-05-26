import os
import sys

# The installed package must be importable for autodoc and notebook execution.
# Run `uv run maturin develop` before building docs.
sys.path.insert(0, os.path.abspath("../src"))

project = "pratmo"
author = "PRATMO Authors"
release = "0.1.0"

extensions = [
    "myst_nb",
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.intersphinx",
    "sphinx_autodoc_typehints",
    "numpydoc",
]

# MyST-NB / MyST-Parser settings
myst_enable_extensions = ["colon_fence", "dollarmath"]
nb_execution_mode = "cache"
nb_execution_timeout = 600
nb_execution_raise_on_error = True

exclude_patterns = [
    "_build",
    "**.ipynb_checkpoints",
    ".jupyter_cache",
    "jupyter_execute",
]

source_suffix = {
    ".md": "myst-nb",
    ".ipynb": "myst-nb",
    ".rst": "restructuredtext",
}

autodoc_default_options = {
    "members": True,
    "undoc-members": False,
    "show-inheritance": True,
}
autodoc_typehints = "description"
autodoc_typehints_format = "short"

napoleon_numpy_docstring = True
napoleon_google_docstring = False

numpydoc_show_class_members = False

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable", None),
}

html_theme = "furo"
