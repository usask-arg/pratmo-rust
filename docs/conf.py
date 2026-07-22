import os

# Autodoc and notebook execution use the installed wheel or editable package.
# Run `uv run maturin develop` before a local documentation build.

# Ask Plotly for portable HTML output that MyST-NB can embed. Notebook kernels
# inherit this setting from Sphinx, so examples can simply return a Figure.
os.environ.setdefault("PLOTLY_RENDERER", "notebook_connected")

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
nb_execution_mode = "force"
nb_execution_timeout = 600
nb_execution_raise_on_error = True
nb_output_stderr = "remove"

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
html_title = "PRATMO Python guide"
html_static_path = ["_static"]
html_logo = "_static/pratmo-mark.svg"
html_favicon = "_static/pratmo-mark.svg"
html_theme_options = {
    "source_repository": "https://github.com/usask-arg/pratmo-rust/",
    "source_branch": "main",
    "source_directory": "docs/",
}
