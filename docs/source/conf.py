# Configuration file for the Sphinx documentation builder.

# -- Project information

project = "Nupstup"
copyright = "2022, Gurobi Optimization"
author = "Gurobi Optimization"

release = "0.1"
version = "0.1.0"

# -- General configuration

extensions = [
    "sphinx.ext.duration",
    "sphinx.ext.doctest",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.intersphinx",
    "sphinx_tabs.tabs",
    "sphinx_toolbox.collapse",
    "sphinx_toolbox.code",
    "sphinx_copybutton",
]

pygments_style = "vs"

intersphinx_mapping = {
    "python": ("https://docs.python.org/3/", None),
    "sphinx": ("https://www.sphinx-doc.org/en/master/", None),
}
intersphinx_disabled_domains = ["std"]

templates_path = ["_templates"]

# -- Include only prompts

copybutton_prompt_text = ">>> "

# -- Make examples importable for use in doctests

doctest_global_setup = """
import sys
sys.path.append("")
"""

doctest_global_cleanup = """
sys.path.pop()
"""

# -- Options for EPUB output
epub_show_urls = "footnote"
