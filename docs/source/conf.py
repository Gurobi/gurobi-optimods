# Configuration file for the Sphinx documentation builder.

import gurobi_optimods

# -- Project information

project = "Gurobi OptiMods"
copyright = "2022, Gurobi Optimization"
author = "Gurobi Optimization"

version = gurobi_optimods.__version__
release = version

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

html_theme = "sphinx_rtd_theme"

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
