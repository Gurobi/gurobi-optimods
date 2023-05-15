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
    "sphinx_copybutton",
    "sphinx_tabs.tabs",
    "sphinx_toolbox.code",
    "sphinx_toolbox.collapse",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.extlinks",
    "sphinx.ext.doctest",
    "sphinx.ext.duration",
    "sphinx.ext.intersphinx",
    "sphinxcontrib.bibtex",
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

extlinks_detect_hardcoded_links = True
extlinks = {
    "pypi": ("https://pypi.org/project/%s/", "%s"),
}

# -- Bibfiles

bibtex_bibfiles = [
    "refs/graphs.bib",
    "refs/qubo.bib",
    "refs/regression.bib",
]
