# Configuration file for the Sphinx documentation builder.

import gurobi_optimods

# -- Project information

project = "Gurobi OptiMods"
copyright = "2023, Gurobi Optimization"
author = "Gurobi Optimization"

version = gurobi_optimods.__version__
release = version

# -- General configuration

extensions = [
    "numpydoc",
    "sphinx_copybutton",
    "sphinx_design",
    "sphinx_tabs.tabs",
    "sphinx_toolbox.code",
    "sphinx_toolbox.collapse",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosectionlabel",
    "sphinx.ext.autosummary",
    "sphinx.ext.doctest",
    "sphinx.ext.duration",
    "sphinx.ext.extlinks",
    "sphinx.ext.intersphinx",
    "sphinxcontrib.bibtex",
]

pygments_style = "vs"

intersphinx_mapping = {
    "python": ("https://docs.python.org/3/", None),
    "sphinx": ("https://www.sphinx-doc.org/en/master/", None),
    "pandas": ("https://pandas.pydata.org/pandas-docs/stable/", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "scipy": ("https://docs.scipy.org/doc/scipy/", None),
    "networkx": ("https://networkx.org/documentation/stable", None),
}

templates_path = ["_templates"]

html_theme = "sphinx_rtd_theme"

autosectionlabel_prefix_document = True

# -- Include only prompts
copybutton_prompt_text = ">>> "

# -- Options for EPUB output
epub_show_urls = "footnote"

extlinks_detect_hardcoded_links = True
extlinks = {
    "pypi": ("https://pypi.org/project/%s/", "%s"),
    "ghsrc": ("https://github.com/Gurobi/gurobi-optimods/tree/main/%s", "%s"),
}

# -- Bibfiles

bibtex_bibfiles = [
    "refs/graphs.bib",
    "refs/portfolio.bib",
    "refs/qubo.bib",
    "refs/regression.bib",
    "refs/workforce.bib",
]

# -- numpydoc magic linking

numpydoc_xref_param_type = True
numpydoc_xref_aliases = {
    "DataFrame": "pandas.DataFrame",
    "DiGraph": "networkx.DiGraph",
    "Graph": "networkx.Graph",
    "LinAlgError": "numpy.linalg.LinAlgError",
    "spmatrix": "scipy.sparse.spmatrix",
}
numpydoc_xref_ignore = {"optional", "or", "of"}


# -- Docstring preprocessing for autodoc

autodoc_typehints = "none"
autodoc_docstring_signature = True
add_module_names = False


def process_signature(app, what, name, obj, options, signature, return_annotation):
    """Replace the create_env keyword argument accepted by decorated mods with
    the parameters added by the decorator"""

    if what in ["function", "method"] and hasattr(obj, "_decorated_mod"):
        if "create_env" not in signature:
            raise ValueError(f"Decorated mod {name} does not accept create_env")
        new_signature = signature.replace(
            "create_env", "verbose=True, logfile=None, solver_params=None"
        )
        print(f"Modified signature of {name}")
        return new_signature, return_annotation

    return signature, return_annotation


boilerplate = """
    **verbose** : :ref:`bool <python:bltin-boolean-values>`, optional
        ``verbose=False`` suppresses all console output

    **logfile** : :class:`python:str`, optional
        Write all mod output to the given file path

    **solver_params** : :class:`python:dict`, optional
        Gurobi parameters to be passed to the solver
"""
boilerplate = boilerplate.split("\n")


def process_docstring(app, what, name, obj, options, lines):
    """Add parameter entries for decorated mods"""

    if what in ["function", "method"] and hasattr(obj, "_decorated_mod"):
        # Find where the last input parameter is listed
        in_paramlist = False
        lineno = None
        for i, line in enumerate(lines):
            if ":Parameters:" in line:
                in_paramlist = True
            elif in_paramlist and (
                ":Returns:" in line or "processed by numpydoc" in line
            ):
                lineno = i - 1
                break

        if lineno is None:
            raise ValueError(f"Failed to find param list for {name}")

        # Insert boilerplate bits
        for line in reversed(boilerplate):
            lines.insert(lineno, line)

    if what == "module":
        lines.append("")
        lines.append(f"The following mods can be imported from ``{name}``:")


def setup(app):
    app.connect("autodoc-process-signature", process_signature)
    app.connect("autodoc-process-docstring", process_docstring)
