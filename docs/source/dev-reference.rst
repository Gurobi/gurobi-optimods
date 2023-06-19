Developer Reference
===================

This page contains further information on the tools/packages we use, to guide
developers who are working on mods.

The ``@optimod`` decorator
--------------------------

The ``@optimod`` decorator should be added to any Mod functions that need to
create Gurobi models. :doc:`adding` provides a template for how to apply the
decorator. The decorator provides the following arguments automatically:

* ``verbose``, enabling users to shut off output;
* ``logfile``, enabling users to send Mod logs to a file; and
* ``solver_params``, enabling users to pass a dictionary of parameters to the
  Gurobi Optimizer.

As a Mod developer, your Mod must take a keyword argument ``create_env``. This
function should be called to create a Gurobi environment inside your Mod, and
all gurobipy ``Model`` objects should use this environment. The environment must
be properly closed before your Mod function returns (this is best achieved by
using context managers).

Adding Citations
----------------

We use
`sphinxcontrib-bibtex <https://sphinxcontrib-bibtex.readthedocs.io/en/latest/index.html>`_
for citations. To add citations in your mod documentation page:

- Add a ``<mod>.bib`` file under ``docs/source/refs``
- Add your new ``<mod>.bib`` file to the ``bibtex_bibfiles`` list in
  ``docs/source/conf.py``
- Use ``:footcite:t:`` or ``:footcite:p:`` as needed to cite references within
  your documentation
- Add the ``.. footbibliography::`` directive at the bottom of your
  documentation page to show references as footnotes
