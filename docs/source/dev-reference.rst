Developer Reference
===================

.. warning::
    This page is under construction

A place for detailed information on the tools/packages we use, to guide
developers who are working on mods.

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
