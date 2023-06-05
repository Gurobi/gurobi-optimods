Contributing to OptiMods
========================

We welcome contributions from the commumnity: bug fixes, doc corrections,
extensions to existing Mods, or development of new Mods. This page describes
conventions to follow when contributing, and the development process we follow
for the project, which is applicable to all contributors. For guidelines on
proposing and implementing new Mods, see :doc:`adding`.

Coding Standards
----------------

We use `black <https://github.com/psf/black>`_ to format all source code.
``pre-commit`` takes care of this formatting when installed (see
:ref:`contributing:development environment`) and runs checks automatically when
you commit changes. To run code quality checks manually, run the following and
ensure no issues are reported::

    pre-commit run --all-files
    python -m flake8 . --select=E9,F63,F7,F82,F811,F401 --show-source

We follow the `numpydoc <https://numpydoc.readthedocs.io/en/latest/format.html>`_
style guide for docstrings of all public API functions.

Testing
-------

All new additions, including bug fixes, Mod extensions, or brand new Mods `must`
include unit tests to verify their functionality. We use
:external+python:std:doc:`unittest <library/unittest>` to test Mod functionality
and :doc:`sphinx doctest <sphinx:usage/extensions/doctest>` to ensure code
snippets provided in the documentation work as expected.

Documentation
-------------

Documentation is a core part of the OptiMods project. All Mods must be
accompanied by a documentation page explaining their use and providing
appropriate background information. If new mod features are added, the
documentation must reflect the changes.

Development Environment
-----------------------

To set up your development environment:

1. Create and activate a Python 3.8 virtual environment using your preferred
   tool. We maintain compatibility with Python 3.8 and above, so developing
   using 3.8 locally will ensure your contributions are compatible.
2. Run ``make develop`` from the top-level of the repository. This command will
   install (using ``pip``) all packages required to run the unit tests, run the
   doc tests, and build the sphinx documentation. It will also install
   pre-commit hooks to enforce coding standards.

To run the tests:

1. Activate your virtual environment.
2. Run ``make test``. This command will run the unit tests of all mod
   implementations, and the doctests for examples in the documentation.

To build and view the docs:

1. Activate your virtual environment and change to the ``docs`` directory.
2. Run ``make livehtml``. This command will build the docs and open up a browser
   window at the index page.
3. Start editing the docs! Any changes to the documentation source files will
   automatically trigger a rebuild of the docs and a refresh of your browser
   window.

Developer Workflow
------------------

We use issues and pull requests to manage and review contributions.

- Minor doc fixes can be submitted directly as pull requests.
- Bugs should first be reported as issues before submitting a pull request to
  fix them. Use the 'Bug report' issue template.
- New Mod proposals should first be submitted as issues for discussion, see
  :doc:`adding` for further details.
