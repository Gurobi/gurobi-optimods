# Contributing Guide

## Bug fixes

Create an issue first using the 'Bug Report' template.

## Propose a new Mod

Create an issue, and use the 'New Mod Proposal' template to gather the required details.

## Implementing a new mod

1. Assign yourself to an open mod issue
2. Fork Gurobi/gurobi-optimods to your own github account
3. Set up your development environment (see below) and ensure that you can run the tests and build the docs
4. Copy the template files (replacing `<mod>` with your mod name):
```
cp _templates/mod.py src/gurobi_optimods/<mod>.py
cp _templates/test_mod.py tests/test_<mod>.py
cp _templates/mod.rst docs/source/mods/<mod>.rst
```
5. Include your mod in the toctree at `docs/source/mods/index.rst` (maintain lexicographic order)
6. Create a new branch, commit the scaffolding files
7. Push the branch to your fork
8. Create a *draft* PR for your branch into gurobi-optimods/main
9. Follow the checklist in the PR
10. When complete, take the PR out of draft and assign a reviewer

## Implementation notes

- Data files should live under `src/gurobi_optimods/data/<mod-name>` to reduce clutter
- Mods should be stateless (see the template example). This means gurobipy environments and models are created within a mod function, and closed before the function returns.
- We use `sphinxcontrib-bibtex` for citations/referencing. To add references for a mod:
  - Add a new `.bib` file under `docs/source/refs`
  - Add your new `.bib` file to the `bibtex_bibfiles` list in `docs/source/conf.py`
  - Use `:footcite:t:` or `:footcite:p:` to cite references within your documentation
  - Add the `.. footbibliography::` directive at the bottom of your documentation page to show references as footnotes (example in the L1 regression mod)
  - For more info, see [the docs](https://sphinxcontrib-bibtex.readthedocs.io/en/latest/index.html) for `sphinxcontrib-bibtex`

### Development environment

To set up your development environment:

1. Create and activate a python >=3.8 virtual environment
2. Run `make develop` from the top directory of this repository

To run the tests

1. Activate your virtual environment
2. Run `make test`. This will run the unit tests of the mod implementations and doctests examples from the docs

To build and view the docs:

1. Activate your virtual environment
2. Change to `docs` directory
3. Run `make livehtml`. This will build the docs and (after a little while) open up a browser window at the index page
4. Any change to the documentation source files will automatically rebuild the docs and trigger an update in your browser

### Motivation

The idea in one sentence: we will create open-source Python repository of
implemented optimization use cases, each with clear, informative, and pretty
documentation that explains how to use it, the mathematical model behind it,
and the implementation in code.

- A plethora of useful optimization models, in- and outside classical OR
- Data driven APIs
- Intuitive for Python users

A good mod:

1. provides background information for the topic and a formal statement
of the problem in the domain language of the target user
2. has a simple interface which shields the user from interacting with
gurobipy, a clean data-in data-out style using sensible data types, and demos
usage via a runnable example in the documentation
3. documentation presents results of an example using familiar packages:
