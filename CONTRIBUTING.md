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
5. Include your mod in the toctree at `docs/source/mods/index.rst`
6. Create a new branch, commit the scaffolding files
7. Push the branch to your fork
8. Create a *draft* PR for your branch into gurobi-optimods/main
9. Follow the checklist in the PR
10. When complete, take the PR out of draft and assign a reviewer

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
