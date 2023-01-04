# Contributing Guide

## Bug fixes

Create an issue first using the 'Bug Report' template.

## Propose a new Mod

Create an issue, and use the 'New Mod Proposal' template to gather the required details.

## Implementing a new mod

1. Assign yourself to an open mod issue
2. Fork Gurobi/gurobi-optimods to your own github account
3. Set up your development environment (see below) and ensure that you can run the tests and build the docs
4. Follow the scaffolding instructions (see below) to set up the required files
5. Create a new branch, commit the scaffolding files
6. Push the branch to your fork
6. Create a *draft* PR for your branch into gurobi-optimods/main
7. Follow the checklist in the PR
8. When complete, take the PR out of draft and assign a reviewer

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

### Scaffolding

To get started implementing a new mod, create the following files (replacing `<mod>` as needed with the name of your mod). Each of the files below has an equivalent under the `_templates` directory which provides further details.

* `src/gurobi_optimods/<mod>.py`: Implementation of the mod goes here. This must be an importable module exposing a function or class which is entry-point for using the mod.
* `tests/test_<mod>.py`: Unit tests for the mod. Enough to exercise some simple cases and ensure that the code is in a working state.
* `docs/source/mods/<mod>.rst`: The documentation page for the mod. The example under `_templates` guides you through the specifics. It must include a working example showing how to use the mod. This should be a simple python script that only reads or creates input data structures and runs the main mod function or class.
* Include your mod in the toctree at `docs/source/mods/index.rst`
* Include your mod in the api docs at `docs/source/api.rst`
