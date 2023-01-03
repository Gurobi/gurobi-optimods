# Contributing Guide

## Bug fixes

Create an issue first using the 'Bug Report' template.

## Propose a new Mod

Create an issue, and use the 'New Mod Proposal' template to gather the required details.

## Implement a mod

1. Assign yourself to an open mod issue
2. Fork the repository to your own public github account
3. Follow the scaffolding instructions below to set up the required file structure
4. Create a new branch, commit the scaffolding files
5. Push the branch to your fork
6. Create a draft PR for your branch into gurobi-optimods/main (the implementation does not have to be complete to create a PR).
7. Follow the PR checklist.
8. When complete, take the PR out of draft and assign a reviewer.

### Development environment

To set up your development environment:

1. Create a python 3.9 virtual environment
2. Run `make develop` from the top directory of this repository

To build and view the docs (note that when building locally, there is no theme applied, so it won't look like the readthedocs page):

1. Activate your virtual environment
2. Change to `docs` directory
3. Run `make livehtml`. This will build the docs and (after a little while) open up a browser window at the index page
4. Any change to the documentation source files will automatically rebuild and update your browser

To run the tests

1. Activate your virtual environment
2. Run `make test`. This will run the tests of the mod implementations and examples in the docs

### Scaffolding

To get started implementing a new mod, create the following files (replacing `<mod>` and `<category>` as needed). Each of the files below has an equivalent under the `_templates` directory which provides further details.

* `src/gurobi_optimods/<mod>.py`: Implementation of the mod goes here. This must be an importable module exposing a function or class which is entry-point for using the mod.
* `tests/test_<mod>.py`: Unit tests for the mod. Enough to exercise some simple cases and know the code is in a working state.
* `docs/source/mods/<mod>.rst`: The documentation page for the mod. The example under `_templates` guides you through the specifics. It should included a working example showing how to use the mod. This should be a simple python script that only reads or creates input data structures and runs the main mod function or class.
