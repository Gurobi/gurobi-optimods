# Contributing Guide

## Bug fixes

Create an issue first using the 'Bug Report' template.

## Propose a new nup

Create an issue, and use the 'New Nup Proposal' template to gather the required details.

## Implement a nup

1. Assign yourself to an open nup issue
2. Fork the repository to your own public github account
3. Follow the scaffolding instructions below to set up the required file structure
4. Create a new branch, commit the scaffolding files
5. Push the branch to your fork
6. Create a draft PR for your branch into nupstup/main (the implementation does not have to be complete to create a PR).
7. Follow the PR checklist.
8. When complete, take the PR out of draft and assign a reviewer.

### Scaffolding

To get started implementing a new nup, create the following files (replacing `<nup>` and `<category>` as needed). Each of the files below has an equivalent under the `_templates` directory which provides further details.

* `src/nupstup/<nup>.py`: Implementation of the nup goes here. This must be an importable module exposing a function or class which is entry-point for using the nup.
* `tests/test_<nup>.py`: Unit tests for the nup. Enough to exercise some simple cases and know the code is in a working state.
* `docs/examples/<nup>/__init__.py`: Empty file (this is needed for some testing machinery).
* `docs/examples/<nup>/nupstup.py`: A working example showing how to use the nup. This should be a simple python script that only reads or creates input data structures and runs the main nup function or class. The script shouldn't do anything with its results; this will live in the documentation file.
* `docs/examples/<nup>/gurobipy.py`: A working example which produces identical results to the `<nup>/nupstup.py` file above, but without using nupstup itself. This will serve as a comparison implementation in the documentation to show users how they would implement this nup themselves. It should not be a carbon copy of the nup code, just the bare essentials to get a result in the specific documented case.
* `docs/examples/test_examples.py`: Add a simple test for the examples. See the test_examples.py file for the pattern. This test should verify that nupstup.py and gurobipy.py produce the same result. It does so by importing both example scripts and running a simple check asserting the results are the same (or numerically close).
* `docs/source/library/<category>/<nup>.rst`: The documentation page for the nup. The example under `_templates` guides you through the specifics. In particular note the literalinclude blocks which will pull your example codes into the docs.
