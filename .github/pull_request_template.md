<!-- Provide a general summary of your proposed changes in the Title field above -->

### Description
<!-- Describe your changes in detail -->

For a new Mod contribution, please include the issue number for the Mod proposal and check through the list below.

### Checklist
<!-- go over following points. check them with an `x` if they are completed, (they turn into clickable checkboxes once the PR is submitted, so no need to do everything at once) -->

- Implementation:
  - [ ] Implementation of the Mod in the `gurobi_optimods` installable package
  - [ ] Tests for the Mod implementation in `tests/`
  - [ ] Docstrings for public API, correctly linked using sphinx-autodoc
- Example codes:
  - [ ] Standalone `<mod-name>.py` script in examples
  - [ ] Any necessary example data importable from `gurobipy_optimods.datasets` if not easily available elsewhere
  - [ ] Simple test in `test_examples.py` which imports both versions and checks that they match
- Documentation page:
  - [ ] Problem specification with domain tab and mathprog tab
  - [ ] Example of the input data format
  - [ ] Code example literalincluded into the docs
  - [ ] Solutions presented in doctests

**Have a nice day!**