<!-- Provide a general summary of your proposed changes in the Title field above -->

### Description
<!-- Describe your changes in detail -->

For a new Nup contribution, please include the issue number for the Nup proposal and check through the list below.

### Checklist
<!-- go over following points. check them with an `x` if they are completed, (they turn into clickable checkboxes once the PR is submitted, so no need to do everything at once) -->

- Implementation:
  - [ ] Implementation of the Nup in the `nupstup` installable package
  - [ ] Tests for the Nup implementation in `tests/`
  - [ ] Docstrings for public API, correctly linked using sphinx-autodoc
- Example codes:
  - [ ] Standalone `xxx_nupstup.py` script in examples
  - [ ] Standalone `xxx_gurobipy.py` script in examples
  - [ ] Any necessary example data importable from `nupstup.datasets` if not easily available elsewhere
  - [ ] Simple test in `test_examples.py` which imports both versions and checks that they match
- Documentation page:
  - [ ] Problem specification with domain tab and mathprog tab
  - [ ] Example of the input data format
  - [ ] Code example with nupstup and gurobipy tab
  - [ ] Solutions presented in doctests

**Have a nice day!**