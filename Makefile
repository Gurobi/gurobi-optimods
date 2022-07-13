.PHONY: test

test: Makefile
	cd docs/examples && python -m pytest
	cd docs && make doctest
	cd docs && make html
