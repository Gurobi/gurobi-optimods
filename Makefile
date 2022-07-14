.PHONY: wheel sdist develop test

wheel:
	python -m build --wheel

sdist:
	python -m build --sdist

develop:
	python -m pip install -e .
	python -m pip install -r docs/requirements.txt
	python -m pip install pytest

test: Makefile
	cd docs && python -m unittest -b examples.test_examples
	cd docs && make doctest
