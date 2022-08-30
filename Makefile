.PHONY: wheel sdist develop test

wheel:
	python -m build --wheel

sdist:
	python -m build --sdist

develop:
	python -m pip install -e .
	python -m pip install -r docs/requirements.txt
	python -m pip install sphinx-autobuild

test: Makefile
	python -m unittest discover -b
	cd docs && python -m unittest discover -b
	cd docs && make doctest
