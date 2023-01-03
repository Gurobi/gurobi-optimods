.PHONY: wheel sdist develop test

wheel:
	python -m build --wheel

sdist:
	python -m build --sdist

develop:
	python -m pip install pip setuptools wheel --upgrade
	python -m pip install -e .
	python -m pip install -r docs/requirements.txt
	python -m pip install sphinx-autobuild pre-commit
	pre-commit install

test:
	python -m unittest discover -b
	cd docs && python -m unittest discover -b
	cd docs && make doctest
