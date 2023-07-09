.PHONY: develop test

develop:
	python -m pip install pip setuptools wheel --upgrade
	python -m pip install -e .[examples]
	python -m pip install -r docs/requirements.txt
	python -m pip install sphinx-autobuild pre-commit
	pre-commit install

test:
	python -m unittest discover -b
	cd docs && make doctest
