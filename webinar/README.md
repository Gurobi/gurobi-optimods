# Simple one-pager

From this directory, start these two processes:

```
python -m jupyter notebook --config jupyter_notebook_config.py
python -m http.server -d pages
```

then navigate to `http://localhost:8000/empty.html`.

Open as a chrome app for a tidier view:

```
/Applications/Google\ Chrome.app/Contents/MacOS/Google\ Chrome --app='http://localhost:8000/empty.html'
```

# Over-complicated Quarto setup

## Requirements

Any virtualenv with

```
pip install jupyter gurobi-optimods[examples]
```

should do the trick. A full `requirements.txt` file is provided with a fixed
set of package versions (just a freeze from the above command in a clean python
3.11 environment).

You also need Quarto (https://quarto.org/) which can be installed with homebrew.

## Run

Start jupyter (config is needed to disable authentication & allow embedding):

```
python -m jupyter notebook --config jupyter_notebook_config.py
```

Start quarto (the port is important, jupyter is configured to allow embedding from localhost:8000):

```
quarto preview slides.qmd --port 8000 --no-browser
```

Open as a chrome app for a tidier view

```
/Applications/Google\ Chrome.app/Contents/MacOS/Google\ Chrome --app='http://localhost:8000'
```

## Editing

`slides.qmd` should be the only file that need additions, it is a quarto revealjs markdown document, see https://quarto.org/docs/presentations/revealjs/

## Todo

- Some formatting fixes are still needed (probably need to shrink the default font size)
- Zoom on the embedded pages sometimes appears off
