# Building the AFRC documentation

The docs are built with [Sphinx](https://www.sphinx-doc.org/) and the Read the Docs theme. Install the requirements (Sphinx comes in with the theme):

```bash
pip install -r requirements.txt
```

Then, from this directory, build the HTML pages:

```bash
make html
```

The output goes to `_build/html/`; open `index.html` there to view it. To treat warnings as errors (as we do before a release), run `sphinx-build -W -b html . _build/html` instead.
