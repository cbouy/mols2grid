This is a short guide to setup a dev environment for mols2grid.

1. Install [uv](https://docs.astral.sh/uv/)
2. Create a new environment and install:
   ```
   uv sync --python 3.11
   ```

We use [poethepoet](https://poethepoet.natn.io/) to define tasks to run in development environments.
You can run all of the checks detailed below using the command
```
uv run poe check
```

You can also get a list of available checks with:
```
uv run poe --help
```

  a. Running tests

To run the test suite, simply execute:
```
uv run poe tests
```
You can select/skip the UI testing by specifying the `webdriver` mark in the pytest
  command: `-m webdriver` to select UI tests only, or `-m "not webdriver"` to skip them.

  b. Building the documentation

Building the HTML files for the documentation and tutorials can be done with the following command:
```
uv run poe docs
```
You can then open the `docs/_build/html/index.html` file with your browser to navigate the docs and
see any changes that you've made.

If you're adding a new module, you will need to update some `.rst` files in the `docs/api/`
folder.

You will find the tutorials notebooks in the `docs/notebooks/` section. These are Jupyter-notebook
files with a twist for Markdown cells: you can use the
[MyST syntax](https://myst-nb.readthedocs.io/en/latest/authoring/jupyter-notebooks.html#syntax)
to format the content. See the [Authoring section](https://myst-parser.readthedocs.io/en/latest/syntax/typography.html)
for more details.

  c. Code formatting and linting

You can check if your code complies with our code style standards with the following command:
```
uv run poe style-check
```

You can automatically format your changes to match with the style used in this project, as well as
fixing any lint errors (unused imports, type annotations...etc.) with the following command:

```
uv run poe style-fix
```

Making a pull request will automatically run the tests and documentation build for you.
Don't forget to update the `CHANGELOG.md` file with your changes.

For versioning, you'll have to update both `package.json` and `mols2grid/_version.py`
files.

The build and deployment process is run automatically when making a release on
GitHub.
To make a prerelease, bump the versions accordingly (`X.Y.Z-rc1` format) and run
the `build` GitHub Action manually.
