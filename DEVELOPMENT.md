# Dev guide

This is a short guide to setup a dev environment for mols2grid.

## Installation

1. Install conda or mamba
2. Create a new environment. Python 3.8+ (prefer 3.10):
   ```
   conda env create --name mols2grid --file docs/environment.yml
   ```
3. Install the package in editable mode:
   ```
   pip install -e .[dev]
   jupyter labextension develop . --overwrite
   ```

## Tests

To run tests locally:
- Install Firefox (needed for UI testing)
- Test your installation:
  ```
  pytest tests/
  ```
- You can select/skip the UI testing by specifying the `webdriver` mark in the pytest
  command: `-m webdriver` to select UI tests only, or `-m "not webdriver"` to skip them.

### Notebook test

The CI pipeline contains an additional test that runs a simple notebook (see [](tests/notebooks/))
to make sure that there are no issues with the widget and the different versions of JupyterLab and
IPywidgets. This test requires the version of the package to **NOT** be a prerelease and thus the
versions in [](package.json) and [](mols2grid/_version.py) should strictly follow the `X.Y.Z` scheme
(i.e. no `X.Y.Z-a1`).

If this pipeline fails, the CI will automatically produce an artifact named
`notebook-updated-snapshots` containing the captures that were produced while running the test.
Review these snapshots and if they seem correct, commit them to the
[](tests/notebooks/tests/mols2grid.test.ts-snapshots) directory.

## Formatting

We use `black` and `isort` for formatting so either install the corresponding extension
from your IDE or install the package with `pip install black isort`. The configuration
is done inside the `pyproject.toml` file.

## Pull requests

Making a pull request will automatically run the tests and documentation build for you.
Don't forget to update the `CHANGELOG.md` file with your changes.

## Release

For versioning, you'll have to update both `package.json` and `mols2grid/_version.py`
files (but this should be done by a maintainer directly).

The build and deployment process is run automatically when making a release on
GitHub.
To make a prerelease, bump the versions accordingly (`X.Y.Z-rc1` format) and run
the `build` GitHub Action manually.
