This is a short guide to setup a dev environment for mols2grid.

1. Install conda or mamba
2. Create a new environment. Python 3.7+ (prefer 3.8):
   ```
   conda env create --name mols2grid --file docs/environment.yml
   ```
3. Install all the package dependencies in editable mode:
   ```
   pip install -e .[dev]
   ```

To run tests locally:
- Install Firefox (needed for UI testing)
- Test your installation:
  ```
  pytest tests/
  ```
- You can select/skip the UI testing by specifying the `webdriver` mark in the pytest
  command: `-m webdriver` to select UI tests only, or `-m "not webdriver"` to skip them.

We use `black` and `isort` for formatting so either install the corresponding extension
from your IDE or install the package with `pip install black isort`. The configuration
is done inside the `pyproject.toml` file.

Making a pull request will automatically run the tests and documentation build for you.
Don't forget to update the `CHANGELOG.md` file with your changes.

For versioning, you'll have to update both `package.json` and `mols2grid/_version.py`
files.

The build and deployment process is run automatically when making a release on
GitHub.
To make a prerelease, bump the versions accordingly (`X.Y.Z-rc1` format) and run
the `build` GitHub Action manually.
