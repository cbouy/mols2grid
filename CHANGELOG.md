# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

---

## [2.0.0] - 2023/07/23

This release is a major change on the UI contributed by @themoenen, refer to
[PR#55](https://github.com/cbouy/mols2grid/pull/55) for the full list of changes:

### Added
- `background_color="white"` parameter added to `display` and `save` to control the background
  color of each cell.
- Property values displayed via subset or tooltip can now be copied by clicking them.

### Changed
- Responsiveness: the grid as well as all UI components are now fully responsive, up until around
  ~415px wide. Smaller usage seems unlikely.
- `n_items_per_page` replaces the `n_rows` and `n_cols` parameters. The responsive CSS assumes
  results to be a multiple of 12, otherwise a gap is displayed at the end of the grid.
- Hover tooltips: now displayed when hovering the *`i`* icon, and anchored by clicking this icon.
- Save CVS: When exporting as CSV we now use a semicolon `;` as delineator instead of a tab, this
  way CSVs are properly previewed on Mac.
- Templates: the `pages`/`table` templates and corresponding functions have been renamed to
  `interactive`/`static` for clarity.
- Parameters `width` and `height` have been renamed to `iframe_width` and `iframe_height` to be more
  descriptive.
- Improved the sorting UI to be more intuitive.
- You can now toggle between text/SMARTS search instead of having to click a dropdown.
- The checkbox menu icon was replaced with a more standard triple dot menu.
- The mols2grid-id is now permanently displayed next to the checkbox in the corner.
- SVGs are now rendered with a transparent background.
- The entire cell is now clickable, instead of just the tiny checkbox.
- Implemented some basic keyboard navigation: ENTER / ESC for selecting / unselecting, arrows and
  TAB for navigation.
- Copy to clipboard: this will now record a tab delineated text which is ready to be pasted into a
  spreadsheet, instead of the less useful dictionary format.
- The `tooltip_trigger` parameter has been removed since callback and hover functionalities don't
  overlap anymore.
- We're now automatically adding `"img"` to the subset instead of throwing an error.
- A smaller default font size (12px) makes the experience out of the box a bit more practical.

### Fixed
- Docstring now mention default values and some inconsistencies have been fixed.
- All UI elements are now neatly aligned and displayed on top so they're accessible without
  scrolling.
- Longer property values are now only truncated in the interactive grid, and broken into multiple
  lines by default in the static grid, as it is mostly meant for printing. A new parameter
  `truncate` lets you override the default truncating behavior.
- The tooltip display zone (around which the tooltip is displayed) is now the entire cell instead
  of only the image, so it never overlaps with any of the cell's data or functionality.
- When you download a CSV or SMILES file without any cells selected, you will now download data for
  all cells instead of an empty document.
- Parameter `gap` in static template didn't work.
- We no longer resize the iframe if a custom `iframe_height` is set by the display function.

---

## [1.1.1] - 2023/03/18

### Added
- Support for `pathlib.Path` objects as input for `display`, `save`, `MolGrid.from_sdf`
  and `sdf_to_dataframe`.

### Changed
- The hover tooltip placement has been changed from `"bottom"` to `"auto"`.
- Code and notebook formatting with `black` and `isort`.
- Switched to `hatchling` for the build process, `tbump` for versioning, and migrated to
  using only the `pyproject.toml` file.
- Refactored tests to use Pytest's `contest.py` file.

### Fixed
- CSV export when sorting the grid was not using the selected molecules.

---

## [1.1.0] - 2022/12/24

### Added
- Predefined JavaScript callbacks in the `mols2grid.callbacks` module. Those can be
  extensively configured:
  - `info`: displays a bigger image alongside some common descriptors for the molecule
  - `show_3d`: displays the molecule in 3D
  - `external_link`: opens a new tab. By default, opens [Leruli.com](https://leruli.com/)
    using the SMILES of the molecule.
- Support for `tuple` of molecules in `display` and `save`.

### Changed
- The `"click"` event is now automatically removed from `tooltip_trigger` when
  specifying a callback.

### Fixed
- Text searches containing any of the following regex characters `-[]{}()*+?.,\^$|#`
  would automatically return an empty grid, preventing searching for CAS numbers and any
  other identifier or text containing the above characters. This has been temporarily
  patched until a proper fix is released in the underlying `list.js` library.
- The link to the KNIME component on the corresponding badges has been fixed.


## [1.0.0] - 2022/09/04
### Added
- Notebooks running in VSCode and Jupyter Lab now support accessing selections from
  Python, executing Python callback functions, and filtering based on other widgets.
### Changed
- Python callbacks can now also be `lambda` functions.
- If ``prerender=True``, substructure highlighting will be automatically disabled by
  default instead of raising an error.
- When exporting a selection to a SMILES file through the GUI, the output no longer
  contains a header.
- Relies on a custom ipywidget to handle communication between the front-end/Javascript
  and the back-end/Python.
- When calling `grid.filter` and other filtering methods, mols2grid will now use the
  filtering code based on ipywidgets, except for Streamlit where it will use the older
  JavaScript version of the code to maintain compatibility.
### Fixed
- Automatically fitting to the content's height in Streamlit.
### Removed
- `mapping` argument for renaming fields, replaced by `rename` in `v0.1.0`.
- `mols2grid.selection`, replaced by `mols2grid.get_selection()` in `v0.1.0`.


## [0.2.4] - 2022/05/29
### Fixed
- Calling `MolGrid.get_selection()` when 2 grids with different names are present should
  now display the selection of the grid itself, and not the selection corresponding to
  indices of the grid that was last interacted with.


## [0.2.3] - 2022/05/10
### Fixed
- Doing a substructure search on molecules with explicit hydrogens should now highlight
  the correct atoms.


## [0.2.2] - 2022/04/04
### Added
- A proper documentation page with tutorials can now be accessed online.
- Added a `single_highlight=False` parameter to only highlight a single match per
  molecule in substructure queries.
- Added a *"Check matching"* button that only selects items that match the current search
  and/or filters.
- Added `custom_css`, `custom_header` and `sort_by` to the "table" template
### Changed
- Compounds matching a substructure search are now aligned to the query molecule before
  rendering the image.
- When doing a substructure search, all matches are now highlighted by default. To only
  show a single one, use `single_highlight=True`.
- The *Check all*, *Uncheck all* and *Invert* selection buttons have been fixed. They now
  actually check/uncheck ALL items, and not just the ones matching the current search. A
  *Check matching* button has been added to reproduce the old behaviour.
- If both `subset` and `tooltip` are `None`, the index and image will be directly
  displayed on the grid while the remaining fields will be in the tooltip. This makes the
  default representation much more readable.
- The default number of columns is now 5 for `template="table"` (same as the other default
  template)
### Fixed
- `template="table"` now correctly displays images when `prerender=True` (Issue #27)
- Displaying the grid with `template="table"` in a notebook now automatically fits to the
  content of the table.


## [0.2.1] - 2022/02/23
### Fixed
- Field names containing spaces are now correctly delt with
- The text search now looks for matches inside the values of the tooltip fields, rather
  than inside the HTML code of the tooltip which included tags and other irrelevant text
- Fixed an encoding bug when saving the grid as an HTML file on French Windows, which uses
  CP-1252 encoding instead of UTF-8


## [0.2.0] - 2022/02/10
### Added
- `cache_selection=True` allows to retrieve the checkbox state when re-displaying a grid,
  as long as they have the same name. Fixes #22
- `prerender=False` moves the rendering of molecule images from Python to the browser and
  only when the molecule is on the current page, giving a performance boost and allowing 
  to process much larger files. Fixes #17
- `substruct_highlight=True` highlight the atoms that matched the substructure query when
  using the SMARTS search (only available when `prerender=False`). Fixes #18
- Added CSV save option. Exports all the data present in `subset` and `tooltip` for the
  current selection
- Support for `.sdf.gz` files
- Added automated tests of the interface, which should prevent future updates from
  breaking things
### Changed
- Python 3.6 is no longer supported
- Molecule images are now generated by the web browser (see `prerender=False` argument)
- The coordinates of the input file are now ignored by default (`use_coords=False`). This
  change was made to comply with generating images from SMILES string with the browser by
  default.
- Python callbacks are now automatically registered in Google Colab
- Javascript callbacks can access RDKit as either `RDKit` or `RDKitModule`
- The "img" field is now available from the callback data
- The `subset` parameter now throws an error if "img" is not present
- Clicking "Check all"/"Uncheck all" should now be faster
- Bumped RDKit JS version to `2021.9.4` to better support moldrawoptions
- Installation now requires `jinja2>=2.11.0` to prevent an error when given a pathlib.Path
  object instead of a string
### Fixed
- Callbacks now work when `selection=False`. Fixes: Issue #22
- Using both `transform` and `style` should now display the labels as expected in the
  tooltip
- Fixed a race condition when clicking checkboxes on different grids
- Fixed the `gap` argument not being properly taken into account
- Automatic resizing of the iframe (used in `mols2Grid.display`) should now work even
  better


## [0.1.0] - 2021/10/11
### Added
- The grid can be filtered using pandas DataFrame's `query` and `loc` logic (mostly
  useful to combine with ipywidgets) with `MolGrid.filter_by_index` and `MolGrid.filter`.
- Selections can now be modified (select and unselect all, or invert) and exported (to
  clipboard or a SMILES file) even without a notebook kernel. Fixes: Issue #16.
- The grid can be sorted according to the selection status and to values in the tooltips.
- Added tracking the selection in multiple grids at the same time (i.e. it's not a
  global object that get's overwritten anymore).
- Added support for executing custom JavaScript code or Python function when clicking on
  a molecule's image through the `callback` argument.
- Added the `mols2grid.make_popup_callback` helper function to create a popup-like window
  as a JavaScript callback.
- Added styling for the whole cell through `style={"__all__": userfunction}`.
- Added `mols2grid.get_selection()` allowing users to specify which grid selection should
  be returned. Without argument, the most recently updated grid is returned.
- Added `mols2grid.list_grids()` to return a list of grid names available.
- Added the `mols2grid.sdf_to_dataframe` function to easily convert an SDF file to a
  pandas DataFrame.
- Added the `custom_css` argument to pass custom CSS for the HTML document.
- Added the `sort_by` argument to change how the grid elements are ordered
### Changed
- The functions in `style` and `transform` are now also applied to tooltips.
- The sizing of the iframe displaying the grid is now fully automated and more precise.
- Reorganized the code to separate the JS, CSS and HTML templates.
### Fixed
- Fixed `mols2grid.save` that returned an error about missing the `output` argument.
- The tooltip is now compatible with the "focus" mode: `tooltip_trigger="focus"`.
- Fixed rendering SVG images in tooltips.
### Deprecated
- Deprecated `mols2grid.selection` in favor of `mols2grid.get_selection()`.
- Deprecated `mapping` in favor of `rename` in the MolGrid class and `mols2grid.display`.


## [0.0.6] - 2021/05/14
### Changed
- Javascript module for RDKit is now sourced from `unpkg.com` and pinned to `v2021.3.2`


## [0.0.5] - 2021/04/08
### Added
- New `transform` parameter that accepts a dictionary of field-function items where each
  function transforms the input value that will be displayed. Fixes: Issue #10
### Fixed
- Running mols2grid could throw an ImportError (instead of ModuleNotFoundError) if the
  `google` module was installed, but not `google.colab`. Solved by PR #11
- Private molecule properties (i.e properties starting with `_`) were not registered when
  reading properties from RDKit molecules (SDF or list of mols).


## [0.0.4] - 2021/04/01
### Changed
- The demo notebook can now be run on Google Colab
### Fixed
- DataFrames with `NaN` values would previously lead to an **empty grid** as `NaN` were
  converted to `nan` (not recognized by JS) instead of `NaN`.
- Selection of molecules in Google Colab now works as expected.
- Saved documents are now displayed properly


## [0.0.3] - 2021/03/31
### Added
- **SMARTS search**: the "🔎" button now lets users choose between a text search or a
  SMARTS search. This relies on RDKit's "MinimalLib" JS wrapper which is still in beta,
  and will likely break at some point. Use at your own risk!
- **Sorting**: added a "Sort by" button that lets users choose in which order the
  molecules should be listed. Default: by index. Fixes: Issue #7
- `MolDrawOptions` **drawing** parameter: this will allow further customization of the
  drawing options.
- **Selection**: added checkboxes to each cell. Clicking on a checkbox will add the
  molecule's corresponding index and SMILES to the `mols2grid.selection` dictionary.
- New **input** formats: dict and record (list of dicts) are automatically converted to a
  pandas DataFrame when used as input to the MolGrid class. The `mols2grid.display`
  function only accepts the dict option (since the list format is already used for lists
  of RDKit molecules).
- New **input** options: `mol_col` parameter. Adds the ability to directly use an RDKit
  mol instead of relying on a SMILES intermediate. This makes using the 2D coordinates of
  the input mol a possibility, instead of systematically generating new ones. It also
  allows for adding annotations and highlights on drawings. Introduces 2 new parameters:
  - `mol_col=None`: Column of the dataframe containing RDKit molecules
  - `use_coords=True`: directly use the coordinates from each molecule, or generate new
    ones
### Changed
- The "mols2grid-id" now keeps track of molecules that could not be read by RDKit. This
  makes relying on the index of the corresponding entry more reliable.
### Fixed
- The "mols2grid-id" field is now correctly set in the internal `DataFrame` and the
  JavaScript `List`.
- Using the search bar will now only search inside the fields listed in `subset` and
  `tooltip` and exclude the `img` field.
- When using the `display` function, the `height` of the iframe is now automatically set
  based on the different parameters, instead of a fixed 600px height. Fixes: Issue #6


## [0.0.2] - 2021/03/23
- First release
