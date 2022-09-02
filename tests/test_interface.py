import os
import json
import pytest
from flaky import flaky
from ast import literal_eval
from base64 import b64encode, b64decode
from types import SimpleNamespace
from selenium import webdriver
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.common.by import By
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.common.action_chains import ActionChains
from selenium.common.exceptions import (NoSuchElementException,
                                        StaleElementReferenceException)
import geckodriver_autoinstaller
from rdkit import __version__ as rdkit_version
from rdkit import Chem, RDConfig
from rdkit.Chem import AllChem
from PIL import Image
from io import BytesIO
from hashlib import md5
from cairosvg import svg2png
import mols2grid
from mols2grid.utils import env
from mols2grid.select import register

geckodriver_autoinstaller.install()
pytestmark = pytest.mark.webdriver
GITHUB_ACTIONS = os.environ.get("GITHUB_ACTIONS")

HEADLESS = True
PAGE_LOAD_TIMEOUT = 10

class selection_available:
    def __init__(self, is_empty=False):
        self.empty = is_empty

    def __call__(self, driver):
        sel = driver.execute_script("return SELECTION.to_dict();")
        sel = literal_eval(sel)
        if sel == {} and self.empty:
            return True
        elif sel != {} and not self.empty:
            return sel
        return False

class FirefoxDriver(webdriver.Firefox):
    def wait_for_img_load(self, max_delay=5, selector="#mols2grid .cell .data-img svg"):
        return (WebDriverWait(self, max_delay)
                .until(EC.presence_of_element_located((By.CSS_SELECTOR,
                                                       selector))))

    def wait_for_selection(self, is_empty=False, max_delay=3):
        return (WebDriverWait(self, max_delay)
                .until(selection_available(is_empty)))

    def wait(self, condition, max_delay=5):
        return (WebDriverWait(self, max_delay,
                ignored_exceptions=[StaleElementReferenceException])
                .until(condition))

    def find_by_id(self, element_id, **kwargs):
        condition = EC.presence_of_element_located((By.ID, element_id))
        return self.wait(condition, **kwargs)

    def find_clickable(self, by, selector, **kwargs):
        condition = EC.element_to_be_clickable((by, selector))
        return self.wait(condition, **kwargs)

    def find_by_css_selector(self, css_selector, **kwargs):
        condition = EC.presence_of_element_located((By.CSS_SELECTOR, css_selector))
        return self.wait(condition, **kwargs)

    def find_by_class_name(self, name, **kwargs):
        condition = EC.presence_of_element_located((By.CLASS_NAME, name))
        return self.wait(condition, **kwargs)

    def find_all_by_class_name(self, name, **kwargs):
        condition = EC.presence_of_all_elements_located((By.CLASS_NAME, name))
        return self.wait(condition, **kwargs)

    def get_svg_md5_hash(self, selector="#mols2grid .cell .data-img"):
        img = self.find_by_css_selector(selector)
        im = svg2png(bytestring=(img.get_attribute("innerHTML")))
        im = Image.open(BytesIO(im))
        return md5(im.tobytes()).hexdigest()

    def get_png_md5_hash(self, selector="#mols2grid .cell .data-img *"):
        img = self.find_by_css_selector(selector)
        im = Image.open(BytesIO(b64decode(img.get_attribute("src")[22:])))
        return md5(im.tobytes()).hexdigest()

    def substructure_query(self, smarts):
        self.find_clickable(By.ID, "searchBtn").click()
        self.find_clickable(By.ID, "smartsSearch").click()
        self.find_by_id("searchbar").send_keys(smarts)
        self.wait_for_img_load()

    def text_search(self, txt):
        self.find_clickable(By.ID, "searchBtn").click()
        self.find_clickable(By.ID, "txtSearch").click()
        self.find_by_id("searchbar").send_keys(txt)
        self.wait_for_img_load()

def determine_scope(fixture_name, config):
    if GITHUB_ACTIONS:
        return "function"
    return "module"

@pytest.fixture(scope=determine_scope)
def driver():
    options = webdriver.FirefoxOptions()
    options.headless = True if GITHUB_ACTIONS else HEADLESS
    driver = FirefoxDriver(options=options)
    driver.set_page_load_timeout(PAGE_LOAD_TIMEOUT)
    yield driver
    driver.quit()

sdf_path = f"{RDConfig.RDDocsDir}/Book/data/solubility.test.sdf"

def get_grid(df, **kwargs):
    kwargs.setdefault("mol_col", "mol")
    return mols2grid.MolGrid(df, **kwargs)

def get_doc(grid, kwargs):
    html = grid.render(**kwargs)
    html = b64encode(html.encode()).decode()
    return "data:text/html;base64,{}".format(html)

@pytest.fixture(scope="module")
def df():
    return mols2grid.sdf_to_dataframe(sdf_path)[:30]

@pytest.fixture(scope="module")
def mols(df):
    return df["mol"][:5]

@pytest.fixture(scope="module")
def grid(df):
    return get_grid(df)

@pytest.fixture(scope="module")
def html_doc(grid):
    return get_doc(grid, dict(n_rows=1, subset=["_Name", "img"]))

# make sure non-parametrized test is ran first
@pytest.mark.order(1)
def test_no_subset_all_visible(driver, grid):
    doc = get_doc(grid, {"tooltip": [], "selection": False})
    driver.get(doc)
    columns = set(grid.dataframe.columns.drop("mol").to_list())
    cell = driver.find_by_css_selector("#mols2grid .cell")
    data_el = cell.find_elements_by_class_name("data")
    classes = [c.replace("data-", "").replace("-copy", "")
               for x in data_el
               for c in x.get_attribute("class").split(" ")
               if c.startswith("data-")]
    classes = set(classes)
    assert classes == columns

def test_smiles_hidden(driver, html_doc):
    driver.get(html_doc)
    el = driver.find_by_css_selector("#mols2grid .cell .data-SMILES")
    assert not el.is_displayed()

@pytest.mark.parametrize("page", [1, 2, 3])
def test_page_click(driver, grid, page):
    doc = get_doc(grid, dict(subset=["img", "_Name"], n_cols=3, n_rows=3))
    driver.get(doc)
    for i in range(2, page+1):
        driver.wait_for_img_load()
        next_page = driver.find_by_css_selector(f'a.page-link[data-i="{i}"]')
        next_page.click()
    first_cell = driver.find_by_class_name('cell')
    mols2grid_id = 9 * (page - 1)
    name = first_cell.find_element_by_class_name('data-_Name')
    ref = grid.dataframe.iloc[mols2grid_id]
    assert name.text == ref["_Name"]

@pytest.mark.parametrize(["name", "css_prop", "value", "expected"], [
    ("gap", "margin-top", 20, "20px"),
    ("border", "border-top-width", "3px solid", "3px"),
    ("border", "border-top-style", "1px dashed", "dashed"),
    ("border", "border-top-color", "1px solid blue", "rgb(0, 0, 255)"),
    ("fontsize", "font-size", "16pt", "21.3333px"),
    ("fontfamily", "font-family", "Consolas", "Consolas"),
    ("textalign", "text-align", "right", "right"),
    ("custom_css", "background-color", ".cell { background-color: black; }", "rgb(0, 0, 0)"),
])
def test_css_properties(driver, grid, name, css_prop, value, expected):
    doc = get_doc(grid, {name: value})
    driver.get(doc)
    computed = driver.execute_script(f"return getComputedStyle(document.querySelector('#mols2grid .cell')).getPropertyValue({css_prop!r});")
    assert computed == expected

def test_text_search(driver, html_doc):
    driver.get(html_doc)
    driver.wait_for_img_load()
    driver.text_search("iodopropane")
    el = driver.find_by_css_selector("#mols2grid .cell .data-SMILES")
    assert el.get_attribute("innerHTML") == "CC(I)C"

@flaky(max_runs=3, min_passes=1)
def test_smarts_search(driver, html_doc):
    driver.get(html_doc)
    driver.wait_for_img_load()
    driver.substructure_query("CC(I)C")
    el = driver.find_by_css_selector("#mols2grid .cell .data-_Name")
    assert el.text == "2-iodopropane"

def test_selection_click(driver, html_doc):
    driver.get(html_doc)
    driver.wait_for_img_load()
    driver.find_clickable(By.CSS_SELECTOR, "input[type='checkbox']").click()
    sel = driver.wait_for_selection(is_empty=False)
    assert sel == {0: "CCC(C)CC"}
    register._clear()

def test_selection_with_cache_check_and_uncheck(driver, df):
    register._init_grid("cached_sel")
    event = SimpleNamespace(new='{0: "CCC(C)CC"}')
    register.selection_updated("cached_sel", event)
    grid = get_grid(df, name="cached_sel", cache_selection=True)
    doc = get_doc(grid, {})
    driver.get(doc)
    driver.wait_for_img_load()
    sel = driver.wait_for_selection(is_empty=False)
    assert sel == {0: "CCC(C)CC"}
    driver.find_clickable(By.CSS_SELECTOR, "input[type='checkbox']").click()
    empty_sel = driver.wait_for_selection(is_empty=True)
    assert empty_sel is True
    register._clear()

def test_selection_check_uncheck_invert(driver, html_doc):
    driver.get(html_doc)
    driver.wait_for_img_load()
    # search
    driver.text_search("iodopropane")
    # check all
    driver.find_clickable(By.ID, "chkboxDropdown").click()
    driver.find_clickable(By.ID, "btn-chkbox-all").click()
    sel = driver.wait_for_selection(is_empty=False)
    assert len(sel) == 30
    # uncheck all
    driver.find_clickable(By.ID, "chkboxDropdown").click()
    driver.find_clickable(By.ID, "btn-chkbox-none").click()
    empty_sel = driver.wait_for_selection(is_empty=True)
    assert empty_sel is True
    # check matching
    driver.find_clickable(By.ID, "chkboxDropdown").click()
    driver.find_clickable(By.ID, "btn-chkbox-match").click()
    sel = driver.wait_for_selection(is_empty=False)
    assert sel == {27: 'CC(I)C'}
    # invert
    driver.find_clickable(By.ID, "chkboxDropdown").click()
    driver.find_clickable(By.ID, "btn-chkbox-invert").click()
    sel = driver.wait_for_selection(is_empty=False)
    assert len(sel) == 29
    register._clear()

@pytest.mark.parametrize("prerender", [True, False])
def test_image_size(driver, df, prerender):
    grid = get_grid(df, size=(100, 100), prerender=prerender)
    doc = get_doc(grid, {"selection": False, "border": "0", "substruct_highlight": False})
    driver.get(doc)
    if not prerender:
        driver.wait_for_img_load()
    img = driver.find_by_css_selector("#mols2grid .cell .data-img *")
    assert img.size == {"height": 100.0, "width": 100.0}

def test_image_use_coords(driver, df):
    mols = df["mol"][:1]
    AllChem.EmbedMolecule(mols[0], randomSeed=0xf00d)
    grid = mols2grid.MolGrid.from_mols(
        mols, use_coords=True, prerender=True, useSVG=False)
    doc = get_doc(grid, {"substruct_highlight": False})
    driver.get(doc)
    md5_hash = driver.get_png_md5_hash()
    if rdkit_version == "2020.03.1":
        assert md5_hash == "aed60ed28347831d24f02dbb5be19007"
    else:
        assert md5_hash == "80ac785821c6aee7923c793930c4723e"

@pytest.mark.parametrize(["coordGen", "prerender", "expected"], [
    (True,  True, "acf7cf7cd5cfaa5eaf4a4f257e290e49"),
    (False, True, "34a44eef3c0d2af6503069fc068e705d"),
    (True, False, "d8227a9948b6f801193085e5a8b257a2"),
    (False, False, "a43327f2cfcff45a2aa87217cb658eb8"),
])
def test_coordgen(driver, mols, coordGen, prerender, expected):
    useSVG = not prerender
    grid = mols2grid.MolGrid.from_mols(mols,
        coordGen=coordGen, prerender=prerender, useSVG=useSVG, use_coords=False)
    doc = get_doc(grid, {"substruct_highlight": False})
    driver.get(doc)
    if not prerender:
        driver.wait_for_img_load()
    if useSVG:
        md5_hash = driver.get_svg_md5_hash()
    else:
        md5_hash = driver.get_png_md5_hash()
    assert md5_hash == expected

@pytest.mark.parametrize(["removeHs", "prerender", "expected"], [
    (True, True, "acf7cf7cd5cfaa5eaf4a4f257e290e49"),
    (False, True, ("2d30379732d312409d1d93ad5b3c3e60"
                   if rdkit_version == "2020.03.1"
                   else "b163a9cb19caee002ff9d739ed513672")),
    (True, False, "d8227a9948b6f801193085e5a8b257a2"),
    (False, False, "a8d830b6127c390dd5759a6ab1af8b3a"),
])
def test_removeHs(driver, df, removeHs, prerender, expected):
    useSVG = not prerender
    mol = df["mol"][0]
    mol.ClearProp("SMILES")
    mols = [Chem.AddHs(mol)]
    grid = mols2grid.MolGrid.from_mols(
        mols, removeHs=removeHs, prerender=prerender, useSVG=useSVG, use_coords=False)
    doc = get_doc(grid, {"n_rows": 1, "substruct_highlight": False})
    driver.get(doc)
    if not prerender:
        driver.wait_for_img_load()
    if useSVG:
        md5_hash = driver.get_svg_md5_hash()
    else:
        md5_hash = driver.get_png_md5_hash()
    assert md5_hash == expected

@pytest.mark.parametrize(["kwargs", "expected"], [
    (dict(addAtomIndices=True), "e01cd9295b091e9e2167568bd3cf139d"),
    (dict(fixedBondLength=10), "1dfa2a5c43022ce66a2ff26a517698cf"),
    (dict(atomColourPalette={6: (0, .8, .8)}), "91d114d421fc6a05e4537d19ad4b8324"),
    (dict(legend="foo"), "49cfe8038b8fc0b0401fd0056426d6b7"),
])
def test_moldrawoptions(driver, df, kwargs, expected):
    grid = get_grid(df, **kwargs)
    doc = get_doc(grid, dict(n_rows=1, n_cols=1, subset=["img"]))
    driver.get(doc)
    driver.wait_for_img_load()
    md5_hash = driver.get_svg_md5_hash()
    assert md5_hash == expected

def test_hover_color(driver, grid):
    doc = get_doc(grid, {"hover_color": "red"})
    driver.get(doc)
    (ActionChains(driver)
        .move_to_element(driver.find_by_css_selector("#mols2grid .cell"))
        .perform())
    color = driver.execute_script(f"return getComputedStyle(document.querySelector('#mols2grid .cell')).getPropertyValue('background-color');")
    assert color == "rgb(255, 0, 0)"

@flaky(max_runs=3, min_passes=1)
def test_tooltip(driver, grid):
    doc = get_doc(grid, {"tooltip": ["_Name"]})
    driver.get(doc)
    driver.wait_for_img_load()
    (ActionChains(driver)
        .move_to_element(driver.find_by_css_selector("#mols2grid .cell .data-img *"))
        .perform())
    tooltip = driver.find_by_css_selector('div.popover[role="tooltip"]')
    el = tooltip.find_element_by_class_name("popover-body")
    assert el.get_attribute("innerHTML") == "<strong>_Name</strong>: 3-methylpentane"

@flaky(max_runs=3, min_passes=1)
def test_tooltip_trigger(driver, grid):
    doc = get_doc(grid, {"tooltip": ["_Name"], "tooltip_trigger": "click"})
    driver.get(doc)
    driver.wait_for_img_load()
    (ActionChains(driver)
        .move_to_element(driver.find_by_css_selector("#mols2grid .cell .data-img *"))
        .perform())
    with pytest.raises(NoSuchElementException):
        driver.find_element_by_css_selector('div.popover[role="tooltip"]')
    driver.find_clickable(By.CSS_SELECTOR, "#mols2grid .cell .data-img *").click()
    tooltip = driver.find_by_css_selector('div.popover[role="tooltip"]')
    el = tooltip.find_element_by_class_name("popover-body")
    assert el.get_attribute("innerHTML") == "<strong>_Name</strong>: 3-methylpentane"

@flaky(max_runs=3, min_passes=1)
def test_tooltip_fmt(driver, grid):
    doc = get_doc(grid, {"tooltip": ["_Name"], "tooltip_fmt": "<em>{value}</em>"})
    driver.get(doc)
    driver.wait_for_img_load()
    (ActionChains(driver)
        .move_to_element(driver.find_by_css_selector("#mols2grid .cell .data-img *"))
        .perform())
    tooltip = driver.find_by_css_selector('div.popover[role="tooltip"]')
    el = tooltip.find_element_by_class_name("popover-body")
    assert el.get_attribute("innerHTML") == "<em>3-methylpentane</em>"

def test_tooltip_not_in_subset(driver, grid):
    doc = get_doc(grid, {"tooltip": ["_Name"], "subset": ["ID", "img"]})
    driver.get(doc)
    driver.wait_for_img_load()
    (ActionChains(driver)
        .move_to_element(driver.find_by_css_selector("#mols2grid .cell .data-img *"))
        .perform())
    tooltip = driver.find_by_css_selector('div.popover[role="tooltip"]')
    el = tooltip.find_element_by_class_name("popover-body")
    assert el.get_attribute("innerHTML") == "<strong>_Name</strong>: 3-methylpentane"

@flaky(max_runs=3, min_passes=1)
def test_style(driver, grid):
    doc = get_doc(grid, {
        "tooltip": ["_Name"],
        "style": {
            "__all__": lambda x: "color: red",
            "_Name": lambda x: "color: blue",
    }})
    driver.get(doc)
    driver.wait_for_img_load()
    el = driver.find_by_css_selector("#mols2grid .cell")
    assert el.value_of_css_property("color") == "rgb(255, 0, 0)"
    el = driver.find_by_css_selector("#mols2grid .cell .data-_Name")
    assert el.value_of_css_property("color") == "rgb(0, 0, 255)"
    (ActionChains(driver)
        .move_to_element(driver.find_by_css_selector("#mols2grid .cell .data-img *"))
        .perform())
    tooltip = driver.find_by_css_selector('div.popover[role="tooltip"]')
    el = tooltip.find_element_by_class_name("popover-body")
    assert el.get_attribute("innerHTML") == '<strong>_Name</strong>: <span style="color: blue">3-methylpentane</span>'

@flaky(max_runs=3, min_passes=1)
def test_transform(driver, grid):
    doc = get_doc(grid, {
        "tooltip": ["_Name"],
        "transform": {"_Name": lambda x: x.upper()}
    })
    driver.get(doc)
    name = driver.find_by_css_selector("#mols2grid .cell .data-_Name")
    assert name.text == "3-METHYLPENTANE"
    driver.wait_for_img_load()
    (ActionChains(driver)
        .move_to_element(driver.find_by_css_selector("#mols2grid .cell .data-img *"))
        .pause(.5)
        .perform())
    tooltip = driver.find_by_css_selector('div.popover[role="tooltip"]')
    el = tooltip.find_element_by_class_name("popover-body")
    assert el.get_attribute("innerHTML") == '<strong>_Name</strong>: 3-METHYLPENTANE'

@flaky(max_runs=3, min_passes=1)
def test_transform_style_tooltip(driver, grid):
    doc = get_doc(grid, {
        "tooltip": ["_Name"],
        "transform": {"_Name": lambda x: "foo"},
        "style": {
            "__all__": lambda x: "background-color: red",
            "_Name": lambda x: "color: green" if x == "foo" else "color: blue",
    }})
    driver.get(doc)
    driver.wait_for_img_load()
    cell = driver.find_by_css_selector("#mols2grid .cell")
    assert cell.value_of_css_property("background-color") == "rgb(255, 0, 0)"
    name = cell.find_element_by_class_name("data-_Name")
    assert name.text == "foo"   
    (ActionChains(driver)
        .move_to_element(driver.find_by_css_selector("#mols2grid .cell .data-img *"))
        .pause(.5)
        .perform())
    tooltip = driver.find_by_css_selector('div.popover[role="tooltip"]')
    el = tooltip.find_element_by_class_name("popover-body")
    assert el.get_attribute("innerHTML") == '<strong>_Name</strong>: <span style="color: blue">foo</span>'

@pytest.mark.parametrize("selection", [True, False])
def test_callback_js(driver, grid, selection):
    doc = get_doc(grid, {
        "subset": ["img", "_Name"],
        "callback": "$('#mols2grid .cell .data-_Name').html('foo')",
        "selection": selection
    })
    driver.get(doc)
    driver.wait_for_img_load()
    driver.find_clickable(By.CSS_SELECTOR, "#mols2grid .cell .data-img").click()
    el = driver.find_by_css_selector("#mols2grid .cell .data-_Name")
    assert el.text == "foo"

def test_sort_by(driver, grid):
    doc = get_doc(grid, {"subset": ["img", "_Name"], "sort_by": "_Name"})
    driver.get(doc)
    el = driver.find_by_css_selector("#mols2grid .cell .data-_Name")
    assert el.text == "1,1,2,2-tetrachloroethane"

def test_sort_button(driver, html_doc):
    driver.get(html_doc)
    driver.wait_for_img_load()
    driver.find_clickable(By.ID, "sortDropdown").click()
    driver.find_clickable(By.CSS_SELECTOR,
                          'button.sort-btn[data-name="data-_Name"]').click()
    el = driver.find_by_css_selector("#mols2grid .cell .data-_Name")
    assert el.text == "1,1,2,2-tetrachloroethane"
    driver.find_clickable(By.ID, "sortDropdown").click()
    driver.find_clickable(By.CSS_SELECTOR,
                          'button.sort-btn[data-name="data-_Name"]').click()
    el = driver.find_by_css_selector("#mols2grid .cell .data-_Name")
    assert el.text == "tetrachloromethane"

@pytest.mark.parametrize(["substruct_highlight", "expected"], [
    (True,  "2d943c1c2ba5e4bc75ed307d9e1dc456"),
    (False, "343b129e85fd56945a7bcfa7c12d63f7")
])
def test_substruct_highlight(driver, grid, substruct_highlight, expected):
    doc = get_doc(grid, {"n_rows": 1, "substruct_highlight": substruct_highlight})
    driver.get(doc)
    driver.wait_for_img_load()
    driver.substructure_query("CC(I)C")
    md5_hash = driver.get_svg_md5_hash()
    assert md5_hash == expected

def test_substruct_clear_removes_highlight(driver, grid):
    doc = get_doc(grid, {"n_rows": 1, "substruct_highlight": True})
    driver.get(doc)
    driver.wait_for_img_load()
    driver.substructure_query("C")
    md5_hash_hl = driver.get_svg_md5_hash()
    driver.find_by_id("searchbar").clear()
    driver.find_by_id("searchbar").send_keys(Keys.BACKSPACE)
    driver.wait_for_img_load()
    md5_hash = driver.get_svg_md5_hash()
    assert md5_hash != md5_hash_hl
    assert md5_hash == "d8227a9948b6f801193085e5a8b257a2"

def test_smarts_to_text_search_removes_highlight(driver, grid):
    doc = get_doc(grid, {"n_rows": 1, "substruct_highlight": True})
    driver.get(doc)
    driver.wait_for_img_load()
    driver.substructure_query("I")
    md5_hash_hl = driver.get_svg_md5_hash()
    driver.text_search("odopropane")
    md5_hash = driver.get_svg_md5_hash()
    assert md5_hash != md5_hash_hl
    assert md5_hash == "343b129e85fd56945a7bcfa7c12d63f7"

def test_filter(driver, grid):
    doc = get_doc(grid, {"subset": ["img", "_Name"],})
    driver.get(doc)
    mask = grid.dataframe["_Name"].str.contains("iodopropane")
    filter_code = env.get_template('js/filter.js').render(
        grid_id = grid._grid_id,
        mask = json.dumps(mask.tolist()))
    driver.wait_for_img_load()
    driver.execute_script(filter_code)
    el = driver.find_by_css_selector("#mols2grid .cell .data-_Name")
    assert el.text == "2-iodopropane"

def test_subset_gives_rows_order(driver, grid):
    subset = ["_Name", "img", "ID"]
    doc = get_doc(grid, {"subset": subset, "n_rows": 1})
    driver.get(doc)
    driver.wait_for_img_load()
    cell = driver.find_by_css_selector("#mols2grid .cell")
    elements = cell.find_elements_by_class_name("data")
    for i, el in enumerate(elements):
        class_list = el.get_attribute("class").split()
        name = [x.replace("data-", "") for x in class_list
                if x.startswith("data-")][0]
        if i == len(elements) - 1:
            assert name == "SMILES"
            assert el.value_of_css_property("display") == "none"
        else:
            assert name == subset[i]

def test_colname_with_spaces(driver, df):
    df = (df
        .rename(columns={"SMILES": "Molecule", "_Name": "Molecule name"})
        .drop(columns="mol"))
    grid = mols2grid.MolGrid(df, smiles_col="Molecule")
    doc = get_doc(grid, dict(
        subset=["Molecule name", "img"],
        tooltip=["Molecule"],
        n_rows=1))
    driver.get(doc)
    driver.wait_for_img_load()
    el = driver.find_by_css_selector("#mols2grid .cell .data")
    assert el.text == "3-methylpentane"

@flaky(max_runs=3, min_passes=1)
def test_custom_header(driver, grid):
    doc = get_doc(grid, {
        "subset": ["img"],
        "custom_header": '<script src="https://unpkg.com/@rdkit/rdkit@2021.9.2/Code/MinimalLib/dist/RDKit_minimal.js"></script>',
        "n_rows": 1})
    driver.get(doc)
    driver.wait_for_img_load()
    val = driver.execute_script("return RDKit.version();")
    assert val == "2021.09.2"

def test_table_template(driver):
    df = mols2grid.sdf_to_dataframe(sdf_path)[:15]
    grid = mols2grid.MolGrid(df, mol_col="mol", prerender=True)
    doc = get_doc(grid, dict(
        template="table", 
        subset=["mols2grid-id", "img"],
        tooltip=["_Name"],
        sort_by="_Name"
    ))
    driver.get(doc)
    el = driver.find_by_css_selector("#mols2grid td.col-0")
    assert el.find_element_by_class_name("data-mols2grid-id").text == "8"
    (ActionChains(driver)
        .move_to_element(el.find_element_by_css_selector(".data-img *"))
        .perform())
    tooltip = driver.find_by_css_selector('div.popover[role="tooltip"]')
    el = tooltip.find_element_by_class_name("popover-body")
    assert el.get_attribute("innerHTML") == "<strong>_Name</strong>: 1,3,5-trimethylbenzene"
    md5_hash = driver.get_svg_md5_hash("#mols2grid td .data-img")
    if rdkit_version == "2020.03.1":
        assert md5_hash == "7ca709f65c41fcfe090e98525920fb40"
    else:
        assert md5_hash == "9a681d0d978a6cb7ed28a95850d775b8"

def test_default_subset_tooltip(driver, grid):
    doc = get_doc(grid, {"n_rows": 1})
    driver.get(doc)
    driver.wait_for_img_load()
    expected_subset = ["mols2grid-id", "img"]
    expected_tooltip = [x for x in grid.dataframe.columns.drop("mol").to_list()
                        if x not in expected_subset]
    cell = driver.find_by_css_selector("#mols2grid .cell")
    data_elements = cell.find_elements_by_class_name("data")
    subset = [c.replace("data-", "").replace("-copy", "")
               for x in data_elements
               for c in x.get_attribute("class").split(" ")
               if c.startswith("data-") and not x.get_attribute("style")]
    assert subset == expected_subset
    tooltip = [c.replace("data-", "").replace("-copy", "")
               for x in data_elements
               for c in x.get_attribute("class").split(" ")
               if c.startswith("data-")
               and x.get_attribute("style") == "display: none;"]
    assert tooltip == expected_tooltip

def test_multi_highlight(driver, html_doc):
    driver.get(html_doc)
    driver.wait_for_img_load()
    driver.substructure_query("Br")
    md5_hash = driver.get_svg_md5_hash()
    assert md5_hash == "814a4c5dfa5226b629031861508505c5"

def test_single_highlight(driver, grid):
    doc = get_doc(grid, {"single_highlight": True})
    driver.get(doc)
    driver.wait_for_img_load()
    driver.substructure_query("Br")
    md5_hash = driver.get_svg_md5_hash()
    assert md5_hash == "ca043475667e4f79a3b87a96d7dd731a"

def test_mol_depiction_aligned_to_query(driver, html_doc):
    driver.get(html_doc)
    driver.wait_for_img_load()
    driver.substructure_query("CCCBr")
    images = driver.wait(EC.presence_of_all_elements_located(
                         (By.CSS_SELECTOR, "#mols2grid .cell .data-img")))
    for img, expected in zip(images, [
        "4b19aa7cf16fe691ea7840a3ce9601b5",
        "4d7a5d8cd6652c1fd8d04a6fffdd2e3a",
    ]):
        im = svg2png(bytestring=(img.get_attribute("innerHTML")))
        im = Image.open(BytesIO(im))
        md5_hash = md5(im.tobytes()).hexdigest()
        assert md5_hash == expected

def test_highlight_with_hydrogens(driver, df):
    mols = []
    for mol in df["mol"][25:]:
        m = Chem.AddHs(mol)
        m.ClearProp("SMILES")
        mols.append(m)
    grid = mols2grid.MolGrid.from_mols(mols, removeHs=False)
    doc = get_doc(grid, dict(substruct_highlight=True, single_highlight=False))
    driver.get(doc)
    driver.wait_for_img_load()
    driver.substructure_query("Cl")
    md5_hash = driver.get_svg_md5_hash()
    assert md5_hash == "1133a4662a00bde52798fa5673251d8c"
