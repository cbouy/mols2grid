import pytest
from flaky import flaky
from base64 import b64encode, b64decode
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

class FirefoxDriver(webdriver.Firefox):
    def wait_for_img_load(self, max_delay=15, selector="#mols2grid .cell .data-img svg"):
        return (WebDriverWait(self, max_delay)
                .until(EC.presence_of_element_located(
                        (By.CSS_SELECTOR, selector))))

    def wait(self, condition, max_delay=5):
        return (WebDriverWait(self, max_delay,
                      ignored_exceptions=[StaleElementReferenceException])
                .until(condition))

    def find_by_id(self, element_id, **kwargs):
        condition = EC.presence_of_element_located((By.ID, element_id))
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


@pytest.fixture(scope="function")
def driver():
    options = webdriver.FirefoxOptions()
    options.headless = True
    driver = FirefoxDriver(options=options)
    driver.set_page_load_timeout(10)
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
    doc = get_doc(grid, {"selection": False})
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
@pytest.mark.parametrize("n_cols", [1, 3])
@pytest.mark.parametrize("n_rows", [1, 3])
def test_page_click(driver, grid, page, n_cols, n_rows):
    doc = get_doc(grid, dict(n_cols=n_cols, n_rows=n_rows))
    driver.get(doc)
    for i in range(2, page+1):
        next_page = driver.find_by_css_selector(f'a.page-link[data-i="{i}"]')
        next_page.click()
    first_cell = driver.find_by_class_name('cell')
    mols2grid_id = n_cols * n_rows * (page - 1)
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
    (ActionChains(driver)
        .send_keys_to_element(driver.find_by_id("searchbar"), "iodopropane")
        .key_up(Keys.TAB)
        .pause(.5)
        .perform())
    el = driver.find_by_css_selector("#mols2grid .cell .data-SMILES")
    assert el.get_attribute("innerHTML") == "CC(I)C"

@flaky(max_runs=3, min_passes=1)
def test_smarts_search(driver, html_doc):
    driver.get(html_doc)
    driver.wait_for_img_load()
    (ActionChains(driver)
        .click(driver.find_by_id("searchBtn"))
        .pause(.5)
        .click(driver.find_by_id("smartsSearch"))
        .send_keys_to_element(driver.find_by_id("searchbar"), "CC(I)C")
        .key_up(Keys.TAB)
        .pause(.5)
        .perform())
    el = driver.find_by_css_selector("#mols2grid .cell .data-_Name")
    assert el.text == "2-iodopropane"

def test_selection_click(driver, html_doc):
    driver.get(html_doc)
    driver.wait_for_img_load()
    (ActionChains(driver)
        .click(driver.find_by_css_selector("input[type='checkbox']"))
        .pause(.5)
        .perform())
    sel = driver.execute_script("return SELECTION.to_dict();")
    sel = eval(sel)
    assert sel == {0: "CCC(C)CC"}
    register._clear()

def test_selection_with_cache_check_and_uncheck(driver, df):
    register._init_grid("cached_sel")
    register._set_selection(0, "CCC(C)CC")
    grid = get_grid(df, name="cached_sel", cache_selection=True)
    doc = get_doc(grid, {})
    driver.get(doc)
    driver.wait_for_img_load()
    sel = driver.execute_script("return SELECTION.to_dict();")
    sel = eval(sel)
    assert sel == {0: "CCC(C)CC"}
    (ActionChains(driver)
        .click(driver.find_by_css_selector("input[type='checkbox']"))
        .pause(.5)
        .perform())
    sel = driver.execute_script("return SELECTION.to_dict();")
    sel = eval(sel)
    assert sel == {}
    register._clear()

def test_selection_check_uncheck_invert(driver, html_doc):
    driver.get(html_doc)
    driver.wait_for_img_load()
    (ActionChains(driver)
        .click(driver.find_by_id("chkboxDropdown"))
        .pause(.5)
        .click(driver.find_by_id("btn-chkbox-all"))
        .perform())
    sel = driver.execute_script("return SELECTION.to_dict();")
    sel = eval(sel)
    assert len(sel) == 30
    (ActionChains(driver)
        .click(driver.find_by_id("chkboxDropdown"))
        .pause(.5)
        .click(driver.find_by_id("btn-chkbox-none"))
        .perform())
    sel = driver.execute_script("return SELECTION.to_dict();")
    sel = eval(sel)
    assert sel == {}
    (ActionChains(driver)
        .click(driver.find_by_css_selector("input[type='checkbox']"))
        .click(driver.find_by_id("chkboxDropdown"))
        .pause(.5)
        .click(driver.find_by_id("btn-chkbox-invert"))
        .perform())
    sel = driver.execute_script("return SELECTION.to_dict();")
    sel = eval(sel)
    assert len(sel) == 29
    register._clear()

def test_check_all_selects_only_not_hidden(driver, html_doc):
    driver.get(html_doc)
    driver.wait_for_img_load()
    (ActionChains(driver)
        .send_keys_to_element(driver.find_by_id("searchbar"), "iodopropane")
        .key_up(Keys.TAB)
        .click(driver.find_by_id("chkboxDropdown"))
        .pause(.3)
        .click(driver.find_by_id("btn-chkbox-all"))
        .perform())
    driver.find_by_id("searchbar").clear()
    (ActionChains(driver)
        .send_keys_to_element(driver.find_by_id("searchbar"), Keys.BACKSPACE)
        .key_up(Keys.TAB)
        .perform())
    sel = driver.execute_script("return SELECTION.to_dict();")
    sel = eval(sel)
    assert sel == {27: 'CC(I)C'}
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
    img = driver.find_by_css_selector("#mols2grid .cell .data-img *")
    im = Image.open(BytesIO(b64decode(img.get_attribute("src")[22:])))
    md5_hash = md5(im.tobytes()).hexdigest()
    if rdkit_version == "2020.03.1":
        assert md5_hash == "aed60ed28347831d24f02dbb5be19007"
    else:
        assert md5_hash == "6909ba43f86003cea9b0fd8d723cddfe"

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
        img = driver.find_by_css_selector("#mols2grid .cell .data-img")
        img_data = img.get_attribute("innerHTML")
        im = svg2png(bytestring=img_data)
        im = Image.open(BytesIO(im))
    else:
        img = driver.find_by_css_selector("#mols2grid .cell .data-img *")
        img_data = img.get_attribute("src")[22:]
        im = Image.open(BytesIO(b64decode(img_data)))
    md5_hash = md5(im.tobytes()).hexdigest()
    assert md5_hash == expected

@pytest.mark.parametrize(["removeHs", "prerender", "expected"], [
    (True, True, "acf7cf7cd5cfaa5eaf4a4f257e290e49"),
    (False, True, ("2d30379732d312409d1d93ad5b3c3e60"
                   if rdkit_version == "2020.03.1"
                   else "b01c2146523b8a006bb34621839668eb")),
    (True, False, "d8227a9948b6f801193085e5a8b257a2"),
    (False, False, "e7b837b4d0c5f1ed8e1cd79734b6845e"),
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
        img = driver.find_by_css_selector("#mols2grid .cell .data-img")
        img_data = img.get_attribute("innerHTML")
        im = svg2png(bytestring=img_data)
        im = Image.open(BytesIO(im))
    else:
        img = driver.find_by_css_selector("#mols2grid .cell .data-img *")
        img_data = img.get_attribute("src")[22:]
        im = Image.open(BytesIO(b64decode(img_data)))
    md5_hash = md5(im.tobytes()).hexdigest()
    assert md5_hash == expected

@pytest.mark.parametrize(["kwargs", "expected"], [
    (dict(addAtomIndices=True), "d9611dd45ff82ab73cf31dfc3a33706d"),
    (dict(fixedBondLength=10), "1dfa2a5c43022ce66a2ff26a517698cf"),
    (dict(atomColourPalette={6: (0, .8, .8)}), "91d114d421fc6a05e4537d19ad4b8324"),
    (dict(legend="foo"), "bb5bc6ae713cfaea4e6374fc35e4be4e"),
])
def test_moldrawoptions(driver, df, kwargs, expected):
    grid = get_grid(df, **kwargs)
    doc = get_doc(grid, dict(n_rows=1, n_cols=1, subset=["img"]))
    driver.get(doc)
    driver.wait_for_img_load()
    img = driver.find_by_css_selector("#mols2grid .cell .data-img")
    im = svg2png(bytestring=(img.get_attribute("innerHTML")))
    im = Image.open(BytesIO(im))
    md5_hash = md5(im.tobytes()).hexdigest()
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
    (ActionChains(driver)
     .click(driver.find_by_css_selector("#mols2grid .cell .data-img *"))
     .perform())
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

def test_callback_js(driver, grid):
    doc = get_doc(grid, {
        "callback": "$('#mols2grid .cell .data-_Name').html('foo')"
    })
    driver.get(doc)
    driver.wait_for_img_load()
    (ActionChains(driver)
     .click(driver.find_by_css_selector("#mols2grid .cell .data-img"))
     .pause(.5)
     .perform())
    el = driver.find_by_css_selector("#mols2grid .cell .data-_Name")
    assert el.text == "foo"

def test_sort_by(driver, grid):
    doc = get_doc(grid, {"sort_by": "_Name"})
    driver.get(doc)
    el = driver.find_by_css_selector("#mols2grid .cell .data-_Name")
    assert el.text == "1,1,2,2-tetrachloroethane"

def test_sort_button(driver, html_doc):
    driver.get(html_doc)
    driver.wait_for_img_load()
    (ActionChains(driver)
        .click(driver.find_by_id("sortDropdown"))
        .pause(.3)
        .click(driver.find_by_css_selector('button.sort-btn[data-name="data-_Name"]'))
        .perform())
    el = driver.find_by_css_selector("#mols2grid .cell .data-_Name")
    assert el.text == "1,1,2,2-tetrachloroethane"
    (ActionChains(driver)
        .click(driver.find_by_id("sortDropdown"))
        .pause(.3)
        .click(driver.find_by_css_selector('button.sort-btn[data-name="data-_Name"]'))
        .perform())
    el = driver.find_by_css_selector("#mols2grid .cell .data-_Name")
    assert el.text == "tetrachloromethane"

@pytest.mark.parametrize(["substruct_highlight", "expected"], [
    (True,  "6a1ae3fc97cbe003c86db3c624e5ff09"),
    (False, "4884c75977e8c3e0fa2e6a1607275406")
])
def test_substruct_highlight(driver, grid, substruct_highlight, expected):
    doc = get_doc(grid, {"n_rows": 1, "substruct_highlight": substruct_highlight})
    driver.get(doc)
    text_box = driver.find_by_id("searchbar")
    driver.wait_for_img_load()
    (ActionChains(driver)
        .click(driver.find_by_id("searchBtn"))
        .pause(.5)
        .click(driver.find_by_id("smartsSearch"))
        .send_keys_to_element(text_box, "CC(I)C")
        .key_up(Keys.TAB)
        .perform())
    img = driver.find_by_css_selector("#mols2grid .cell .data-img")
    im = svg2png(bytestring=(img.get_attribute("innerHTML")))
    im = Image.open(BytesIO(im))
    md5_hash = md5(im.tobytes()).hexdigest()
    assert md5_hash == expected

def test_substruct_clear_removes_highlight(driver, grid):
    doc = get_doc(grid, {"n_rows": 1, "substruct_highlight": True})
    driver.get(doc)
    driver.wait_for_img_load()
    (ActionChains(driver)
        .click(driver.find_by_id("searchBtn"))
        .pause(.5)
        .click(driver.find_by_id("smartsSearch"))
        .send_keys_to_element(driver.find_by_id("searchbar"), "C")
        .key_up(Keys.TAB)
        .perform())
    img = driver.find_by_css_selector("#mols2grid .cell .data-img")
    im = svg2png(bytestring=(img.get_attribute("innerHTML")))
    im = Image.open(BytesIO(im))
    md5_hash_hl = md5(im.tobytes()).hexdigest()
    driver.find_by_id("searchbar").clear()
    (ActionChains(driver)
        .send_keys_to_element(driver.find_by_id("searchbar"), Keys.BACKSPACE)
        .key_up(Keys.TAB)
        .perform())
    img = driver.find_by_css_selector("#mols2grid .cell .data-img")
    im = svg2png(bytestring=(img.get_attribute("innerHTML")))
    im = Image.open(BytesIO(im))
    md5_hash = md5(im.tobytes()).hexdigest()
    assert md5_hash != md5_hash_hl
    assert md5_hash == "d8227a9948b6f801193085e5a8b257a2"

def test_filter(driver, grid):
    doc = get_doc(grid, {})
    driver.get(doc)
    mask = grid.dataframe["_Name"].str.contains("iodopropane")
    ids = grid.dataframe.loc[mask]["mols2grid-id"].to_list()
    filter_code = env.get_template('js/filter.js').render(
        grid_id = grid._grid_id,
        ids = ids)
    driver.wait_for_img_load()
    driver.execute_script(filter_code)
    el = driver.find_by_css_selector("#mols2grid .cell .data-_Name")
    assert el.text == "2-iodopropane"
