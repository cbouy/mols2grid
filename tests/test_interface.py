import pytest
from flaky import flaky
from functools import partial
from base64 import b64encode, b64decode
from selenium import webdriver
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions
from selenium.webdriver.common.by import By
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.common.action_chains import ActionChains
from selenium.common.exceptions import NoSuchElementException
import geckodriver_autoinstaller
from rdkit import Chem, RDConfig
from rdkit.Chem import AllChem
from PIL import Image
from io import BytesIO
from hashlib import md5
from cairosvg import svg2png
import mols2grid
from mols2grid.select import register

geckodriver_autoinstaller.install()

sdf_path = f"{RDConfig.RDDocsDir}/Book/data/solubility.test.sdf"

get_grid = partial(mols2grid.MolGrid.from_sdf, sdf_path)

@pytest.fixture(scope="module")
def driver():
    options = webdriver.FirefoxOptions()
    options.headless = True
    d = webdriver.Firefox(options=options)
    yield d
    d.quit()

@pytest.fixture(scope="module")
def df():
    return mols2grid.sdf_to_dataframe(sdf_path)[:30]

@pytest.fixture(scope="module")
def mols(df):
    return df["mol"][:5]

@pytest.fixture(scope="module")
def grid(df):
    return mols2grid.MolGrid(df, mol_col="mol")

@pytest.fixture(scope="module")
def html_doc(grid):
    html = grid.render()
    html = b64encode(html.encode()).decode()
    return "data:text/html;base64,{}".format(html)

def get_doc(grid, kwargs):
    html = grid.render(**kwargs)
    html = b64encode(html.encode()).decode()
    return "data:text/html;base64,{}".format(html)

def wait_for_img_load(driver):
    WebDriverWait(driver, 8).until(
        expected_conditions.presence_of_element_located((By.CSS_SELECTOR,
            "#mols2grid .cell .data-img svg")))

# make sure non-parametrized test is ran first
@pytest.mark.order(1)
def test_no_subset_all_visible(driver, grid):
    doc = get_doc(grid, {"selection": False})
    driver.get(doc)
    columns = set(grid.dataframe.columns.drop("mol").to_list())
    cell = driver.find_element_by_css_selector("#mols2grid .cell")
    data_el = cell.find_elements_by_class_name("data")
    classes = [c.replace("data-", "").replace("-copy", "")
               for x in data_el
               for c in x.get_attribute("class").split(" ")
               if c.startswith("data-")]
    classes = set(classes)
    assert classes == columns

def test_smiles_hidden(driver, grid):
    doc = get_doc(grid, {"subset": ["_Name", "img"]})
    driver.get(doc)
    el = driver.find_element_by_css_selector("#mols2grid .cell .data-SMILES")
    assert not el.is_displayed()

@pytest.mark.parametrize("page", [1, 2, 3])
@pytest.mark.parametrize("n_cols", [1, 3])
@pytest.mark.parametrize("n_rows", [1, 3])
def test_page_click(driver, grid, page, n_cols, n_rows):
    doc = get_doc(grid, dict(n_cols=n_cols, n_rows=n_rows))
    driver.get(doc)
    for i in range(2, page+1):
        next_page = driver.find_element_by_css_selector(f'a.page-link[data-i="{i}"]')
        next_page.click()
    first_cell = driver.find_element_by_class_name('cell')
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

def test_text_search(driver, grid):
    doc = get_doc(grid, {"n_rows": 1})
    driver.get(doc)
    text_box = driver.find_element_by_css_selector("#searchbar")
    wait_for_img_load(driver)
    (ActionChains(driver)
        .send_keys_to_element(text_box, "iodopropane")
        .key_up(Keys.TAB)
        .perform())
    el = driver.find_element_by_css_selector("#mols2grid .cell .data-SMILES")
    assert el.text == "CC(I)C"

@flaky(max_runs=3, min_passes=1)
def test_smarts_search(driver, grid):
    doc = get_doc(grid, {"subset": ["_Name", "img"]})
    driver.get(doc)
    text_box = driver.find_element_by_id("searchbar")
    wait_for_img_load(driver)
    (ActionChains(driver)
        .click(driver.find_element_by_id("searchBtn"))
        .pause(.5)
        .click(driver.find_element_by_id("smartsSearch"))
        .send_keys_to_element(text_box, "CC(I)C")
        .key_up(Keys.TAB)
        .perform())
    el = driver.find_element_by_css_selector("#mols2grid .cell .data-_Name")
    assert el.text == "2-iodopropane"

def test_selection_click(driver, html_doc):
    driver.get(html_doc)
    chkbox = driver.find_element_by_css_selector("input[type='checkbox']")
    wait_for_img_load(driver)
    (ActionChains(driver)
        .click(chkbox)
        .pause(.5)
        .perform())
    sel = driver.execute_script("return SELECTION.to_dict();")
    sel = eval(sel)
    assert sel == {0: "CCC(C)CC"}

def test_selection_with_cache_check_and_uncheck(driver):
    register._init_grid("cached_sel")
    register._set_selection(0, "CCC(C)CC")
    grid = get_grid(name="cached_sel", cache_selection=True)
    doc = get_doc(grid, {})
    driver.get(doc)
    wait_for_img_load(driver)
    sel = driver.execute_script("return SELECTION.to_dict();")
    sel = eval(sel)
    assert sel == {0: "CCC(C)CC"}
    (ActionChains(driver)
        .click(driver.find_element_by_css_selector("input[type='checkbox']"))
        .pause(.5)
        .perform())
    sel = driver.execute_script("return SELECTION.to_dict();")
    sel = eval(sel)
    assert sel == {}
    register._clear()

def test_selection_check_uncheck_invert(driver, grid):
    doc = get_doc(grid, {"subset": ["_Name", "img"]})
    driver.get(doc)
    wait_for_img_load(driver)
    (ActionChains(driver)
        .click(driver.find_element_by_id("chkboxDropdown"))
        .pause(.5)
        .click(driver.find_element_by_id("btn-chkbox-all"))
        .perform())
    sel = driver.execute_script("return SELECTION.to_dict();")
    sel = eval(sel)
    assert len(sel) == 30
    (ActionChains(driver)
        .click(driver.find_element_by_id("chkboxDropdown"))
        .pause(.5)
        .click(driver.find_element_by_id("btn-chkbox-none"))
        .perform())
    sel = driver.execute_script("return SELECTION.to_dict();")
    sel = eval(sel)
    assert sel == {}
    (ActionChains(driver)
        .click(driver.find_element_by_css_selector("input[type='checkbox']"))
        .click(driver.find_element_by_id("chkboxDropdown"))
        .pause(.5)
        .click(driver.find_element_by_id("btn-chkbox-invert"))
        .perform())
    sel = driver.execute_script("return SELECTION.to_dict();")
    sel = eval(sel)
    assert len(sel) == 29

@pytest.mark.parametrize("prerender", [True, False])
def test_image_size(driver, prerender):
    grid = get_grid(size=(100, 100), prerender=prerender)
    doc = get_doc(grid, {"selection": False, "border": "0", "substruct_highlight": False})
    driver.get(doc)
    if not prerender:
        wait_for_img_load(driver)
    img = driver.find_element_by_css_selector("#mols2grid .cell .data-img *")
    assert img.size == {"height": 100.0, "width": 100.0}

def test_image_use_coords(driver, df):
    mols = df["mol"][:1]
    AllChem.EmbedMolecule(mols[0], randomSeed=0xf00d)
    grid = mols2grid.MolGrid.from_mols(
        mols, use_coords=True, prerender=True, useSVG=False)
    doc = get_doc(grid, {"substruct_highlight": False})
    driver.get(doc)
    img = driver.find_element_by_css_selector("#mols2grid .cell .data-img *")
    im = Image.open(BytesIO(b64decode(img.get_attribute("src")[22:])))
    md5_hash = md5(im.tobytes()).hexdigest()
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
        wait_for_img_load(driver)
    img = driver.find_element_by_css_selector("#mols2grid .cell .data-img *")
    if useSVG:
        im = svg2png(bytestring=(img
                                 .find_element_by_xpath("..")
                                 .get_attribute("innerHTML")))
        im = Image.open(BytesIO(im))
    else:
        im = Image.open(BytesIO(b64decode(img.get_attribute("src")[22:])))
    md5_hash = md5(im.tobytes()).hexdigest()
    assert md5_hash == expected

@pytest.mark.parametrize(["removeHs", "prerender", "expected"], [
    (True, True, "acf7cf7cd5cfaa5eaf4a4f257e290e49"),
    (False, True, "b01c2146523b8a006bb34621839668eb"),
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
        wait_for_img_load(driver)
    img = driver.find_element_by_css_selector("#mols2grid .cell .data-img *")
    if useSVG:
        im = svg2png(bytestring=(img
                                 .find_element_by_xpath("..")
                                 .get_attribute("innerHTML")))
        im = Image.open(BytesIO(im))
    else:
        im = Image.open(BytesIO(b64decode(img.get_attribute("src")[22:])))
    md5_hash = md5(im.tobytes()).hexdigest()
    assert md5_hash == expected

def test_hover_color(driver, grid):
    doc = get_doc(grid, {"hover_color": "red"})
    driver.get(doc)
    cell = driver.find_element_by_css_selector("#mols2grid .cell")
    (ActionChains(driver)
     .move_to_element(cell)
     .perform())
    color = driver.execute_script(f"return getComputedStyle(document.querySelector('#mols2grid .cell')).getPropertyValue('background-color');")
    assert color == "rgb(255, 0, 0)"

@flaky(max_runs=3, min_passes=1)
def test_tooltip(driver, grid):
    doc = get_doc(grid, {"tooltip": ["_Name"]})
    driver.get(doc)
    wait_for_img_load(driver)
    img = driver.find_element_by_css_selector("#mols2grid .cell .data-img *")
    (ActionChains(driver)
     .move_to_element(img)
     .perform())
    tooltip = driver.find_element_by_css_selector('div.popover[role="tooltip"]')
    el = tooltip.find_element_by_class_name("popover-body")
    assert el.get_attribute("innerHTML") == "<strong>_Name</strong>: 3-methylpentane"

@flaky(max_runs=3, min_passes=1)
def test_tooltip_trigger(driver, grid):
    doc = get_doc(grid, {"tooltip": ["_Name"], "tooltip_trigger": "click"})
    driver.get(doc)
    wait_for_img_load(driver)
    img = driver.find_element_by_css_selector("#mols2grid .cell .data-img *")
    (ActionChains(driver)
     .move_to_element(img)
     .perform())
    with pytest.raises(NoSuchElementException):
        driver.find_element_by_css_selector('div.popover[role="tooltip"]')
    (ActionChains(driver)
     .click(img)
     .perform())
    tooltip = driver.find_element_by_css_selector('div.popover[role="tooltip"]')
    el = tooltip.find_element_by_class_name("popover-body")
    assert el.get_attribute("innerHTML") == "<strong>_Name</strong>: 3-methylpentane"

@flaky(max_runs=3, min_passes=1)
def test_tooltip_fmt(driver, grid):
    doc = get_doc(grid, {"tooltip": ["_Name"], "tooltip_fmt": "<em>{value}</em>"})
    driver.get(doc)
    wait_for_img_load(driver)
    img = driver.find_element_by_css_selector("#mols2grid .cell .data-img *")
    (ActionChains(driver)
     .move_to_element(img)
     .perform())
    tooltip = driver.find_element_by_css_selector('div.popover[role="tooltip"]')
    el = tooltip.find_element_by_class_name("popover-body")
    assert el.get_attribute("innerHTML") == "<em>3-methylpentane</em>"

@flaky(max_runs=3, min_passes=1)
def test_style(driver, grid):
    doc = get_doc(grid, {
        "tooltip": ["_Name"],
        "style": {
            "__all__": lambda x: "color: red",
            "_Name": lambda x: "color: blue",
    }})
    driver.get(doc)
    wait_for_img_load(driver)
    el = driver.find_element_by_css_selector("#mols2grid .cell")
    assert el.value_of_css_property("color") == "rgb(255, 0, 0)"
    el = driver.find_element_by_css_selector("#mols2grid .cell .data-_Name")
    assert el.value_of_css_property("color") == "rgb(0, 0, 255)"
    img = driver.find_element_by_css_selector("#mols2grid .cell .data-img *")
    (ActionChains(driver)
     .move_to_element(img)
     .perform())
    tooltip = driver.find_element_by_css_selector('div.popover[role="tooltip"]')
    el = tooltip.find_element_by_class_name("popover-body")
    assert el.get_attribute("innerHTML") == '<strong>_Name</strong>: <span style="color: blue">3-methylpentane</span>'

@flaky(max_runs=3, min_passes=1)
def test_transform(driver, grid):
    doc = get_doc(grid, {
        "tooltip": ["_Name"],
        "transform": {"_Name": lambda x: x.upper()}
    })
    driver.get(doc)
    name = driver.find_element_by_css_selector("#mols2grid .cell .data-_Name")
    assert name.text == "3-METHYLPENTANE"
    wait_for_img_load(driver)
    img = driver.find_element_by_css_selector("#mols2grid .cell .data-img *")
    (ActionChains(driver)
     .move_to_element(img)
     .pause(.5)
     .perform())
    tooltip = driver.find_element_by_css_selector('div.popover[role="tooltip"]')
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
    wait_for_img_load(driver)
    cell = driver.find_element_by_css_selector("#mols2grid .cell")
    assert cell.value_of_css_property("background-color") == "rgb(255, 0, 0)"
    name = cell.find_element_by_class_name("data-_Name")
    assert name.text == "foo"   
    img = driver.find_element_by_css_selector("#mols2grid .cell .data-img *")
    (ActionChains(driver)
     .move_to_element(img)
     .pause(.5)
     .perform())
    tooltip = driver.find_element_by_css_selector('div.popover[role="tooltip"]')
    el = tooltip.find_element_by_class_name("popover-body")
    assert el.get_attribute("innerHTML") == '<strong>_Name</strong>: <span style="color: blue">foo</span>'

def test_callback_js(driver, grid):
    doc = get_doc(grid, {
        "callback": "$('#mols2grid .cell .data-_Name').html('foo')"
    })
    driver.get(doc)
    wait_for_img_load(driver)
    (ActionChains(driver)
     .click(driver.find_element_by_css_selector("#mols2grid .cell .data-img"))
     .pause(.5)
     .perform())
    el = driver.find_element_by_css_selector("#mols2grid .cell .data-_Name")
    assert el.text == "foo"

def test_sort_by(driver, grid):
    doc = get_doc(grid, {"sort_by": "_Name"})
    driver.get(doc)
    el = driver.find_element_by_css_selector("#mols2grid .cell .data-_Name")
    assert el.text == "1,1,2,2-tetrachloroethane"

@pytest.mark.parametrize(["substruct_highlight", "expected"], [
    (True,  "6a1ae3fc97cbe003c86db3c624e5ff09"),
    (False, "4884c75977e8c3e0fa2e6a1607275406")
])
def test_substruct_highlight(driver, grid, substruct_highlight, expected):
    doc = get_doc(grid, {"n_rows": 1, "substruct_highlight": substruct_highlight})
    driver.get(doc)
    text_box = driver.find_element_by_id("searchbar")
    wait_for_img_load(driver)
    (ActionChains(driver)
        .click(driver.find_element_by_id("searchBtn"))
        .pause(.5)
        .click(driver.find_element_by_id("smartsSearch"))
        .send_keys_to_element(text_box, "CC(I)C")
        .key_up(Keys.TAB)
        .perform())
    img = driver.find_element_by_css_selector("#mols2grid .cell .data-img")
    im = svg2png(bytestring=(img.get_attribute("innerHTML")))
    im = Image.open(BytesIO(im))
    md5_hash = md5(im.tobytes()).hexdigest()
    assert md5_hash == expected

def test_substruct_clear_removes_highlight(driver, grid):
    doc = get_doc(grid, {"n_rows": 1, "substruct_highlight": True})
    driver.get(doc)
    text_box = driver.find_element_by_id("searchbar")
    wait_for_img_load(driver)
    (ActionChains(driver)
        .click(driver.find_element_by_id("searchBtn"))
        .pause(.5)
        .click(driver.find_element_by_id("smartsSearch"))
        .send_keys_to_element(text_box, "C")
        .key_up(Keys.TAB)
        .perform())
    img = driver.find_element_by_css_selector("#mols2grid .cell .data-img")
    im = svg2png(bytestring=(img.get_attribute("innerHTML")))
    im = Image.open(BytesIO(im))
    md5_hash_hl = md5(im.tobytes()).hexdigest()
    text_box.clear()
    (ActionChains(driver)
        .send_keys_to_element(text_box, Keys.BACKSPACE)
        .key_up(Keys.TAB)
        .perform())
    img = driver.find_element_by_css_selector("#mols2grid .cell .data-img")
    im = svg2png(bytestring=(img.get_attribute("innerHTML")))
    im = Image.open(BytesIO(im))
    md5_hash = md5(im.tobytes()).hexdigest()
    assert md5_hash != md5_hash_hl
    assert md5_hash == "d8227a9948b6f801193085e5a8b257a2"
