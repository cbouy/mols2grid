import pytest
from functools import partial
from base64 import b64encode
from selenium import webdriver
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions
from selenium.webdriver.common.by import By
import geckodriver_autoinstaller
from rdkit import RDConfig
import mols2grid
from mols2grid.dispatch import _SIGNATURE

geckodriver_autoinstaller.install()

DEFAULTS = {name: _SIGNATURE[meth][name].default
            for meth in ["render", "to_pages", "display"]
            for name in _SIGNATURE[meth]}

get_grid = partial(mols2grid.MolGrid.from_sdf, 
                   f"{RDConfig.RDDocsDir}/Book/data/solubility.test.sdf")

def get_doc(grid, kwargs):
    html = grid.render(**kwargs)
    html = b64encode(html.encode()).decode()
    return "data:text/html;base64,{}".format(html)

@pytest.fixture(scope="module")
def grid():
    return get_grid()

@pytest.fixture(scope="module")
def html_doc(grid):
    html = grid.render()
    html = b64encode(html.encode()).decode()
    return "data:text/html;base64,{}".format(html)

@pytest.fixture(scope="module")
def driver():
    options = webdriver.FirefoxOptions()
    options.headless = True
    d = webdriver.Firefox(options=options)
    yield d
    d.quit()

@pytest.mark.parametrize("page", [1, 2, 3])
@pytest.mark.parametrize("n_cols", [DEFAULTS["n_cols"], 6])
@pytest.mark.parametrize("n_rows", [DEFAULTS["n_rows"], 6])
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

def test_style(driver, grid):
    doc = get_doc(grid, {"style": {
        "__all__": lambda x: "color: red",
        "_Name": lambda x: "color: blue",
        }})
    driver.get(doc)
    wait = WebDriverWait(driver, 8)
    el = wait.until(
        expected_conditions.presence_of_element_located((By.CSS_SELECTOR, "#mols2grid .cell .data-img svg"))
    )
    el = driver.find_element_by_css_selector("#mols2grid .cell")
    assert el.value_of_css_property("color") == "rgb(255, 0, 0)"
    el = driver.find_element_by_css_selector("#mols2grid .cell .data-_Name")
    assert el.value_of_css_property("color") == "rgb(0, 0, 255)"

def test_no_subset_all_visible(html_doc):
    pass

def test_m2g_id_in_subset(driver, grid):
    pass

def test_smiles_hidden(driver, grid):
    pass

def test_text_search(html_doc):
    pass

def test_smarts_search(html_doc):
    pass

def test_selection_click(driver, grid):
    pass

def test_selection_with_cache(driver, grid):
    pass

def test_selection_check_uncheck_all(html_doc):
    pass

def test_selection_invert(html_doc):
    pass

def test_selection_clipboard(html_doc):
    pass

def test_image_args(driver):
    pass

def test_tooltip(driver, grid):
    pass

def test_hover(driver, grid):
    pass

def test_transform(driver, grid):
    pass

def test_callback_js(driver, grid):
    pass

def test_sort_by(driver, grid):
    pass

def test_substruct_highlight(driver, grid):
    pass
