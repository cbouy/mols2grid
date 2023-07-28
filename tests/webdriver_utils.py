from ast import literal_eval
from base64 import b64decode
from io import BytesIO

import imagehash
from cairosvg import svg2png
from PIL import Image
from selenium import webdriver
from selenium.common.exceptions import StaleElementReferenceException
from selenium.webdriver.common.action_chains import ActionChains
from selenium.webdriver.common.by import By
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.support.ui import WebDriverWait


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
    def wait_for_img_load(
        self, max_delay=5, selector="#mols2grid .m2g-cell .data-img svg"
    ):
        return WebDriverWait(self, max_delay).until(
            EC.presence_of_element_located((By.CSS_SELECTOR, selector))
        )

    def wait_for_selection(self, is_empty=False, max_delay=3):
        return WebDriverWait(self, max_delay).until(selection_available(is_empty))

    def wait(self, condition, max_delay=5):
        return WebDriverWait(
            self, max_delay, ignored_exceptions=[StaleElementReferenceException]
        ).until(condition)

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

    def get_svg_hash(self, *args, **kwargs):
        im = next(self.get_imgs_from_svgs(*args, **kwargs))
        return imagehash.average_hash(im, hash_size=16)

    def get_imgs_from_svgs(self, selector="#mols2grid .m2g-cell .data-img"):
        condition = EC.presence_of_all_elements_located((By.CSS_SELECTOR, selector))
        svgs = self.wait(condition)
        for svg in svgs:
            im = svg2png(bytestring=(svg.get_attribute("innerHTML")))
            yield Image.open(BytesIO(im))

    def get_png_hash(self, selector="#mols2grid .m2g-cell .data-img *"):
        img = self.find_by_css_selector(selector)
        im = Image.open(BytesIO(b64decode(img.get_attribute("src")[22:])))
        return imagehash.average_hash(im, hash_size=16)

    def substructure_query(self, smarts):
        self.find_clickable(By.CSS_SELECTOR, "#mols2grid .m2g-search-smarts").click()
        self.find_by_css_selector("#mols2grid .m2g-searchbar").send_keys(smarts)
        self.wait_for_img_load()

    def text_search(self, txt):
        self.find_clickable(By.CSS_SELECTOR, "#mols2grid .m2g-search-text").click()
        self.find_by_css_selector("#mols2grid .m2g-searchbar").send_keys(txt)
        self.wait_for_img_load()

    def clear_search(self):
        self.find_by_css_selector("#mols2grid .m2g-searchbar").clear()
        self.find_by_css_selector("#mols2grid .m2g-searchbar").send_keys(Keys.BACKSPACE)
        self.wait_for_img_load()

    def sort_grid(self, field):
        self.find_clickable(By.CSS_SELECTOR, "#mols2grid .m2g-sort").click()
        self.find_clickable(
            By.CSS_SELECTOR, f'#mols2grid .m2g-sort option[value="data-{field}"]'
        ).click()
        self.wait_for_img_load()

    def invert_sort(self):
        self.find_clickable(By.CSS_SELECTOR, "#mols2grid .m2g-sort .m2g-order").click()
        self.wait_for_img_load()

    def grid_action(self, action):
        self.find_clickable(By.CSS_SELECTOR, "#mols2grid .m2g-actions").click()
        self.find_clickable(
            By.CSS_SELECTOR, f'#mols2grid .m2g-actions option[value="{action}"]'
        ).click()

    def get_tooltip_content(self, pause=0, selector=".m2g-cell .m2g-info"):
        (
            ActionChains(self)
            .move_to_element(self.find_by_css_selector(f"#mols2grid {selector}"))
            .pause(pause)
            .perform()
        )
        tooltip = self.find_by_css_selector('div.popover[role="tooltip"]')
        el = tooltip.find_element(By.CLASS_NAME, "popover-body")
        return el.get_attribute("innerHTML")

    def trigger_callback(
        self, selector="#mols2grid .m2g-cell .m2g-callback", pause=0.2
    ):
        self.wait_for_img_load()
        el = self.find_clickable(By.CSS_SELECTOR, selector)
        (ActionChains(self).move_to_element(el).pause(pause).click().perform())

    def click_checkbox(self, is_empty=False):
        self.find_clickable(By.CSS_SELECTOR, ".m2g-cb").click()
        return self.wait_for_selection(is_empty=is_empty)
