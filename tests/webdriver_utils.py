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
    def wait_for_img_load(self, max_delay=5, selector="#mols2grid .cell .data-img svg"):
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

    def trigger_callback(self, selector="#mols2grid .cell .data-img *", pause=0.2):
        self.wait_for_img_load()
        el = self.find_clickable(By.CSS_SELECTOR, selector)
        (ActionChains(self).move_to_element(el).pause(pause).click().perform())

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

    def get_imgs_from_svgs(self, selector="#mols2grid .cell .data-img"):
        condition = EC.presence_of_all_elements_located((By.CSS_SELECTOR, selector))
        svgs = self.wait(condition)
        for svg in svgs:
            im = svg2png(bytestring=(svg.get_attribute("innerHTML")))
            yield Image.open(BytesIO(im))

    def get_png_hash(self, selector="#mols2grid .cell .data-img *"):
        img = self.find_by_css_selector(selector)
        im = Image.open(BytesIO(b64decode(img.get_attribute("src")[22:])))
        return imagehash.average_hash(im, hash_size=16)

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
