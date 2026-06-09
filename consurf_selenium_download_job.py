###################################################################################################
############################ MOREIRA LAB - DATA DRIVEN MOLECULAR DESIGN ###########################
###################################################################################################

from pathlib import Path
import argparse, time
from selenium import webdriver
from selenium.webdriver.chrome.options import Options

p = argparse.ArgumentParser()
p.add_argument("--number", required=True)
p.add_argument("--job-name", required=True)
p.add_argument("--wait", type=int, default=90)
args = p.parse_args()

out = Path("consurf/results_manual") / args.job_name
out.mkdir(parents=True, exist_ok=True)

opts = Options()
opts.add_experimental_option("prefs", {
    "download.default_directory": str(out.resolve()),
    "download.prompt_for_download": False,
    "download.directory_upgrade": True,
    "safebrowsing.enabled": True,
})
try:
    driver = webdriver.Chrome(options=opts)
except Exception:
    from selenium.webdriver.chrome.service import Service
    from webdriver_manager.chrome import ChromeDriverManager
    driver = webdriver.Chrome(service=Service(ChromeDriverManager().install()), options=opts)

url = f"https://consurf.tau.ac.il/final_output/?number={args.number}"
driver.get(url)
time.sleep(args.wait)
(out / "selenium_page.html").write_text(driver.page_source, encoding="utf-8", errors="ignore")
driver.save_screenshot(str(out / "selenium_page.png"))

clicked = False
for txt in ["Download all files", "download all files", "Download All Files"]:
    try:
        el = driver.find_element("xpath", f"//*[contains(normalize-space(.), '{txt}')]")
        driver.execute_script("arguments[0].click();", el)
        clicked = True
        time.sleep(30)
        break
    except Exception:
        pass

if not clicked:
    elems = driver.find_elements("xpath", "//*[contains(translate(., 'DOWNLOADRESULTFILES', 'downloadresultfiles'), 'download') or contains(translate(., 'DOWNLOADRESULTFILES', 'downloadresultfiles'), 'result')]")
    for el in elems:
        try:
            driver.execute_script("arguments[0].click();", el)
            clicked = True
            time.sleep(20)
            break
        except Exception:
            pass

print("clicked_download:", clicked)
print("saved_to:", out)
print("files:")
for f in out.iterdir():
    print(" -", f.name)
