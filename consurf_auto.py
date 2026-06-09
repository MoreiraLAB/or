###################################################################################################
############################ MOREIRA LAB - DATA DRIVEN MOLECULAR DESIGN ###########################
###################################################################################################

from __future__ import annotations

"""
Robust ConSurf submitter for the current https://consurf.tau.ac.il/ front-end.

This version is intentionally specific to the HTML observed in June 2026:
  - file input:          #browse / input[name='pdb_FILE']
  - chain select:        #PDB_chain / select[name='PDB_chain']
  - job title:           #JOB_TITLE
  - email:               #user_email
  - default submit btn:  #startJob
  - blocking popup:      #popup.active, close button #popup .close

Important: the chain select appears only after the PDB upload has been processed by
ConSurf. This script waits until #PDB_chain contains options, then sets its value by
JavaScript and dispatches change/input events.
"""

import argparse
import json
import re
import time
import requests
from urllib.parse import urljoin
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd

CONSURF_URL = "https://consurf.tau.ac.il/"
ROOT = Path(__file__).resolve().parent


def ensure_dirs(*paths: Path) -> None:
    for p in paths:
        Path(p).mkdir(parents=True, exist_ok=True)


def _setup_driver(headless: bool = False, download_dir: Optional[Path] = None):
    from selenium import webdriver
    from selenium.webdriver.chrome.options import Options

    opts = Options()
    if headless:
        opts.add_argument("--headless=new")
    opts.add_argument("--no-sandbox")
    opts.add_argument("--disable-dev-shm-usage")
    opts.add_argument("--window-size=1400,1000")
    if download_dir:
        ensure_dirs(download_dir)
        opts.add_experimental_option(
            "prefs",
            {
                "download.default_directory": str(download_dir.resolve()),
                "download.prompt_for_download": False,
                "download.directory_upgrade": True,
                "safebrowsing.enabled": True,
            },
        )
    # Prefer Selenium Manager (Selenium >= 4.6). Fall back to webdriver-manager only if needed.
    try:
        driver = webdriver.Chrome(options=opts)
    except Exception:
        from selenium.webdriver.chrome.service import Service
        from webdriver_manager.chrome import ChromeDriverManager
        driver = webdriver.Chrome(service=Service(ChromeDriverManager().install()), options=opts)
    return driver


def export_pdb(pdb_dir: Path = ROOT / "structural_complexes", outdir: Path = ROOT / "consurf" / "input", chain: str = "A") -> Path:
    """Create a PDB manifest for ConSurf structure submission."""
    ensure_dirs(outdir)
    rows = []
    for p in sorted(Path(pdb_dir).glob("*.pdb")):
        rows.append(
            {
                "complex": p.stem,
                "pdb": str(p.resolve()),
                "chain": chain,
                "job_name": p.stem,
            }
        )
    manifest = outdir / "consurf_pdb_manifest.csv"
    pd.DataFrame(rows).to_csv(manifest, index=False)
    return manifest


def _dismiss_alerts(driver) -> List[str]:
    messages = []
    for _ in range(3):
        try:
            alert = driver.switch_to.alert
            msg = alert.text
            messages.append(msg)
            alert.accept()
            time.sleep(0.5)
        except Exception:
            break
    return messages


def _close_popup(driver) -> bool:
    """Close ConSurf popup such as 'remark removed' if present."""
    closed = False
    try:
        driver.execute_script(
            """
            const pop = document.querySelector('#popup');
            if (pop) { pop.classList.remove('active'); pop.style.display='none'; }
            const overlay = document.querySelector('.popup.active');
            if (overlay) { overlay.classList.remove('active'); overlay.style.display='none'; }
            """
        )
        closed = True
    except Exception:
        pass
    try:
        close = driver.find_element("css selector", "#popup .close")
        driver.execute_script("arguments[0].click();", close)
        closed = True
    except Exception:
        pass
    return closed


def _set_text_by_id(driver, element_id: str, value: str) -> bool:
    try:
        el = driver.find_element("id", element_id)
        driver.execute_script(
            """
            arguments[0].focus();
            arguments[0].value = arguments[1];
            arguments[0].dispatchEvent(new Event('input', {bubbles:true}));
            arguments[0].dispatchEvent(new Event('change', {bubbles:true}));
            arguments[0].blur();
            """,
            el,
            value,
        )
        return True
    except Exception:
        return False


def _upload_pdb(driver, pdb_file: Path) -> bool:
    selectors = ["#browse", "input[name='pdb_FILE']", "input[type='file']"]
    for sel in selectors:
        try:
            inp = driver.find_element("css selector", sel)
            inp.send_keys(str(Path(pdb_file).resolve()))
            return True
        except Exception:
            continue
    return False


def _wait_for_chain_options(driver, timeout: int = 60) -> List[str]:
    """Wait until ConSurf has parsed the uploaded PDB and populated #PDB_chain."""
    end = time.time() + timeout
    last_options: List[str] = []
    while time.time() < end:
        _dismiss_alerts(driver)
        _close_popup(driver)
        try:
            opts = driver.execute_script(
                """
                const s = document.querySelector('#PDB_chain, select[name="PDB_chain"]');
                if (!s) return [];
                return Array.from(s.options).map(o => ({value:o.value, text:o.textContent.trim()}));
                """
            )
            last_options = [f"{o.get('value','')}::{o.get('text','')}" for o in opts or []]
            non_empty = [o for o in (opts or []) if (o.get("value") or "").strip()]
            if non_empty:
                return last_options
        except Exception:
            pass
        time.sleep(1)
    return last_options


def _select_chain(driver, chain: str = "A", timeout: int = 60) -> bool:
    options = _wait_for_chain_options(driver, timeout=timeout)
    # Direct JS on real select is most reliable; visual dropdown is just decoration.
    try:
        ok = driver.execute_script(
            """
            const chain = arguments[0];
            const s = document.querySelector('#PDB_chain, select[name="PDB_chain"]');
            if (!s) return false;
            let found = false;
            for (const opt of Array.from(s.options)) {
                const txt = (opt.textContent || '').trim();
                const val = (opt.value || '').trim();
                if (val === chain || txt === ('Chain ' + chain) || txt === chain) {
                    s.value = opt.value;
                    found = true;
                    break;
                }
            }
            if (!found) return false;
            s.dispatchEvent(new Event('input', {bubbles:true}));
            s.dispatchEvent(new Event('change', {bubbles:true}));
            if (typeof IsChainInDB === 'function') { try { IsChainInDB(' db_link'); } catch(e) {} }
            const nice = document.querySelector('.dropdown-select .current');
            if (nice) nice.textContent = 'Chain ' + chain;
            return s.value === chain;
            """,
            chain,
        )
        if ok:
            time.sleep(1)
            return True
    except Exception:
        pass

    # Fallback: click custom dropdown text.
    try:
        driver.execute_script("document.querySelector('.dropdown-select')?.click();")
        time.sleep(0.5)
        opts = driver.find_elements("xpath", f"//*[normalize-space(text())='Chain {chain}']")
        for o in opts:
            try:
                driver.execute_script("arguments[0].click();", o)
                time.sleep(1)
                val = driver.execute_script("return document.querySelector('#PDB_chain')?.value || '';")
                if val == chain:
                    return True
            except Exception:
                continue
    except Exception:
        pass
    return False


def _save_debug(driver, job_name: str, label: str) -> None:
    d = ROOT / "consurf" / "debug"
    ensure_dirs(d)
    safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", job_name)
    try:
        (d / f"{safe}_{label}.html").write_text(driver.page_source, encoding="utf-8", errors="ignore")
    except Exception:
        pass
    try:
        driver.save_screenshot(str(d / f"{safe}_{label}.png"))
    except Exception:
        pass




def _get_run_number_from_page(driver) -> str:
    """Extract ConSurf Run_Number/job number from DOM, URL, or page source."""
    candidates = []
    try:
        val = driver.execute_script("return document.querySelector('#Run_Number')?.value || '';" )
        if val:
            candidates.append(str(val))
    except Exception:
        pass
    try:
        url = driver.current_url or ""
        m = re.search(r"(?:number|Run_Number|run_number)=([0-9]+)", url, re.I)
        if m:
            candidates.append(m.group(1))
    except Exception:
        pass
    try:
        html = driver.page_source or ""
        patterns = [
            "id=[\"']Run_Number[\"'][^>]*value=[\"']([0-9]+)[\"']",
            "name=[\"']Run_Number[\"'][^>]*value=[\"']([0-9]+)[\"']",
            "progress/\\?number=([0-9]+)",
        ]
        for pat in patterns:
            m = re.search(pat, html, re.I)
            if m:
                candidates.append(m.group(1))
    except Exception:
        pass
    for c in candidates:
        if c and str(c).isdigit():
            return str(c)
    return ""


def _progress_url(run_number: str) -> str:
    return f"{CONSURF_URL.rstrip('/')}/progress/?number={run_number}"


def submit_structure(
    pdb_file: Path,
    chain: str = "A",
    email: Optional[str] = None,
    job_name: Optional[str] = None,
    click_submit: bool = False,
    headless: bool = False,
    keep_browser: bool = True,
    wait_after_upload: int = 60,
) -> Dict[str, str]:
    job_name = job_name or Path(pdb_file).stem
    out: Dict[str, str] = {
        "input_file": str(Path(pdb_file).resolve()),
        "job_name": job_name,
        "chain": chain,
        "email": email or "",
        "server": CONSURF_URL,
        "status": "opened",
        "job_url": "",
        "final_url": "",
        "clicked_submit": str(bool(click_submit)),
        "uploaded": "False",
        "job_name_filled": "False",
        "email_filled": "False",
        "chain_filled": "False",
        "alerts": "",
        "note": "",
        "run_number": "",
        "progress_url": "",
    }
    driver = _setup_driver(headless=headless, download_dir=ROOT / "consurf" / "downloads")
    try:
        driver.get(CONSURF_URL)
        time.sleep(4)
        out["run_number"] = _get_run_number_from_page(driver)
        if out["run_number"]:
            out["progress_url"] = _progress_url(out["run_number"])
        _dismiss_alerts(driver)
        _close_popup(driver)

        out["uploaded"] = str(_upload_pdb(driver, Path(pdb_file)))
        # The popup may appear immediately after upload ('remark removed'); close/accept repeatedly.
        time.sleep(2)
        alerts = _dismiss_alerts(driver)
        _close_popup(driver)
        time.sleep(2)
        alerts.extend(_dismiss_alerts(driver))
        _close_popup(driver)

        out["job_name_filled"] = str(_set_text_by_id(driver, "JOB_TITLE", job_name))
        if email:
            out["email_filled"] = str(_set_text_by_id(driver, "user_email", email))

        out["chain_filled"] = str(_select_chain(driver, chain=chain, timeout=wait_after_upload))
        alerts.extend(_dismiss_alerts(driver))
        _close_popup(driver)
        out["alerts"] = " | ".join([a for a in alerts if a])

        _save_debug(driver, job_name, "before_submit")

        if click_submit:
            if out["chain_filled"] != "True":
                out["status"] = "not_submitted_chain_not_selected"
                out["note"] = "Chain was not selected; browser left open for manual inspection."
            else:
                _close_popup(driver)
                time.sleep(0.5)
                try:
                    btn = driver.find_element("id", "startJob")
                    driver.execute_script("arguments[0].scrollIntoView({block:'center'});", btn)
                    time.sleep(0.5)
                    driver.execute_script("arguments[0].click();", btn)
                    time.sleep(8)
                    alerts.extend(_dismiss_alerts(driver))
                    rn = _get_run_number_from_page(driver) or out.get("run_number", "")
                    if rn:
                        out["run_number"] = rn
                        out["progress_url"] = _progress_url(rn)
                    out["status"] = "submitted_or_ready"
                except Exception as e:
                    out["status"] = "submit_click_failed"
                    out["note"] = str(e)
        else:
            out["status"] = "ready_for_user_review"
            out["note"] = "Browser left open so the user can verify required fields and click Submit manually."

        out["final_url"] = driver.current_url
        rn = _get_run_number_from_page(driver) or out.get("run_number", "")
        if rn:
            out["run_number"] = rn
            out["progress_url"] = _progress_url(rn)
        if re.search(r"job|result|output|status|run|progress", driver.current_url, re.I):
            out["job_url"] = driver.current_url
        elif out.get("progress_url"):
            out["job_url"] = out["progress_url"]
        _save_debug(driver, job_name, "after_submit" if click_submit else "ready")

        if headless or not keep_browser:
            driver.quit()
        return out
    except Exception as e:
        out["status"] = "error"
        out["note"] = repr(e)
        try:
            _save_debug(driver, job_name, "error")
        except Exception:
            pass
        if headless or not keep_browser:
            try:
                driver.quit()
            except Exception:
                pass
        return out


def submit_pdb_manifest(manifest: Path, email: str, chain: str = "A", limit: Optional[int] = None, click_submit: bool = False, headless: bool = False) -> List[Dict[str, str]]:
    df = pd.read_csv(manifest)
    if limit:
        df = df.head(limit)
    logs = []
    for _, row in df.iterrows():
        logs.append(
            submit_structure(
                Path(row["pdb"]),
                chain=str(row.get("chain", chain) or chain),
                email=email,
                job_name=str(row.get("job_name", Path(row["pdb"]).stem)),
                click_submit=click_submit,
                headless=headless,
                keep_browser=not headless,
            )
        )
    out = ROOT / "consurf" / "submission_log.json"
    ensure_dirs(out.parent)
    existing = []
    if out.exists():
        try:
            existing = json.loads(out.read_text(encoding="utf-8"))
        except Exception:
            existing = []
    by_job = {str(x.get("job_name", "")): x for x in existing if isinstance(x, dict)}
    for item in logs:
        by_job[str(item.get("job_name", ""))] = item
    merged = list(by_job.values())
    out.write_text(json.dumps(merged, indent=2), encoding="utf-8")
    return logs



def download_current_url(url: str, job_name: str, outdir: Optional[Path] = None) -> Path:
    """
    Download a ConSurf result/progress page and all linked files that look like
    downloadable result assets. Works also if the URL is still a progress page:
    it saves the page and any available links, but does not wait for completion.
    """
    outdir = outdir or (ROOT / "consurf" / "results_manual" / re.sub(r"[^A-Za-z0-9_.-]+", "_", job_name))
    ensure_dirs(outdir)
    session = requests.Session()
    r = session.get(url, timeout=60)
    (outdir / "page.html").write_text(r.text, encoding="utf-8", errors="ignore")
    (outdir / "source_url.txt").write_text(url + "\n", encoding="utf-8")

    links = re.findall(r"href=[\"']([^\"']+)[\"']", r.text, flags=re.I)
    wanted = []
    for link in links:
        full = urljoin(url, link)
        low = full.lower()
        label = link.lower()
        if "stand_alone_consurf" in low or "consurfdb.tau.ac.il" in low:
            continue
        if any(k in low or k in label for k in [
            "download", "zip", "grades", "score", "result", "output",
            "with_conservation_scores", "consurf_grades",
            "pdb", "cif", "txt", "csv", "dat", "pdf", "msa", "fasta", "aln",
            "pymol", "chimera", "jsmol"
        ]):
            wanted.append(full)

    manifest_rows = []
    seen = set()
    for i, full in enumerate(wanted):
        if full in seen:
            continue
        seen.add(full)
        try:
            rr = session.get(full, timeout=60)
            name = Path(full.split("?")[0]).name
            if not name or "." not in name:
                ctype = rr.headers.get("content-type", "").lower()
                ext = ".html" if "html" in ctype else ".txt"
                name = f"linked_{i:03d}{ext}"
            safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", name)
            target = outdir / safe
            target.write_bytes(rr.content)
            manifest_rows.append({"url": full, "status_code": rr.status_code, "file": str(target)})
        except Exception as e:
            manifest_rows.append({"url": full, "status_code": "", "file": "", "error": repr(e)})

    pd.DataFrame(manifest_rows).to_csv(outdir / "download_manifest.csv", index=False)
    return outdir


def collect_results(input_dir: Optional[Path] = None) -> None:
    """Collect ConSurf Download-all outputs using the archive-aware parser."""
    try:
        from consurf_results_parser import collect_results as _collect
        _collect(input_dir or (ROOT / "consurf" / "results_manual"))
        print("ConSurf collection completed.")
    except Exception as e:
        raise RuntimeError(f"Could not collect ConSurf results via consurf_results_parser.py: {e}")


def plot_results() -> None:
    """Plot parsed ConSurf outputs using the archive-aware parser."""
    try:
        from consurf_results_parser import plot_consurf as _plot
        _plot()
        print("ConSurf plotting completed.")
    except Exception as e:
        raise RuntimeError(f"Could not plot ConSurf results via consurf_results_parser.py: {e}")



def _looks_complete(html: str, final_url: str = "") -> bool:
    h = (html or "").lower()
    u = (final_url or "").lower()
    complete_markers = [
        "download all files", "consurf_grades", "with_conservation_scores",
        "colored sequence", "conservation scores", "final results"
    ]
    if any(m in h for m in complete_markers):
        return True
    if "results" in u and "progress" not in u:
        return True
    return False


def _extract_downloaded_archives(folder: Path) -> List[str]:
    extracted = []
    try:
        from consurf_results_parser import extract_archive
    except Exception:
        return extracted
    for p in Path(folder).rglob("*"):
        if p.is_file() and p.name.lower().endswith((".tar.gz", ".tgz", ".zip")):
            try:
                target = Path(folder) / re.sub(r"[^A-Za-z0-9_.-]+", "_", p.name.replace(".tar.gz", "").replace(".tgz", "").replace(".zip", ""))
                out = extract_archive(p, outdir=target)
                extracted.append(str(out))
            except Exception:
                pass
    return extracted


def watch_job(job_number_or_url: str, job_name: str, wait_minutes: int = 180, poll_seconds: int = 60, collect: bool = False, plot: bool = False) -> Dict[str, str]:
    """Poll a ConSurf progress/results URL until outputs appear, then download them."""
    token = str(job_number_or_url)
    if re.fullmatch(r"\d+", token):
        url = _progress_url(token)
        run_number = token
    else:
        url = token
        m = re.search(r"number=([0-9]+)", url)
        run_number = m.group(1) if m else ""
    deadline = time.time() + wait_minutes * 60
    status = {"job_name": job_name, "run_number": run_number, "url": url, "status": "watching", "download_dir": "", "note": ""}
    while time.time() < deadline:
        try:
            r = requests.get(url, timeout=60, allow_redirects=True)
            final_url = r.url
            html = r.text
            status.update({"last_url": final_url, "http_status": str(r.status_code)})
            if _looks_complete(html, final_url):
                outdir = download_current_url(final_url, job_name)
                extracted = _extract_downloaded_archives(outdir)
                status.update({"status": "downloaded", "download_dir": str(outdir), "extracted": ";".join(extracted)})
                if collect:
                    collect_results()
                if plot:
                    plot_results()
                return status
        except Exception as e:
            status["note"] = repr(e)
        time.sleep(poll_seconds)
    status["status"] = "timeout"
    status["note"] = f"No completed ConSurf results detected within {wait_minutes} minutes. Check {url} manually."
    return status


def watch_submission_log(log_file: Path = ROOT / "consurf" / "submission_log.json", wait_minutes: int = 180, poll_seconds: int = 60, collect: bool = False, plot: bool = False) -> List[Dict[str, str]]:
    if not Path(log_file).exists():
        raise FileNotFoundError(log_file)
    logs = json.loads(Path(log_file).read_text(encoding="utf-8"))
    results = []
    for item in logs:
        if not isinstance(item, dict):
            continue
        job_name = str(item.get("job_name") or item.get("input_file", "consurf_job"))
        token = item.get("run_number") or item.get("progress_url") or item.get("job_url") or item.get("final_url")
        if not token:
            results.append({"job_name": job_name, "status": "skipped_no_job_number_or_url"})
            continue
        results.append(watch_job(str(token), job_name, wait_minutes=wait_minutes, poll_seconds=poll_seconds, collect=False, plot=False))
    out = ROOT / "consurf" / "watch_log.json"
    ensure_dirs(out.parent)
    out.write_text(json.dumps(results, indent=2), encoding="utf-8")
    if collect:
        collect_results()
    if plot:
        plot_results()
    return results


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--export-pdb", action="store_true")
    ap.add_argument("--submit-pdb", action="store_true")
    ap.add_argument("--email", default=None)
    ap.add_argument("--chain", default="A")
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--click-submit", action="store_true")
    ap.add_argument("--headless", action="store_true")
    ap.add_argument("--download-current-url", default=None, help="ConSurf progress/results URL to download.")
    ap.add_argument("--job-name", default=None, help="Job name used for downloaded ConSurf files.")
    ap.add_argument("--collect", action="store_true", help="Parse downloaded ConSurf files through consurf_pipeline.py.")
    ap.add_argument("--plot", action="store_true", help="Plot parsed ConSurf scores through consurf_pipeline.py.")
    ap.add_argument("--watch-job", default=None, help="ConSurf job number or progress/results URL to poll until completion, then download.")
    ap.add_argument("--watch-submission-log", action="store_true", help="Poll every submitted job stored in consurf/submission_log.json.")
    ap.add_argument("--wait-minutes", type=int, default=180, help="Maximum minutes to wait for ConSurf completion.")
    ap.add_argument("--poll-seconds", type=int, default=60, help="Polling interval in seconds for ConSurf job monitoring.")
    args = ap.parse_args()

    manifest = ROOT / "consurf" / "input" / "consurf_pdb_manifest.csv"
    if args.export_pdb:
        manifest = export_pdb(chain=args.chain)
        print(manifest)
    if args.submit_pdb:
        if not manifest.exists():
            manifest = export_pdb(chain=args.chain)
        logs = submit_pdb_manifest(
            manifest,
            email=args.email or "",
            chain=args.chain,
            limit=args.limit,
            click_submit=args.click_submit,
            headless=args.headless,
        )
        print(json.dumps(logs, indent=2))

    if args.download_current_url:
        job = args.job_name or "consurf_job"
        out = download_current_url(args.download_current_url, job)
        print(f"Downloaded ConSurf page/assets to: {out}")

    if args.watch_job:
        job = args.job_name or args.watch_job
        res = watch_job(args.watch_job, job_name=job, wait_minutes=args.wait_minutes, poll_seconds=args.poll_seconds, collect=args.collect, plot=args.plot)
        print(json.dumps(res, indent=2))
        return

    if args.watch_submission_log:
        res = watch_submission_log(wait_minutes=args.wait_minutes, poll_seconds=args.poll_seconds, collect=args.collect, plot=args.plot)
        print(json.dumps(res, indent=2))
        return

    if args.collect:
        collect_results()

    if args.plot:
        plot_results()


if __name__ == "__main__":
    main()
