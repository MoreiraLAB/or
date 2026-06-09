import re, requests, time
from pathlib import Path
from urllib.parse import urljoin

ROOT = Path(".")
debug = ROOT / "consurf" / "debug"
outroot = ROOT / "consurf" / "results_manual"
outroot.mkdir(parents=True, exist_ok=True)

htmls = list(debug.rglob("*.html"))
numbers = {}

for h in htmls:
    txt = h.read_text(errors="ignore")
    m = re.search(r'id=["\']Run_Number["\'][^>]*value=["\'](\d+)["\']', txt)
    if m:
        job = h.name.replace("_before_submit.html", "").replace("_ready.html", "").replace("_after_submit.html", "")
        numbers[job] = m.group(1)

print("Jobs found:", len(numbers))

for job, num in numbers.items():
    url = f"https://consurf.tau.ac.il/progress/?number={num}"
    print(job, num, url)

    try:
        r = requests.get(url, timeout=180, allow_redirects=True)
    except Exception as e:
        print("  progress page failed:", e)
        continue

    folder = outroot / job
    folder.mkdir(parents=True, exist_ok=True)
    (folder / "progress_or_results.html").write_text(r.text, encoding="utf-8", errors="ignore")

    links = re.findall(r'href=["\']([^"\']+)["\']', r.text)
    download_links = []
    for link in links:
        full = urljoin(r.url, link)
        low = full.lower()
        if "consurfdb.tau.ac.il" in low:
            continue
        # Skip global/site links, especially the standalone software archive.
        if "stand_alone_consurf" in low or "consurfdb.tau.ac.il" in low:
            continue

        # Keep only job-specific/result-specific files.
        job_specific = (num in full) or any(x in low for x in [
            "consurf_grades",
            "with_conservation_scores",
            "msa_fasta",
            "msa_aa_variety",
            "colored_seq",
            "pymol",
            "chimera",
            "results",
            "output"
        ])

        if job_specific:
            download_links.append(full)

    saved = 0
    for i, link in enumerate(sorted(set(download_links))):
        try:
            rr = requests.get(link, timeout=180)
            name = Path(link.split("?")[0]).name or f"download_{i}.dat"
            if "." not in name:
                name = f"download_{i}.dat"
            target = folder / name
            target.write_bytes(rr.content)
            saved += 1
        except Exception as e:
            print("  failed:", link, e)

    print("  saved links:", saved)
    time.sleep(3)
