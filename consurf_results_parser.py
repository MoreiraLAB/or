###################################################################################################
############################ MOREIRA LAB - DATA DRIVEN MOLECULAR DESIGN ###########################
###################################################################################################

from __future__ import annotations

"""
Parse and integrate ConSurf "Download all files" outputs.

The current ConSurf server produces a tar.gz/zip archive containing files such as:
    <job>_consurf_grades.txt
    <job>_With_Conservation_Scores.pdb
    msa_fasta.aln
    msa_aa_variety_percentage.csv
    consurf_colored_seq.pdf

This module:
  1. extracts ConSurf archives;
  2. parses *_consurf_grades.txt;
  3. maps positions to receptor/domain/Ballesteros-Weinstein labels when metadata are available;
  4. merges conservation scores into interface-contact tables when present;
  5. creates publication-ready summary CSVs and plots.

Expected main output:
    processed_results/consurf_scores.csv
    processed_results/interface_contacts_with_consurf.csv
    summary/consurf_domain_summary.csv
    images/consurf_profile_<job>.png/svg
    images/consurf_domain_summary.png/svg
"""

import argparse
import gzip
import json
import re
import shutil
import tarfile
import zipfile
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent

try:
    from opioid_metadata import ROOT as META_ROOT
    from opioid_metadata import load_bw_numbering, bw_label, receptor_domain, parse_complex_name, ensure_dirs
except Exception:
    META_ROOT = ROOT

    def ensure_dirs(*paths):
        for p in paths:
            Path(p).mkdir(parents=True, exist_ok=True)

    def parse_complex_name(name: str):
        stem = Path(name).stem
        receptor = stem.split("-")[0] if "-" in stem else ""
        partner = stem.split("-", 1)[1] if "-" in stem else ""
        return receptor, partner, ""

    def load_bw_numbering():
        return None

    def bw_label(receptor, position, bw=None):
        return ""

    def receptor_domain(receptor, position, bw=None):
        return ""


GRADE_RE = re.compile(
    r"^\s*(?P<pos>\d+)\s+"
    r"(?P<seq>[A-Za-z\-])\s+"
    r"(?P<atom>[A-Za-z]{3}:\-?\d+[A-Za-z]?:[A-Za-z0-9])\s+"
    r"(?P<score>[-+]?\d+(?:\.\d+)?)\s+"
    r"(?P<color>[1-9]\*?)\s+"
    r"(?P<ci_low>[-+]?\d+(?:\.\d+)?)\s*,\s*"
    r"(?P<ci_high>[-+]?\d+(?:\.\d+)?)\s+"
    r"(?P<ci_colors>[1-9]\s*,\s*[1-9])\s+"
    r"(?P<be>[be])\s*"
    r"(?P<fs>[fs]?)\s+"
    r"(?P<msa>\d+/\d+)\s*"
    r"(?P<variety>.*)$"
)

ATOM_RE = re.compile(r"(?P<resname>[A-Za-z]{3}):(?P<resnum>-?\d+[A-Za-z]?):(?P<chain>[A-Za-z0-9])")


def infer_receptor_job(path: Path) -> Tuple[str, str]:
    """Infer receptor and job name from a ConSurf file or folder."""
    text = str(path)
    job = path.stem
    # Remove common suffixes.
    job = re.sub(r"_consurf.*$", "", job)
    job = re.sub(r"_With_.*$", "", job)
    job = re.sub(r"_ATOMS.*$", "", job)
    job = re.sub(r"_jmol.*$", "", job)
    for rec in ["DOR", "KOR", "MOR", "NOP"]:
        if re.search(rf"(^|[/_\-]){rec}([/_\-]|$)", text, re.I):
            return rec, job
    if "-" in job:
        return job.split("-")[0], job
    return "", job


def extract_archive(archive: Path, outdir: Optional[Path] = None) -> Path:
    """Extract .tar.gz/.tgz/.zip ConSurf output into a stable folder."""
    archive = Path(archive)
    if outdir is None:
        safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", archive.name)
        safe = re.sub(r"\.(tar\.gz|tgz|zip)$", "", safe, flags=re.I)
        outdir = ROOT / "consurf" / "results_manual" / safe
    ensure_dirs(outdir)

    if archive.suffix.lower() == ".zip":
        with zipfile.ZipFile(archive) as zf:
            zf.extractall(outdir)
    elif archive.name.lower().endswith((".tar.gz", ".tgz")):
        with tarfile.open(archive, "r:gz") as tf:
            # Python <3.12 compatible, but avoid absolute path members.
            for member in tf.getmembers():
                if member.name.startswith("/") or ".." in Path(member.name).parts:
                    continue
                tf.extract(member, outdir)
    else:
        raise ValueError(f"Unsupported archive type: {archive}")
    return outdir


def find_grade_files(input_path: Path) -> List[Path]:
    """Find ConSurf *_consurf_grades.txt files under folder or archive."""
    input_path = Path(input_path)
    if input_path.is_file() and input_path.name.lower().endswith((".tar.gz", ".tgz", ".zip")):
        input_path = extract_archive(input_path)
    if input_path.is_file():
        return [input_path] if "consurf_grades" in input_path.name.lower() else []
    return sorted([p for p in input_path.rglob("*") if p.is_file() and "consurf_grades" in p.name.lower() and p.suffix.lower() == ".txt"])


def parse_grades_file(path: Path, receptor: Optional[str] = None, job_name: Optional[str] = None) -> pd.DataFrame:
    """Parse a ConSurf *_consurf_grades.txt table."""
    path = Path(path)
    inferred_receptor, inferred_job = infer_receptor_job(path)
    receptor = receptor or inferred_receptor
    job_name = job_name or inferred_job

    rows = []
    for line in path.read_text(errors="ignore").splitlines():
        m = GRADE_RE.match(line)
        if not m:
            continue
        d = m.groupdict()
        atom = ATOM_RE.match(d["atom"])
        resname = pdb_resnum = chain = ""
        if atom:
            resname = atom.group("resname")
            pdb_resnum = atom.group("resnum")
            chain = atom.group("chain")

        color = d["color"]
        grade_low_conf = color.endswith("*")
        grade = int(color.replace("*", ""))

        ci_grade_low, ci_grade_high = [int(x.strip()) for x in d["ci_colors"].split(",")]
        msa_count, msa_total = [int(x) for x in d["msa"].split("/")]

        rows.append(
            {
                "job_name": job_name,
                "complex": job_name,
                "receptor": receptor,
                "position": int(d["pos"]),
                "aa": d["seq"],
                "atom": d["atom"],
                "atom_resname": resname,
                "pdb_resnum": pdb_resnum,
                "chain": chain,
                "score": float(d["score"]),
                "grade": grade,
                "grade_low_confidence": grade_low_conf,
                "ci_low": float(d["ci_low"]),
                "ci_high": float(d["ci_high"]),
                "ci_grade_low": ci_grade_low,
                "ci_grade_high": ci_grade_high,
                "buried_exposed": d["be"],
                "functional_structural": d["fs"],
                "msa_count": msa_count,
                "msa_total": msa_total,
                "msa_fraction": msa_count / msa_total if msa_total else None,
                "residue_variety": d["variety"].strip(),
                "source_file": str(path),
            }
        )
    df = pd.DataFrame(rows)
    if not df.empty:
        try:
            bw = load_bw_numbering()
            df["bw"] = df.apply(lambda r: bw_label(r["receptor"], int(r["position"]), bw), axis=1)
            df["domain"] = df.apply(lambda r: receptor_domain(r["receptor"], int(r["position"]), bw), axis=1)
        except Exception:
            df["bw"] = ""
            df["domain"] = ""
    return df


def collect_results(input_dir: Path = ROOT / "consurf" / "results_manual",
                    out: Path = ROOT / "processed_results" / "consurf_scores.csv") -> pd.DataFrame:
    """Collect all ConSurf grade files from a folder/archive into one CSV."""
    ensure_dirs(Path(out).parent)
    grade_files = find_grade_files(Path(input_dir))
    frames = [parse_grades_file(p) for p in grade_files]
    frames = [f for f in frames if not f.empty]
    if frames:
        df = pd.concat(frames, ignore_index=True)
    else:
        df = pd.DataFrame(
            columns=[
                "job_name", "complex", "receptor", "position", "aa", "atom", "atom_resname",
                "pdb_resnum", "chain", "score", "grade", "grade_low_confidence",
                "ci_low", "ci_high", "ci_grade_low", "ci_grade_high", "buried_exposed",
                "functional_structural", "msa_count", "msa_total", "msa_fraction",
                "residue_variety", "bw", "domain", "source_file"
            ]
        )
    df.to_csv(out, index=False)
    write_domain_summary(df)
    integrate_interface_contacts(df)
    return df


def write_domain_summary(df: pd.DataFrame, out: Path = ROOT / "summary" / "consurf_domain_summary.csv") -> Optional[pd.DataFrame]:
    """Summarize conservation scores/grades by job/receptor/domain."""
    if df.empty:
        return None
    ensure_dirs(Path(out).parent)
    group_cols = [c for c in ["job_name", "receptor", "domain"] if c in df.columns]
    if not group_cols or "domain" not in group_cols:
        return None
    d = df.copy()
    d["domain"] = d["domain"].replace("", "unmapped").fillna("unmapped")
    summary = (
        d.groupby(group_cols, dropna=False)
        .agg(
            n_residues=("position", "count"),
            mean_score=("score", "mean"),
            median_score=("score", "median"),
            mean_grade=("grade", "mean"),
            median_grade=("grade", "median"),
            high_conf_fraction=("grade_low_confidence", lambda x: 1 - float(pd.Series(x).mean())),
            buried_fraction=("buried_exposed", lambda x: float((pd.Series(x) == "b").mean())),
            exposed_fraction=("buried_exposed", lambda x: float((pd.Series(x) == "e").mean())),
        )
        .reset_index()
    )
    summary.to_csv(out, index=False)
    return summary


def integrate_interface_contacts(df: pd.DataFrame,
                                 contacts: Path = ROOT / "processed_results" / "interface_contacts.csv",
                                 out: Path = ROOT / "processed_results" / "interface_contacts_with_consurf.csv") -> Optional[pd.DataFrame]:
    """Merge conservation scores into interface contact table if present."""
    if df.empty or not Path(contacts).exists():
        return None
    c = pd.read_csv(contacts)
    if c.empty:
        return None

    # Try several likely column combinations from the reconstructed pipeline.
    candidates = [
        ("receptor", "receptor_resnum", "receptor", "position"),
        ("receptor", "resnum", "receptor", "position"),
        ("complex", "receptor_resnum", "complex", "position"),
        ("job_name", "receptor_resnum", "job_name", "position"),
    ]
    merged = None
    for left_rec, left_pos, right_rec, right_pos in candidates:
        if left_rec in c.columns and left_pos in c.columns and right_rec in df.columns and right_pos in df.columns:
            tmp = c.copy()
            tmp[left_pos] = pd.to_numeric(tmp[left_pos], errors="coerce")
            right = df.copy()
            right[right_pos] = pd.to_numeric(right[right_pos], errors="coerce")
            keep_cols = [right_rec, right_pos, "score", "grade", "grade_low_confidence", "bw", "domain", "buried_exposed", "functional_structural"]
            keep_cols = [x for x in keep_cols if x in right.columns]
            merged = tmp.merge(right[keep_cols], left_on=[left_rec, left_pos], right_on=[right_rec, right_pos], how="left", suffixes=("", "_consurf"))
            break

    if merged is None:
        return None
    ensure_dirs(Path(out).parent)
    merged.to_csv(out, index=False)
    return merged


def plot_consurf(scores: Path = ROOT / "processed_results" / "consurf_scores.csv") -> bool:
    """Generate per-job profiles and domain summaries."""
    scores = Path(scores)
    if not scores.exists():
        return False
    df = pd.read_csv(scores)
    if df.empty:
        return False
    ensure_dirs(ROOT / "images")

    for job, g in df.groupby("job_name"):
        plt.figure(figsize=(12, 3.2))
        y = g["grade"] if "grade" in g.columns and g["grade"].notna().any() else g["score"]
        plt.scatter(g["position"], y, s=14)
        plt.xlabel("Residue position")
        plt.ylabel("ConSurf grade (9 = conserved)" if y.name == "grade" else "ConSurf score")
        plt.title(f"ConSurf conservation profile: {job}")
        plt.tight_layout()
        safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(job))
        plt.savefig(ROOT / "images" / f"consurf_profile_{safe}.png", dpi=300)
        plt.savefig(ROOT / "images" / f"consurf_profile_{safe}.svg")
        plt.close()

    domain_file = ROOT / "summary" / "consurf_domain_summary.csv"
    if domain_file.exists():
        s = pd.read_csv(domain_file)
        if not s.empty:
            # Use receptor/domain mean when multiple jobs exist.
            agg = s.groupby(["receptor", "domain"], dropna=False)["mean_grade"].mean().reset_index()
            order = ["TM1", "ICL1", "TM2", "ECL1", "TM3", "ICL2", "TM4", "ECL2", "TM5", "ICL3", "TM6", "ECL3", "TM7", "H8", "unmapped"]
            agg["domain"] = pd.Categorical(agg["domain"].fillna("unmapped"), categories=order, ordered=True)
            agg = agg.sort_values(["receptor", "domain"])
            plt.figure(figsize=(12, 4))
            labels = agg["receptor"].astype(str) + ":" + agg["domain"].astype(str)
            plt.bar(range(len(agg)), agg["mean_grade"])
            plt.xticks(range(len(agg)), labels, rotation=90, fontsize=7)
            plt.ylabel("Mean ConSurf grade")
            plt.title("Mean conservation by receptor/domain")
            plt.tight_layout()
            plt.savefig(ROOT / "images" / "consurf_domain_summary.png", dpi=300)
            plt.savefig(ROOT / "images" / "consurf_domain_summary.svg")
            plt.close()
    return True


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--extract", help="Extract a ConSurf .tar.gz/.tgz/.zip archive.")
    ap.add_argument("--input-dir", default=str(ROOT / "consurf" / "results_manual"), help="Folder containing extracted ConSurf results.")
    ap.add_argument("--collect", action="store_true", help="Parse all *_consurf_grades.txt files into processed_results/consurf_scores.csv.")
    ap.add_argument("--plot", action="store_true", help="Create ConSurf plots from processed_results/consurf_scores.csv.")
    args = ap.parse_args()

    if args.extract:
        out = extract_archive(Path(args.extract))
        print(out)

    if args.collect:
        df = collect_results(Path(args.input_dir))
        print(json.dumps({"rows": int(len(df)), "output": str(ROOT / "processed_results" / "consurf_scores.csv")}, indent=2))

    if args.plot:
        print(json.dumps({"plot": bool(plot_consurf())}, indent=2))


if __name__ == "__main__":
    main()
