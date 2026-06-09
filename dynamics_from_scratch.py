###################################################################################################
############################ MOREIRA LAB - DATA DRIVEN MOLECULAR DESIGN ###########################
###################################################################################################

from __future__ import annotations
"""
Python replacement for the legacy R/Bio3D Normal Mode Analysis layer.

This module computes NMA-like dynamics directly from CA elastic-network models
(Gaussian Network Model style) and produces the Figure S3 family:
  S3A: mean fluctuation fold-change by receptor subdomain
  S3B: Bhattacharyya coefficient between monomer and complex fluctuation profiles
  S3C: MDS embedding of domain-level dynamics signatures

It is intentionally metadata-driven: receptor domain definitions come from
metadata/weinstein_numbering_opioids.csv and complex names from the original OR
PDB naming convention.
"""
import argparse
import math
import re
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa
from sklearn.manifold import MDS

from opioid_metadata import ROOT, ensure_dirs, load_bw_numbering, receptor_domain, bw_label, parse_complex_name

DYNAMICS_DIR = ROOT / "dynamics_complexes"
OUTDIR = ROOT / "processed_results" / "dynamics_from_scratch"
IMAGES = ROOT / "images"
SUMMARY = ROOT / "summary"
DOMAIN_ORDER = ["TM1", "ICL1", "TM2", "ECL1", "TM3", "ICL2", "TM4", "ECL2", "TM5", "ICL3", "TM6", "ECL3", "TM7", "H8", "Other"]


def parse_ca(pdb_path: Path, chain_id: str = "A") -> pd.DataFrame:
    parser = PDBParser(QUIET=True)
    st = parser.get_structure(pdb_path.stem, str(pdb_path))
    rows = []
    for model in st:
        for chain in model:
            if chain.id != chain_id:
                continue
            for res in chain:
                if not is_aa(res, standard=True) or "CA" not in res:
                    continue
                atom = res["CA"]
                rows.append({
                    "chain": chain.id,
                    "resnum": int(res.id[1]),
                    "icode": res.id[2].strip(),
                    "resname": res.resname,
                    "x": float(atom.coord[0]),
                    "y": float(atom.coord[1]),
                    "z": float(atom.coord[2]),
                })
    return pd.DataFrame(rows)


def infer_template_from_partner(partner: str | None) -> str:
    if not partner:
        return "6DDF"
    p = partner.upper()
    if "6PWC" in p:
        return "6PWC"
    if "6U1N" in p:
        return "6U1N"
    if "6OIJ" in p:
        return "6DDF"  # monomer templates available are 6DDF/6PWC/6U1N in this repository
    return "6DDF"


def monomer_for_complex(complex_pdb: Path, receptor: str, partner: str | None) -> Path | None:
    template = infer_template_from_partner(partner)
    candidates = [
        DYNAMICS_DIR / f"{receptor}_{template}.pdb",
        DYNAMICS_DIR / f"{receptor}_6DDF.pdb",
        DYNAMICS_DIR / f"{receptor}_6PWC.pdb",
        DYNAMICS_DIR / f"{receptor}_6U1N.pdb",
    ]
    for c in candidates:
        if c.exists():
            return c
    return None


def gnm_fluctuations(coords: np.ndarray, cutoff: float = 10.0) -> Tuple[np.ndarray, np.ndarray]:
    """Return per-node mean-square fluctuations and normalized cross-correlation.

    This is a compact GNM implementation: off-diagonal Kirchhoff elements are -1
    for CA pairs within cutoff, diagonal elements are contact counts. The Moore-
    Penrose pseudoinverse gives covariance up to an arbitrary gamma/kBT scaling.
    """
    n = coords.shape[0]
    if n < 3:
        return np.full(n, np.nan), np.full((n, n), np.nan)
    diff = coords[:, None, :] - coords[None, :, :]
    dist2 = np.sum(diff * diff, axis=2)
    contact = (dist2 <= cutoff * cutoff) & (dist2 > 0)
    kirchhoff = -contact.astype(float)
    np.fill_diagonal(kirchhoff, -kirchhoff.sum(axis=1))
    cov = np.linalg.pinv(kirchhoff, rcond=1e-8)
    fluc = np.diag(cov).copy()
    denom = np.sqrt(np.outer(fluc, fluc))
    with np.errstate(divide="ignore", invalid="ignore"):
        corr = cov / denom
    corr[~np.isfinite(corr)] = 0.0
    return fluc, corr


def normalize_profile(values: Iterable[float]) -> np.ndarray:
    v = np.asarray(list(values), dtype=float)
    v = v[np.isfinite(v)]
    if len(v) == 0:
        return np.array([])
    v = v - np.nanmin(v)
    if np.nansum(v) <= 0:
        v = np.ones_like(v)
    return v / np.nansum(v)


def bhattacharyya(p_values: Iterable[float], q_values: Iterable[float]) -> float:
    p = normalize_profile(p_values)
    q = normalize_profile(q_values)
    m = min(len(p), len(q))
    if m == 0:
        return np.nan
    return float(np.sum(np.sqrt(p[:m] * q[:m])))


def domain_pair_correlations(df: pd.DataFrame, corr: np.ndarray, domain_col="domain") -> pd.DataFrame:
    rows = []
    domains = [d for d in DOMAIN_ORDER if d in set(df[domain_col])]
    for d1 in domains:
        idx1 = df.index[df[domain_col] == d1].to_numpy()
        for d2 in domains:
            idx2 = df.index[df[domain_col] == d2].to_numpy()
            if len(idx1) == 0 or len(idx2) == 0:
                continue
            block = corr[np.ix_(idx1, idx2)]
            rows.append({"domain_1": d1, "domain_2": d2, "mean_correlation": float(np.nanmean(block))})
    return pd.DataFrame(rows)


def analyze_one(complex_pdb: Path, cutoff: float, bw: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    receptor, partner, partner_class = parse_complex_name(complex_pdb.name)
    if receptor not in {"DOR", "KOR", "MOR", "NOP"} or not partner:
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
    mono_pdb = monomer_for_complex(complex_pdb, receptor, partner)
    if not mono_pdb:
        print(f"No monomer found for {complex_pdb.name}; skipping dynamics.")
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()

    comp = parse_ca(complex_pdb, "A")
    mono = parse_ca(mono_pdb, "A")
    if comp.empty or mono.empty:
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()

    common = sorted(set(comp.resnum) & set(mono.resnum))
    comp = comp[comp.resnum.isin(common)].sort_values("resnum").reset_index(drop=True)
    mono = mono[mono.resnum.isin(common)].sort_values("resnum").reset_index(drop=True)
    if len(comp) < 10:
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()

    cfluc, ccorr = gnm_fluctuations(comp[["x", "y", "z"]].to_numpy(), cutoff)
    mfluc, mcorr = gnm_fluctuations(mono[["x", "y", "z"]].to_numpy(), cutoff)
    with np.errstate(divide="ignore", invalid="ignore"):
        fold = cfluc / mfluc
    res = comp[["chain", "resnum", "resname"]].copy()
    res["complex"] = complex_pdb.stem
    res["monomer"] = mono_pdb.stem
    res["receptor"] = receptor
    res["partner"] = partner
    res["partner_class"] = partner_class
    res["domain"] = res.resnum.apply(lambda x: receptor_domain(receptor, int(x), bw))
    res["bw"] = res.resnum.apply(lambda x: bw_label(receptor, int(x), bw))
    res["fluctuation_complex"] = cfluc
    res["fluctuation_monomer"] = mfluc
    res["fold_change"] = fold

    domains = []
    for dom, g in res.groupby("domain", sort=False):
        domains.append({
            "complex": complex_pdb.stem,
            "receptor": receptor,
            "partner": partner,
            "partner_class": partner_class,
            "domain": dom,
            "n_residues": len(g),
            "mean_fluctuation_complex": float(np.nanmean(g.fluctuation_complex)),
            "mean_fluctuation_monomer": float(np.nanmean(g.fluctuation_monomer)),
            "mean_fold_change": float(np.nanmean(g.fold_change)),
            "bhattacharyya_coefficient": bhattacharyya(g.fluctuation_complex, g.fluctuation_monomer),
        })
    domdf = pd.DataFrame(domains)
    corrdf = domain_pair_correlations(res, ccorr)
    if not corrdf.empty:
        corrdf.insert(0, "complex", complex_pdb.stem)
        corrdf.insert(1, "receptor", receptor)
        corrdf.insert(2, "partner", partner)
        corrdf.insert(3, "partner_class", partner_class)
    return res, domdf, corrdf


def run_all(cutoff: float = 10.0, limit: int | None = None) -> Dict[str, Path]:
    ensure_dirs(OUTDIR, IMAGES, SUMMARY)
    bw = load_bw_numbering()
    residue_frames: List[pd.DataFrame] = []
    domain_frames: List[pd.DataFrame] = []
    corr_frames: List[pd.DataFrame] = []
    pdbs = sorted(p for p in DYNAMICS_DIR.glob("*.pdb") if "-" in p.stem)
    if limit:
        pdbs = pdbs[:limit]
    for pdb in pdbs:
        try:
            r, d, c = analyze_one(pdb, cutoff, bw)
            if not r.empty:
                residue_frames.append(r)
            if not d.empty:
                domain_frames.append(d)
            if not c.empty:
                corr_frames.append(c)
        except Exception as e:
            print(f"WARNING: dynamics failed for {pdb.name}: {e}")
    residues = pd.concat(residue_frames, ignore_index=True) if residue_frames else pd.DataFrame()
    domains = pd.concat(domain_frames, ignore_index=True) if domain_frames else pd.DataFrame()
    corrs = pd.concat(corr_frames, ignore_index=True) if corr_frames else pd.DataFrame()

    residue_csv = OUTDIR / "residue_fluctuations.csv"
    domain_csv = OUTDIR / "domain_dynamics_metrics.csv"
    corr_csv = OUTDIR / "domain_correlation_metrics.csv"
    residues.to_csv(residue_csv, index=False)
    domains.to_csv(domain_csv, index=False)
    corrs.to_csv(corr_csv, index=False)
    make_s3_figures(domains, corrs)
    return {"residue_fluctuations": residue_csv, "domain_metrics": domain_csv, "correlations": corr_csv}


def make_s3_figures(domains: pd.DataFrame, corrs: pd.DataFrame) -> None:
    if domains.empty:
        return
    # S3A: average fluctuation fold changes by domain, grouped by partner class
    data = domains.copy()
    data["domain"] = pd.Categorical(data["domain"], DOMAIN_ORDER, ordered=True)
    pivot = data.pivot_table(index="domain", columns="partner_class", values="mean_fold_change", aggfunc="mean").dropna(how="all")
    plt.figure(figsize=(10, 4.8))
    x = np.arange(len(pivot.index))
    width = 0.38
    cols = list(pivot.columns)
    for i, col in enumerate(cols):
        off = (i - (len(cols)-1)/2) * width
        plt.bar(x + off, pivot[col].values, width=width, label=col)
    plt.xticks(x, pivot.index, rotation=45, ha="right")
    plt.ylabel("Mean fluctuation fold-change")
    plt.title("Figure S3A-style dynamics: average fluctuation fold-change")
    plt.legend()
    plt.tight_layout()
    plt.savefig(IMAGES / "figure_s3A_fluctuation_foldchange.png", dpi=300)
    plt.savefig(IMAGES / "figure_s3A_fluctuation_foldchange.svg")
    plt.close()

    # S3B: Bhattacharyya coefficients by domain
    pivot = data.pivot_table(index="domain", columns="partner_class", values="bhattacharyya_coefficient", aggfunc="mean").dropna(how="all")
    plt.figure(figsize=(10, 4.8))
    x = np.arange(len(pivot.index))
    for i, col in enumerate(pivot.columns):
        off = (i - (len(pivot.columns)-1)/2) * width
        plt.bar(x + off, pivot[col].values, width=width, label=col)
    plt.xticks(x, pivot.index, rotation=45, ha="right")
    plt.ylabel("Bhattacharyya coefficient")
    plt.ylim(0, 1.05)
    plt.title("Figure S3B-style dynamics: monomer-complex similarity")
    plt.legend()
    plt.tight_layout()
    plt.savefig(IMAGES / "figure_s3B_bhattacharyya.png", dpi=300)
    plt.savefig(IMAGES / "figure_s3B_bhattacharyya.svg")
    plt.close()

    # S3C: MDS on per-complex domain signature
    sig = data.pivot_table(index=["complex", "receptor", "partner", "partner_class"], columns="domain", values="bhattacharyya_coefficient", aggfunc="mean")
    sig = sig.fillna(sig.mean(numeric_only=True)).fillna(0)
    if len(sig) >= 3 and sig.shape[1] >= 2:
        emb = MDS(n_components=2, dissimilarity="euclidean", random_state=7).fit_transform(sig.to_numpy())
        mds = sig.reset_index()[["complex", "receptor", "partner", "partner_class"]].copy()
        mds["MDS1"] = emb[:, 0]
        mds["MDS2"] = emb[:, 1]
        mds.to_csv(OUTDIR / "dynamics_mds_coordinates.csv", index=False)
        plt.figure(figsize=(7, 6))
        for cls, g in mds.groupby("partner_class"):
            plt.scatter(g.MDS1, g.MDS2, s=70, alpha=0.75, label=cls)
            for _, row in g.iterrows():
                plt.text(row.MDS1, row.MDS2, row.receptor, fontsize=6, alpha=.8)
        plt.xlabel("MDS1")
        plt.ylabel("MDS2")
        plt.title("Figure S3C-style dynamics: MDS of domain BC signatures")
        plt.legend()
        plt.tight_layout()
        plt.savefig(IMAGES / "figure_s3C_dynamics_mds.png", dpi=300)
        plt.savefig(IMAGES / "figure_s3C_dynamics_mds.svg")
        plt.close()

    # Correlated motions heatmap by domain pairs
    if not corrs.empty:
        p = corrs.pivot_table(index="domain_1", columns="domain_2", values="mean_correlation", aggfunc="mean")
        p = p.reindex(index=[d for d in DOMAIN_ORDER if d in p.index], columns=[d for d in DOMAIN_ORDER if d in p.columns])
        plt.figure(figsize=(7, 6))
        plt.imshow(p.values, vmin=-1, vmax=1, aspect="auto")
        plt.xticks(range(len(p.columns)), p.columns, rotation=90)
        plt.yticks(range(len(p.index)), p.index)
        plt.colorbar(label="Mean CA cross-correlation")
        plt.title("Domain-level correlated motions from CA elastic network")
        plt.tight_layout()
        plt.savefig(IMAGES / "figure_s3_domain_correlated_motions.png", dpi=300)
        plt.savefig(IMAGES / "figure_s3_domain_correlated_motions.svg")
        plt.close()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cutoff", type=float, default=10.0, help="CA contact cutoff for the elastic network, in Å")
    ap.add_argument("--limit", type=int, default=None, help="Debug: only process the first N complexes")
    args = ap.parse_args()
    print(run_all(cutoff=args.cutoff, limit=args.limit))


if __name__ == "__main__":
    main()
