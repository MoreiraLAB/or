###################################################################################################
############################ MOREIRA LAB - DATA DRIVEN MOLECULAR DESIGN ###########################
###################################################################################################

from __future__ import annotations
from pathlib import Path
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.manifold import MDS
from opioid_metadata import ROOT, ensure_dirs

DOMAINS = ["TM1","TM2","TM3","TM4","TM5","TM6","TM7","H8"]
PROCESSED = ROOT / "processed_results" / "dynamics"
SUMMARY = ROOT / "summary"
IMAGES = ROOT / "images"


def _find_csv(kind: str, family: str) -> Path | None:
    if not PROCESSED.exists():
        return None
    patterns = [f"*{family}*{kind}*.csv", f"*{kind}*{family}*.csv"]
    for pat in patterns:
        hits = sorted(PROCESSED.glob(pat))
        hits = [h for h in hits if "manifest" not in h.name]
        if hits:
            return hits[0]
    return None


def _standardize_bc(path: Path, family: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    # The R exporter should already add Complex; fall back to first column if needed.
    if "Complex" not in df.columns:
        first = df.columns[0]
        if not np.issubdtype(df[first].dtype, np.number):
            df = df.rename(columns={first: "Complex"})
        else:
            df.insert(0, "Complex", [f"{family}_{i+1}" for i in range(len(df))])
    numeric = [c for c in df.columns if c != "Complex" and pd.api.types.is_numeric_dtype(df[c])]
    if len(numeric) >= 8:
        rename = {numeric[i]: f"BC_{DOMAINS[i]}" for i in range(8)}
        df = df.rename(columns=rename)
    keep = ["Complex"] + [f"BC_{d}" for d in DOMAINS if f"BC_{d}" in df.columns]
    df = df[keep].copy()
    df["DXR"] = df["Complex"].str.extract(r"^(DOR|MOR|KOR|NOP)", expand=False)
    df["Partner"] = df["Complex"].str.replace(r"^(DOR|MOR|KOR|NOP)-", "", regex=True)
    df["family"] = family
    return df


def _standardize_fluc(path: Path, family: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    # Expected columns: Fluc_TM1..Fluc_H8, Complex, DXR, Partner
    fluc_cols = [c for c in df.columns if c.startswith("Fluc_")]
    if len(fluc_cols) < 8:
        nums = [c for c in df.columns if pd.api.types.is_numeric_dtype(df[c])]
        for i, c in enumerate(nums[:8]):
            df = df.rename(columns={c: f"Fluc_{DOMAINS[i]}"})
    if "Complex" not in df.columns:
        # Try row names/first nonnumeric; otherwise infer placeholder.
        nonnum = [c for c in df.columns if not pd.api.types.is_numeric_dtype(df[c])]
        if nonnum:
            df = df.rename(columns={nonnum[0]: "Complex"})
        else:
            df["Complex"] = [f"{family}_{i+1}" for i in range(len(df))]
    if "DXR" not in df.columns:
        df["DXR"] = df["Complex"].str.extract(r"^(DOR|MOR|KOR|NOP)", expand=False)
    if "Partner" not in df.columns:
        df["Partner"] = df["Complex"].str.replace(r"^(DOR|MOR|KOR|NOP)-", "", regex=True)
    keep = [f"Fluc_{d}" for d in DOMAINS] + ["Complex","DXR","Partner"]
    df = df[[c for c in keep if c in df.columns]].copy()
    df["family"] = family
    return df


def load_legacy_dynamics() -> dict[str, pd.DataFrame]:
    out = {}
    for family in ["gprot", "arrestin"]:
        bc = _find_csv("bc_scores", family)
        fluc = _find_csv("fluc_fold_change", family)
        if bc:
            out[f"{family}_bc"] = _standardize_bc(bc, family)
        if fluc:
            out[f"{family}_fluc"] = _standardize_fluc(fluc, family)
    return out


def write_standardized_tables() -> dict[str, Path]:
    ensure_dirs(SUMMARY, IMAGES)
    tables = load_legacy_dynamics()
    paths = {}
    for name, df in tables.items():
        p = SUMMARY / f"{name}_standardized.csv"
        df.to_csv(p, index=False)
        paths[name] = p
    if "gprot_bc" in tables and "arrestin_bc" in tables:
        pd.concat([tables["gprot_bc"], tables["arrestin_bc"]], ignore_index=True).to_csv(SUMMARY / "dynamics_bc_scores_all.csv", index=False)
    if "gprot_fluc" in tables and "arrestin_fluc" in tables:
        pd.concat([tables["gprot_fluc"], tables["arrestin_fluc"]], ignore_index=True).to_csv(SUMMARY / "dynamics_fluc_fold_change_all.csv", index=False)
    return paths


def _savefig(path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(path.with_suffix(".png"), dpi=300)
    plt.savefig(path.with_suffix(".svg"))
    plt.close()


def _plot_heatmap(df: pd.DataFrame, cols: list[str], title: str, outfile: Path, vmin=None, vmax=None):
    labels = df["Complex"].tolist()
    vals = df[cols].to_numpy(float)
    plt.figure(figsize=(7.5, max(5, len(labels)*0.18)))
    im = plt.imshow(vals, aspect="auto", vmin=vmin, vmax=vmax)
    plt.colorbar(im, label="Value")
    plt.xticks(range(len(cols)), [c.replace("Fluc_", "").replace("BC_", "") for c in cols], rotation=0)
    step = max(1, len(labels)//35)
    plt.yticks(range(0, len(labels), step), [labels[i] for i in range(0, len(labels), step)], fontsize=6)
    plt.title(title)
    _savefig(outfile)


def figure_s3_from_legacy():
    paths = write_standardized_tables()
    tables = load_legacy_dynamics()
    if not tables:
        print("No legacy dynamic CSV tables found. Run: python call.py --export-rdata")
        return None
    fluc_cols = [f"Fluc_{d}" for d in DOMAINS]
    bc_cols = [f"BC_{d}" for d in DOMAINS]
    for family in ["gprot", "arrestin"]:
        if f"{family}_fluc" in tables:
            _plot_heatmap(tables[f"{family}_fluc"], fluc_cols, f"Figure S3A-style fluctuation fold change ({family})", IMAGES / f"figure_s3a_fluc_fold_change_{family}")
        if f"{family}_bc" in tables:
            _plot_heatmap(tables[f"{family}_bc"], bc_cols, f"Figure S3B-style Bhattacharyya coefficients ({family})", IMAGES / f"figure_s3b_bc_scores_{family}", vmin=0.80, vmax=0.96)
    # S3C MDS from all BC rows, as described in supplementary figure.
    bc_tables = [tables[k] for k in ["gprot_bc", "arrestin_bc"] if k in tables]
    if bc_tables:
        all_bc = pd.concat(bc_tables, ignore_index=True)
        X = all_bc[bc_cols].to_numpy(float)
        # Convert similarity BC to distance. Clip for numerical stability.
        D = 1.0 - np.clip(X @ X.T / (np.linalg.norm(X, axis=1)[:,None] * np.linalg.norm(X, axis=1)[None,:]), 0, 1)
        try:
            coords = MDS(n_components=2, dissimilarity="precomputed", random_state=1, normalized_stress="auto").fit_transform(D)
        except TypeError:
            coords = MDS(n_components=2, dissimilarity="precomputed", random_state=1).fit_transform(D)
        all_bc["PCoord1"] = coords[:,0]
        all_bc["PCoord2"] = coords[:,1]
        all_bc.to_csv(SUMMARY / "dynamics_bc_mds_coordinates.csv", index=False)
        plt.figure(figsize=(8,6))
        for rec, g in all_bc.groupby("DXR"):
            plt.scatter(g.PCoord1, g.PCoord2, s=45, alpha=0.75, label=rec)
            for _, r in g.iterrows():
                plt.text(r.PCoord1, r.PCoord2, r.Complex, fontsize=5, alpha=.7)
        plt.xlabel("PCoord1")
        plt.ylabel("PCoord2")
        plt.title("Figure S3C-style MDS of BC dynamic profiles")
        plt.legend(title="Receptor", fontsize=7)
        plt.grid(alpha=.2)
        _savefig(IMAGES / "figure_s3c_bc_mds")
    return True


if __name__ == "__main__":
    figure_s3_from_legacy()
