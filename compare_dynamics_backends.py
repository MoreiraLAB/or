###################################################################################################
############################ MOREIRA LAB - DATA DRIVEN MOLECULAR DESIGN ###########################
###################################################################################################

from __future__ import annotations
from pathlib import Path
import argparse
import pandas as pd
import numpy as np

from opioid_metadata import ROOT, ensure_dirs

DOMAINS = ['TM1','TM2','TM3','TM4','TM5','TM6','TM7','H8']

def read_domain(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    df = pd.read_csv(path)
    return df[df['domain'].isin(DOMAINS)].copy()

def compare(a: pd.DataFrame, b: pd.DataFrame, label_a: str, label_b: str) -> pd.DataFrame:
    keys = ['complex','domain']
    cols = ['mean_fold_change','bhattacharyya_coefficient']
    m = a[keys+cols].merge(b[keys+cols], on=keys, suffixes=(f'_{label_a}', f'_{label_b}'))
    rows = []
    for metric in cols:
        x = m[f'{metric}_{label_a}'].astype(float)
        y = m[f'{metric}_{label_b}'].astype(float)
        valid = np.isfinite(x) & np.isfinite(y)
        if valid.sum() < 3:
            corr = np.nan; rmse = np.nan
        else:
            corr = float(np.corrcoef(x[valid], y[valid])[0,1])
            rmse = float(np.sqrt(np.mean((x[valid]-y[valid])**2)))
        rows.append({'metric': metric, 'n': int(valid.sum()), 'pearson_r': corr, 'rmse': rmse})
    return pd.DataFrame(rows)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--bio3d', default='processed_results/dynamics_bio3d/bio3d_domain_dynamics_metrics.csv')
    ap.add_argument('--python', default='processed_results/dynamics_from_scratch/domain_dynamics_metrics.csv')
    ap.add_argument('--out', default='summary/dynamics_backend_comparison.csv')
    args = ap.parse_args()
    ensure_dirs(ROOT/'summary')
    bio = read_domain(ROOT/args.bio3d)
    py = read_domain(ROOT/args.python)
    out = compare(bio, py, 'bio3d', 'python')
    out_path = ROOT/args.out
    out.to_csv(out_path, index=False)
    print(out.to_string(index=False))
    print('Wrote', out_path)

if __name__ == '__main__':
    main()
