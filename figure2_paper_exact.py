###################################################################################################
############################ MOREIRA LAB - DATA DRIVEN MOLECULAR DESIGN ###########################
###################################################################################################

#!/usr/bin/env python3
from __future__ import annotations
from pathlib import Path
import argparse, re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm

ROOT = Path(__file__).resolve().parent
RECEPTORS = ['DOR','KOR','MOR','NOP']
SUBFAMILY_ORDER = ['Gi/o','Gs','Gq/11','G12/13']
SUBFAMILY_CODE = {'Gi/o':1, 'Gs':2, 'Gq/11':3, 'G12/13':4}


def partner_subfamily(partner: str) -> str | None:
    p = str(partner)
    if re.search(r'Gs', p, re.I):
        return 'Gs'
    if re.search(r'Gq|G11', p, re.I):
        return 'Gq/11'
    if re.search(r'G12|G13', p, re.I):
        return 'G12/13'
    if re.search(r'Gi|Go|Gob|Gz', p, re.I):
        return 'Gi/o'
    return None


def _find_contacts() -> Path:
    candidates = [
        'processed_results/interface_contacts.csv',
        'processed_results/contacts.csv',
        'summary/interface_contacts.csv',
        'processed_results/receptor_partner_contacts.csv',
    ]
    for c in candidates:
        p = ROOT / c
        if p.exists():
            return p
    raise FileNotFoundError('No contact/interface table found. Run python call.py first.')


def _first_col(df, names, required=True):
    lower = {c.lower(): c for c in df.columns}
    for n in names:
        if n.lower() in lower:
            return lower[n.lower()]
    if required:
        raise KeyError(f'Missing any of columns {names}. Available: {list(df.columns)}')
    return None


def load_contacts() -> pd.DataFrame:
    df = pd.read_csv(_find_contacts())
    complex_col = _first_col(df, ['complex','complex_name','model'], required=False)
    rec_col = _first_col(df, ['receptor','dxr','OR'], required=False)
    partner_col = _first_col(df, ['partner','gprotein','g_protein','arrestin'], required=False)
    bw_col = _first_col(df, ['bw','receptor_bw','bw_label','receptor_bw_label','receptor_position_bw'], required=False)
    partner_domain_col = _first_col(df, ['partner_domain','domain_partner','partner_region','gprotein_domain','g_protein_domain','arrestin_domain'], required=False)
    partner_bw_col = _first_col(df, ['partner_bw','partner_label','partner_position','partner_residue_label','partner_resnum'], required=False)

    if complex_col and (not rec_col or not partner_col):
        tmp = df[complex_col].astype(str).str.extract(r'^(DOR|KOR|MOR|NOP)-(.+)$')
        if not rec_col:
            df['receptor'] = tmp[0]
            rec_col = 'receptor'
        if not partner_col:
            df['partner'] = tmp[1]
            partner_col = 'partner'

    if not bw_col:
        # fallback to residue number if BW absent
        bw_col = _first_col(df, ['receptor_resnum','resnum','position'], required=True)
    if not partner_domain_col:
        partner_domain_col = partner_bw_col or _first_col(df, ['partner_resnum','partner_residue','partner_aa'], required=True)

    out = pd.DataFrame({
        'receptor': df[rec_col].astype(str),
        'partner': df[partner_col].astype(str),
        'receptor_bw': df[bw_col].astype(str),
        'partner_axis': df[partner_domain_col].astype(str),
    })
    out['subfamily'] = out['partner'].map(partner_subfamily)
    out = out[out['receptor'].isin(RECEPTORS) & out['subfamily'].notna()].drop_duplicates()
    return out


def common_interactions(df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for sf, g in df.groupby('subfamily'):
        for pair, h in g.groupby(['receptor_bw','partner_axis']):
            recs = set(h['receptor'])
            if all(r in recs for r in RECEPTORS):
                rows.append({'receptor_bw': pair[0], 'partner_axis': pair[1], 'subfamily': sf, 'code': SUBFAMILY_CODE[sf]})
    return pd.DataFrame(rows)


def _sort_bw(vals):
    def key(x):
        s = str(x).replace('H8.','8.').replace('ICL','')
        m = re.search(r'(\d+)(?:\.(\d+))?', s)
        if m:
            return (int(m.group(1)), int(m.group(2) or 0), str(x))
        return (999, 999, str(x))
    return sorted(vals, key=key)


def _sort_partner(vals):
    order = ['HN','hns1','S3','h4s6','H5','S1','S3','S6','finger','C-loop','lariat','bottom','middle']
    def key(x):
        sx = str(x)
        for i,o in enumerate(order):
            if o.lower() in sx.lower():
                return (i, sx)
        return (99, sx)
    return sorted(vals, key=key)


def plot_figure2(common: pd.DataFrame) -> Path:
    outdir = ROOT/'images'; outdir.mkdir(exist_ok=True)
    (ROOT/'summary').mkdir(exist_ok=True)
    common.to_csv(ROOT/'summary/figure2_common_interactions.csv', index=False)
    if common.empty:
        raise ValueError('No common interactions across DOR/KOR/MOR/NOP were found. Check contact table/BW/domain columns.')
    yvals = _sort_bw(common['receptor_bw'].unique())
    xvals = _sort_partner(common['partner_axis'].unique())
    mat = np.zeros((len(yvals), len(xvals)), dtype=int)
    ymap = {v:i for i,v in enumerate(yvals)}
    xmap = {v:i for i,v in enumerate(xvals)}
    for _, r in common.iterrows():
        mat[ymap[r['receptor_bw']], xmap[r['partner_axis']]] = int(r['code'])

    cmap = ListedColormap(['white', '#e67e22', '#f1c40f', '#7cb342', '#3498db'])
    norm = BoundaryNorm([-0.5,0.5,1.5,2.5,3.5,4.5], cmap.N)
    plt.figure(figsize=(max(10, 0.28*len(xvals)), max(7, 0.18*len(yvals))))
    im = plt.imshow(mat, aspect='auto', cmap=cmap, norm=norm)
    plt.xticks(range(len(xvals)), xvals, rotation=90, fontsize=6)
    plt.yticks(range(len(yvals)), yvals, fontsize=6)
    plt.title('Figure 2-style common OR–G-protein interaction signature')
    plt.xlabel('Partner domain / residue label')
    plt.ylabel('Receptor BW / residue label')
    handles = [plt.Line2D([0],[0], marker='s', linestyle='', markersize=10, markerfacecolor=cmap(SUBFAMILY_CODE[sf]), markeredgecolor='black', label=sf) for sf in SUBFAMILY_ORDER]
    plt.legend(handles=handles, title='G-protein subfamily', bbox_to_anchor=(1.02,1), loc='upper left')
    plt.tight_layout()
    png = outdir/'figure2_paper_exact.png'
    svg = outdir/'figure2_paper_exact.svg'
    plt.savefig(png, dpi=300)
    plt.savefig(svg)
    plt.close()
    return png


def run_all():
    df = load_contacts()
    common = common_interactions(df)
    return plot_figure2(common)


def main():
    ap = argparse.ArgumentParser()
    args = ap.parse_args()
    print(run_all())

if __name__ == '__main__':
    main()
