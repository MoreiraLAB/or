###################################################################################################
############################ MOREIRA LAB - DATA DRIVEN MOLECULAR DESIGN ###########################
###################################################################################################

from __future__ import annotations
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from opioid_metadata import ROOT, ensure_dirs
import legacy_dynamics_figures

IMAGES = ROOT/'images'
SUMMARY = ROOT/'summary'
PROCESSED = ROOT/'processed_results'

def savefig(path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(path.with_suffix('.png'), dpi=300)
    plt.savefig(path.with_suffix('.svg'))
    plt.close()

def figure_s1_interhelical(csv_path: Path = SUMMARY/'interhelical_distances.csv'):
    if not csv_path.exists():
        print(f'Missing {csv_path}; run interhelical_bio3d.R or python fallback first.')
        return None
    df = pd.read_csv(csv_path)
    needed={'tm3_tm6','tm3_tm7','partner_class'}
    if not needed.issubset(df.columns):
        print(f'{csv_path} lacks required columns {needed}')
        return None
    for cls,label in [('gprotein','G-protein'),('arrestin','Arrestin')]:
        sub=df[df.partner_class==cls].copy()
        if sub.empty: continue
        plt.figure(figsize=(7,6))
        for rec,g in sub.groupby('receptor'):
            plt.scatter(g.tm3_tm6, g.tm3_tm7, s=90, alpha=0.65, label=rec)
            for _,r in g.iterrows():
                plt.text(r.tm3_tm6, r.tm3_tm7, str(r.partner).replace('ARR','Arr'), fontsize=6, alpha=.8)
        plt.xlabel('TM3-TM6 distance (Å)')
        plt.ylabel('TM3-TM7 distance (Å)')
        plt.title(f'Figure S1-style interhelical distances: {label} complexes')
        plt.grid(alpha=.25)
        plt.legend(title='Receptor')
        savefig(IMAGES/f'figure_s1_interhelical_{cls}')
    return True

def figure_s2_interface_percentages(csv_path: Path = PROCESSED/'interface_percentages.csv'):
    if not csv_path.exists():
        print(f'Missing {csv_path}; run structural analysis first.')
        return None
    df=pd.read_csv(csv_path)
    if df.empty: return None
    order=['A','G','I','L','M','V','F','Y','W','N','C','Q','P','S','T','D','E','R','H','K']
    for side in ['receptor','partner']:
        sub=df[df.side==side]
        if sub.empty: continue
        piv=sub.pivot_table(index='residue', columns='complex', values='percentage', aggfunc='sum').reindex(order).fillna(0)
        # keep most informative complexes if too many; full CSV is retained.
        cols=list(piv.columns)
        if len(cols)>70: cols=cols[:70]
        piv=piv[cols]
        plt.figure(figsize=(max(10, len(cols)*0.18), 6))
        plt.imshow(piv.values, aspect='auto')
        plt.yticks(range(len(piv.index)), piv.index)
        step=max(1, len(cols)//35)
        plt.xticks(range(0,len(cols),step), [cols[i] for i in range(0,len(cols),step)], rotation=90, fontsize=6)
        plt.colorbar(label='Interface percentage')
        plt.title(f'Figure S2-style amino-acid interface percentages: {side} side')
        savefig(IMAGES/f'figure_s2_interface_percentages_{side}')
    return True

def figure2_interaction_signature(csv_path: Path = PROCESSED/'interface_contacts.csv'):
    if not csv_path.exists(): return None
    df=pd.read_csv(csv_path)
    if df.empty: return None
    gp=df[df.partner_class=='gprotein']
    if gp.empty: return None
    gp['subfamily']=gp.partner.str.extract(r'^(Gi|Go|Gob|Gz|Gs|Gq|G11|G12|G14|G15)', expand=False).fillna(gp.partner)
    mat=gp.groupby(['receptor_bw','partner_domain','subfamily']).size().reset_index(name='n')
    pivot=mat.pivot_table(index='receptor_bw', columns='partner_domain', values='n', aggfunc='sum').fillna(0)
    # sort BW labels roughly
    pivot=pivot.loc[sorted(pivot.index, key=lambda x: [int(t) if t.isdigit() else 999 for t in str(x).replace('.',' ').split()[:2]])]
    plt.figure(figsize=(12, max(5, len(pivot)*0.22)))
    plt.imshow(pivot.values>0, aspect='auto')
    plt.yticks(range(len(pivot.index)), pivot.index, fontsize=6)
    plt.xticks(range(len(pivot.columns)), pivot.columns, rotation=90, fontsize=7)
    plt.title('Figure 2-style receptor BW vs partner-domain interaction signature')
    plt.colorbar(label='Contact present')
    savefig(IMAGES/'figure2_interaction_signature')
    return True

def figure_s3_from_dynamic_exports(dyn_dir: Path = PROCESSED/'dynamics'):
    """Generate the paper-style Figure S3 panels from exported legacy .RData tables.

    Expected RData objects:
      - gprot_bc_scores / arrestin_bc_scores: n_complexes x 8 matrix
      - gprot_fluc_fold_change / arrestin_fluc_fold_change: data.frame with
        Fluc_TM1..Fluc_H8, Complex, DXR, Partner

    The 8 domains are fixed in the order used by the original pipeline:
    TM1, TM2, TM3, TM4, TM5, TM6, TM7, H8.
    """
    return legacy_dynamics_figures.figure_s3_from_legacy()

def all_figures():
    ensure_dirs(IMAGES)
    return {
        's1': figure_s1_interhelical(),
        's2': figure_s2_interface_percentages(),
        'fig2': figure2_interaction_signature(),
        's3': figure_s3_from_dynamic_exports(),
    }

if __name__ == '__main__':
    print(all_figures())
