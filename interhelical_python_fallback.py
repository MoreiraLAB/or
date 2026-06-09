###################################################################################################
############################ MOREIRA LAB - DATA DRIVEN MOLECULAR DESIGN ###########################
###################################################################################################

from __future__ import annotations
from pathlib import Path
import numpy as np, pandas as pd
from pdb_analysis import residues_by_chain, ca_coord
from opioid_metadata import ROOT, load_bw_numbering, parse_complex_name, partner_class, ensure_dirs

def anchor(bw, rec, tm, dec):
    row=bw[(bw.receptor.str.upper()==rec.upper()) & (bw.TM==tm)]
    if row.empty: return None
    return int(row.iloc[0]['x']) + (dec-50)

def dist(a,b):
    if a is None or b is None: return np.nan
    return float(np.linalg.norm(a-b))

def angle(a,b,c,d):
    if any(x is None for x in [a,b,c,d]): return np.nan
    v1=b-a; v2=d-c
    den=np.linalg.norm(v1)*np.linalg.norm(v2)
    return float(np.degrees(np.arccos(np.clip(np.dot(v1,v2)/den,-1,1)))) if den else np.nan

def run(pdb_dir=ROOT/'structural_complexes', out=ROOT/'summary/interhelical_distances.csv'):
    bw=load_bw_numbering(); rows=[]
    for pdb in sorted(Path(pdb_dir).glob('*.pdb')):
        rec, partner, pbase=parse_complex_name(pdb.name)
        if rec not in set(bw.receptor): continue
        ch=residues_by_chain(pdb).get('A') or []
        res={int(r.id[1]): r for r in ch}
        r350=anchor(bw,rec,'TM3',50); r630=anchor(bw,rec,'TM6',30); r753=anchor(bw,rec,'TM7',53); r847=anchor(bw,rec,'H8',47)
        def xyz(n): return ca_coord(res[n]) if n in res else None
        tm3=bw[(bw.receptor==rec)&(bw.TM=='TM3')].iloc[0]; tm6=bw[(bw.receptor==rec)&(bw.TM=='TM6')].iloc[0]; tm7=bw[(bw.receptor==rec)&(bw.TM=='TM7')].iloc[0]
        rows.append(dict(complex=pdb.stem,receptor=rec,partner=partner,partner_class=partner_class(partner),tm3_anchor=r350,tm6_anchor=r630,tm7_anchor=r753,h8_anchor=r847,tm3_tm6=dist(xyz(r350),xyz(r630)),tm3_tm7=dist(xyz(r350),xyz(r753)),tm7_h8=dist(xyz(r753),xyz(r847)),tm3_tm6_angle=angle(xyz(int(tm3.start)),xyz(int(tm3.end)),xyz(int(tm6.start)),xyz(int(tm6.end))),tm3_tm7_angle=angle(xyz(int(tm3.start)),xyz(int(tm3.end)),xyz(int(tm7.start)),xyz(int(tm7.end)))))
    df=pd.DataFrame(rows); ensure_dirs(Path(out).parent); df.to_csv(out,index=False); return out
if __name__=='__main__': print(run())
