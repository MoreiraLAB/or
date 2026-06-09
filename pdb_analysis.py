###################################################################################################
############################ MOREIRA LAB - DATA DRIVEN MOLECULAR DESIGN ###########################
###################################################################################################

from __future__ import annotations
from pathlib import Path
import math
from collections import defaultdict
import numpy as np
import pandas as pd
from Bio.PDB import PDBParser, NeighborSearch
from Bio.PDB.Polypeptide import is_aa
from opioid_metadata import AA1, AA_GROUPS, bw_label, receptor_domain, parse_complex_name, partner_class, load_bw_numbering, load_partner_domain_maps, ensure_dirs, ROOT
from residue_numbering import pdb_to_metadata_residue

HEAVY = {'C','N','O','S','P'}

def atoms_for_residue(res):
    for a in res.get_atoms():
        if a.element in HEAVY and not a.is_disordered():
            yield a

def ca_coord(res):
    if 'CA' in res:
        return np.array(res['CA'].coord, dtype=float)
    atoms = list(atoms_for_residue(res))
    if atoms:
        return np.mean([a.coord for a in atoms], axis=0)
    return None

def residues_by_chain(pdb_path: Path):
    s = PDBParser(QUIET=True).get_structure(pdb_path.stem, str(pdb_path))
    model = next(s.get_models())
    chains = list(model.get_chains())
    return {ch.id: [r for r in ch.get_residues() if is_aa(r, standard=True)] for ch in chains}

def interface_contacts(pdb_path: Path, cutoff: float = 8.0) -> pd.DataFrame:
    """Fast local interface definition: receptor chain A C-alpha/centroid within cutoff of partner residue C-alpha/centroid.
    This replaces the old CoCoMaps/InterProSurf web scraping with a reproducible local calculation.
    """
    receptor, partner, partner_base = parse_complex_name(pdb_path.name)
    chains = residues_by_chain(pdb_path)
    if len(chains) < 2:
        return pd.DataFrame()
    chain_ids = list(chains)
    rec_chain = 'A' if 'A' in chains else chain_ids[0]
    partner_chains = [c for c in chain_ids if c != rec_chain]
    rec_res = chains[rec_chain]
    part_res = [r for c in partner_chains for r in chains[c]]
    part_points=[]
    for pr in part_res:
        xyz=ca_coord(pr)
        if xyz is not None:
            part_points.append((pr, xyz))
    bw_df = load_bw_numbering()
    pdomains = load_partner_domain_maps()
    rows=[]
    for rr in rec_res:
        rxyz=ca_coord(rr)
        if rxyz is None: continue
        hits=[]
        for pr, pxyz in part_points:
            d=float(np.linalg.norm(rxyz-pxyz))
            if d <= cutoff:
                hits.append((d,pr))
        for dist, pr in sorted(hits, key=lambda x:x[0])[:6]:
            rnum = int(rr.id[1]); pnum = int(pr.id[1])
            raa3 = rr.resname; paa3 = pr.resname
            partner_chain = pr.get_parent().id
            pnum_meta = pdb_to_metadata_residue(pdb_path.name, partner_chain, pnum)
            rows.append({
                'complex': pdb_path.stem, 'receptor': receptor, 'partner': partner, 'partner_base': partner_base,
                'partner_class': partner_class(partner), 'receptor_chain': rec_chain, 'partner_chain': partner_chain,
                'receptor_resnum': rnum, 'receptor_resname': raa3, 'receptor_aa': AA1.get(raa3, 'X'),
                'receptor_bw': bw_label(receptor, rnum, bw_df), 'receptor_domain': receptor_domain(receptor, rnum, bw_df),
                'partner_resnum': pnum, 'partner_resnum_metadata': pnum_meta,
                'partner_resname': paa3, 'partner_aa': AA1.get(paa3, 'X'),
                'partner_domain': pdomains.get((partner_base, pnum_meta), pdomains.get((partner, pnum_meta), 'Not defined domain')),
                'distance_angstrom': round(dist,3)
            })
    return pd.DataFrame(rows)

def salt_bridges_hbonds(pdb_path: Path) -> pd.DataFrame:
    # Conservative geometric approximations for automation; detailed chemistry remains auditable from atom-level contacts.
    contacts = interface_contacts(pdb_path, cutoff=8.0)
    if contacts.empty: return contacts
    neg={'ASP','GLU'}; pos={'ARG','LYS','HIS'}; donors={'ARG','LYS','HIS','ASN','GLN','SER','THR','TYR','TRP'}; acceptors={'ASP','GLU','ASN','GLN','SER','THR','TYR','CYS','MET'}
    contacts['salt_bridge_like'] = contacts.apply(lambda r: (r.receptor_resname in neg and r.partner_resname in pos) or (r.receptor_resname in pos and r.partner_resname in neg), axis=1)
    contacts['hbond_like'] = contacts.apply(lambda r: (r.receptor_resname in donors and r.partner_resname in acceptors) or (r.receptor_resname in acceptors and r.partner_resname in donors), axis=1)
    return contacts

def interface_percentages(contacts: pd.DataFrame) -> pd.DataFrame:
    rows=[]
    if contacts.empty: return pd.DataFrame()
    for side, prefix in [('receptor','receptor'), ('partner','partner')]:
        cols = ['complex','receptor','partner','partner_class',f'{prefix}_resnum',f'{prefix}_resname',f'{prefix}_aa']
        unique = contacts[cols].drop_duplicates()
        for (complex_name, rec, partner, pclass), g in unique.groupby(['complex','receptor','partner','partner_class']):
            total=len(g)
            for aa3, gg in g.groupby(f'{prefix}_resname'):
                rows.append({'complex':complex_name,'receptor':rec,'partner':partner,'partner_class':pclass,'side':side,
                             'residue':AA1.get(aa3,'X'),'residue3':aa3,'aa_group':AA_GROUPS.get(aa3,'Other'),
                             'count':len(gg),'total':total,'percentage':100*len(gg)/total if total else 0})
    return pd.DataFrame(rows)

def run_structural(structural_dir: Path = ROOT/'structural_complexes', outdir: Path = ROOT/'processed_results', images_dir: Path = ROOT/'images'):
    ensure_dirs(outdir, images_dir)
    all_contacts=[]
    for pdb in sorted(structural_dir.glob('*.pdb')):
        df = salt_bridges_hbonds(pdb)
        if not df.empty:
            all_contacts.append(df)
    contacts = pd.concat(all_contacts, ignore_index=True) if all_contacts else pd.DataFrame()
    contacts.to_csv(outdir/'interface_contacts.csv', index=False)
    pct = interface_percentages(contacts)
    pct.to_csv(outdir/'interface_percentages.csv', index=False)
    if not contacts.empty:
        summary = contacts.groupby(['complex','receptor','partner','partner_class']).agg(
            contacts=('distance_angstrom','count'),
            salt_bridges_like=('salt_bridge_like','sum'),
            hbonds_like=('hbond_like','sum')
        ).reset_index()
        summary.to_csv(outdir/'interface_summary.csv', index=False)
    return {'contacts': outdir/'interface_contacts.csv', 'percentages': outdir/'interface_percentages.csv'}

if __name__ == '__main__':
    print(run_structural())
