###################################################################################################
############################ MOREIRA LAB - DATA DRIVEN MOLECULAR DESIGN ###########################
###################################################################################################

from __future__ import annotations
from pathlib import Path
import re
import pandas as pd

ROOT = Path(__file__).resolve().parent
METADATA = ROOT / 'metadata'

AA_GROUPS = {
    'ALA':'Nonpolar aliphatic','GLY':'Nonpolar aliphatic','ILE':'Nonpolar aliphatic','LEU':'Nonpolar aliphatic','MET':'Nonpolar aliphatic','VAL':'Nonpolar aliphatic',
    'PHE':'Nonpolar aromatic','TYR':'Nonpolar aromatic','TRP':'Nonpolar aromatic',
    'ASN':'Polar uncharged','CYS':'Polar uncharged','GLN':'Polar uncharged','PRO':'Polar uncharged','SER':'Polar uncharged','THR':'Polar uncharged',
    'ASP':'Acid negative','GLU':'Acid negative',
    'ARG':'Basic positive','HIS':'Basic positive','LYS':'Basic positive',
}
AA1 = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}


def read_semicolon_csv(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep=';', dtype=str, keep_default_na=False, encoding='utf-8-sig')


def load_bw_numbering(path: Path | None = None) -> pd.DataFrame:
    path = path or METADATA / 'weinstein_numbering_opioids.csv'
    df = read_semicolon_csv(path)
    for c in ['start','end','x']:
        df[c] = pd.to_numeric(df[c], errors='coerce').astype('Int64')
    return df


def bw_label(receptor: str, residue_number: int, bw_df: pd.DataFrame | None = None) -> str:
    bw_df = bw_df if bw_df is not None else load_bw_numbering()
    sub = bw_df[(bw_df['receptor'].str.upper() == receptor.upper())]
    for _, r in sub.iterrows():
        start, end, x = int(r.start), int(r.end), int(r.x)
        if start <= residue_number <= end:
            decimal = 50 + (residue_number - x)
            return f"{str(r.TM).replace('TM','')}.{decimal}"
    return str(residue_number)


def receptor_domain(receptor: str, residue_number: int, bw_df: pd.DataFrame | None = None) -> str:
    bw_df = bw_df if bw_df is not None else load_bw_numbering()
    sub = bw_df[(bw_df['receptor'].str.upper() == receptor.upper())]
    for _, r in sub.iterrows():
        if int(r.start) <= residue_number <= int(r.end):
            return str(r.TM)
    # approximate loop assignment between TMs
    ranges = [(str(r.TM), int(r.start), int(r.end)) for _, r in sub.iterrows()]
    ranges = sorted(ranges, key=lambda x: x[1])
    for (tm_a, s_a, e_a), (tm_b, s_b, e_b) in zip(ranges, ranges[1:]):
        if e_a < residue_number < s_b:
            if tm_a == 'TM3' and tm_b == 'TM4': return 'ICL2'
            if tm_a == 'TM5' and tm_b == 'TM6': return 'ICL3'
            if tm_a == 'TM1' and tm_b == 'TM2': return 'ICL1'
            return f'{tm_a}-{tm_b} loop'
    if ranges and residue_number > ranges[-1][2]:
        return 'H8'
    return 'Not defined domain'


def parse_partner_domains(path: Path) -> dict[tuple[str,int], str]:
    df = read_semicolon_csv(path)
    mapping = {}
    for _, row in df.iterrows():
        partner = row.iloc[0]
        domain = row.iloc[1]
        for val in row.iloc[2:].tolist():
            if str(val).strip().isdigit():
                mapping[(str(partner), int(val))] = str(domain)
    return mapping


def load_partner_domain_maps() -> dict[tuple[str,int], str]:
    m = {}
    for fname in ['g_proteins_final_opioids.csv', 'arrestins_final_opioids.csv']:
        p = METADATA / fname
        if p.exists():
            m.update(parse_partner_domains(p))
    return m


def parse_complex_name(name: str) -> tuple[str, str, str]:
    stem = Path(name).stem
    m = re.match(r'^(DOR|KOR|MOR|NOP)[-_](.+)$', stem, flags=re.I)
    if not m:
        return stem.split('-')[0].upper(), '', ''
    receptor = m.group(1).upper()
    partner = m.group(2)
    partner_base = partner.split('_')[0]
    return receptor, partner, partner_base


def partner_class(partner: str) -> str:
    p = partner.lower()
    if p.startswith('arr'):
        return 'arrestin'
    if p.startswith('g'):
        return 'gprotein'
    return 'other'


def ensure_dirs(*paths: Path):
    for p in paths:
        p.mkdir(parents=True, exist_ok=True)
