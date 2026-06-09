###################################################################################################
############################ MOREIRA LAB - DATA DRIVEN MOLECULAR DESIGN ###########################
###################################################################################################

from __future__ import annotations
from pathlib import Path
import re, argparse, shutil
import pandas as pd
import matplotlib.pyplot as plt
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa
from opioid_metadata import ROOT, AA1, bw_label, receptor_domain, load_bw_numbering, parse_complex_name, ensure_dirs

CONSURF_URL = 'https://consurf.tau.ac.il/'


def pdb_to_fasta(pdb_path: Path, chain: str = 'A') -> str:
    s = PDBParser(QUIET=True).get_structure(pdb_path.stem, str(pdb_path))
    residues = []
    for ch in s[0]:
        if ch.id == chain:
            for r in ch.get_residues():
                if is_aa(r, standard=True):
                    residues.append(AA1.get(r.resname, 'X'))
    return ''.join(residues)


def export_fasta(pdb_dir=ROOT/'structural_complexes', outdir=ROOT/'consurf/input', chain: str = 'A'):
    ensure_dirs(Path(outdir)); rows = []
    for pdb in sorted(Path(pdb_dir).glob('*.pdb')):
        rec, partner, _ = parse_complex_name(pdb.name)
        seq = pdb_to_fasta(pdb, chain)
        if not seq:
            continue
        f = Path(outdir) / f'{pdb.stem}_chain{chain}.fasta'
        f.write_text(f'>{pdb.stem}|chain {chain}|{rec}\n{seq}\n')
        rows.append({'complex': pdb.stem, 'receptor': rec, 'chain': chain, 'fasta': str(f), 'length': len(seq)})
    manifest = Path(outdir) / 'consurf_fasta_manifest.csv'
    pd.DataFrame(rows).to_csv(manifest, index=False)
    return rows


def export_pdb(pdb_dir=ROOT/'structural_complexes', outdir=ROOT/'consurf/input_pdb', chain: str = 'A', all_complexes: bool = False):
    """Prepare PDB files for ConSurf current server.

    By default this exports only one representative structure per opioid receptor
    (DOR, KOR, MOR and NOP), because ConSurf conservation is receptor-sequence
    based and therefore does not need to be recalculated for every receptor-partner
    complex. The resulting receptor-level scores are later mapped back to all
    interface/contact tables by receptor and residue position.

    Set all_complexes=True only if you explicitly want one ConSurf job for every
    structural complex.
    """
    pdb_dir = Path(pdb_dir)
    outdir = Path(outdir)
    ensure_dirs(outdir)

    all_pdbs = sorted(pdb_dir.glob('*.pdb'))
    if not all_complexes:
        selected = {}
        receptor_order = ['DOR', 'KOR', 'MOR', 'NOP']
        for rec in receptor_order:
            for pdb in all_pdbs:
                r, partner, _ = parse_complex_name(pdb.name)
                if r == rec:
                    selected[rec] = pdb
                    break
        pdbs = [selected[r] for r in receptor_order if r in selected]
    else:
        pdbs = all_pdbs

    rows = []
    for pdb in pdbs:
        rec, partner, _ = parse_complex_name(pdb.name)
        dest = outdir / pdb.name
        shutil.copy2(pdb, dest)
        rows.append({
            'complex': pdb.stem,
            'receptor': rec,
            'partner': partner,
            'chain': chain,
            'pdb': str(dest),
            'structure_file': str(dest),
            'job_name': pdb.stem,
            'conservation_scope': 'all_complexes' if all_complexes else 'receptor_representative',
        })

    manifest = outdir / 'consurf_pdb_manifest.csv'
    pd.DataFrame(rows).to_csv(manifest, index=False)
    ensure_dirs(ROOT/'consurf/input')
    pd.DataFrame(rows).to_csv(ROOT/'consurf/input/consurf_pdb_manifest.csv', index=False)
    print(f'ConSurf PDB manifest written with {len(rows)} job(s). all_complexes={all_complexes}')
    return rows

def parse_consurf_file(path: Path, receptor: str | None = None) -> pd.DataFrame:
    text = path.read_text(errors='ignore')
    rows = []
    for line in text.splitlines():
        if not line.strip() or line.lstrip().startswith(('#', '%', ';')):
            continue
        toks = re.split(r'\s+', line.strip())
        if len(toks) >= 3 and toks[0].isdigit():
            pos = int(toks[0]); aa = toks[1]
            grade = None; score = None
            for t in toks[2:]:
                if re.fullmatch(r'[1-9]', t):
                    grade = int(t)
                elif re.fullmatch(r'-?\d+\.\d+', t):
                    score = float(t)
                    break
            rows.append({'position': pos, 'aa': aa, 'score': score, 'grade': grade})
    df = pd.DataFrame(rows)
    if not df.empty and receptor:
        bw = load_bw_numbering()
        df['receptor'] = receptor
        df['bw'] = df.position.apply(lambda x: bw_label(receptor, int(x), bw))
        df['domain'] = df.position.apply(lambda x: receptor_domain(receptor, int(x), bw))
    return df


def collect_results(input_dir=ROOT/'consurf/results', out=ROOT/'processed_results/consurf_scores.csv'):
    ensure_dirs(Path(out).parent)
    frames = []
    for p in Path(input_dir).rglob('*'):
        if p.is_file() and p.suffix.lower() in ['.txt', '.csv', '.grades', '.dat']:
            rec = None
            for r in ['DOR', 'KOR', 'MOR', 'NOP']:
                if r in p.name.upper() or r in str(p.parent).upper():
                    rec = r
            df = parse_consurf_file(p, rec)
            if not df.empty:
                df['source_file'] = str(p)
                frames.append(df)
    all_df = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame(columns=['position','aa','score','grade','receptor','bw','domain','source_file'])
    all_df.to_csv(out, index=False)
    return all_df


def plot_consurf(scores=ROOT/'processed_results/consurf_scores.csv'):
    scores = Path(scores)
    if not scores.exists():
        return None
    df = pd.read_csv(scores)
    if df.empty:
        return None
    ensure_dirs(ROOT/'images')
    for rec, g in df.groupby('receptor'):
        plt.figure(figsize=(12, 3))
        y = g['grade'] if 'grade' in g and g['grade'].notna().any() else g['score']
        plt.scatter(g['position'], y, s=12)
        plt.xlabel('Residue position')
        plt.ylabel('ConSurf grade/score')
        plt.title(f'ConSurf conservation profile: {rec}')
        plt.tight_layout()
        plt.savefig(ROOT/'images'/f'consurf_profile_{rec}.png', dpi=300)
        plt.savefig(ROOT/'images'/f'consurf_profile_{rec}.svg')
        plt.close()
    contacts = ROOT/'processed_results/interface_contacts.csv'
    if contacts.exists():
        c = pd.read_csv(contacts)
        if 'receptor_resnum' in c:
            merged = c.merge(df, left_on=['receptor','receptor_resnum'], right_on=['receptor','position'], how='left')
            merged.to_csv(ROOT/'processed_results/interface_contacts_with_consurf.csv', index=False)
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--export-fasta', action='store_true')
    ap.add_argument('--export-pdb', action='store_true')
    ap.add_argument('--chain', default='A')
    ap.add_argument('--all-complexes', action='store_true', help='Export every complex instead of one representative per receptor.')
    ap.add_argument('--collect', action='store_true')
    ap.add_argument('--plot', action='store_true')
    args = ap.parse_args()
    if args.export_fasta:
        print(export_fasta(chain=args.chain))
    if args.export_pdb:
        print(export_pdb(chain=args.chain, all_complexes=args.all_complexes))
    if args.collect:
        print(collect_results())
    if args.plot:
        print(plot_consurf())

if __name__ == '__main__':
    main()
