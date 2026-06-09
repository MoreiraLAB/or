#!/usr/bin/env python3
from __future__ import annotations
from pathlib import Path
import argparse, subprocess, shutil, textwrap

ROOT = Path(__file__).resolve().parent
DEFAULT_COMPLEXES = {
    'gio': ['DOR-Gi1.pdb', 'DOR-Gi1_6DDF.pdb'],
    'gs': ['KOR-Gssh.pdb', 'KOR-Gssh_3SN6.pdb', 'KOR-Gssh_6DDF.pdb'],
    'gq': ['DOR-Gq_6DDF.pdb', 'DOR-Gq.pdb'],
    'g12': ['KOR-G12_6DDF.pdb', 'KOR-G12.pdb'],
}

PYMOL_SCRIPT_TEMPLATE = r'''
reinitialize
load {pdb}, complex
hide everything
bg_color white
set ray_opaque_background, off
set antialias, 2
set orthoscopic, on
set ray_trace_mode, 1

# Generic chain-based split. Adjust chain names manually if needed.
select receptor, chain A
select partner, not receptor

show surface, receptor or partner
color grey80, receptor or partner
set transparency, 0.25, receptor or partner

# Receptor colours approximate Figure 3 legend. These selections use residue ranges as fallback;
# if BW/domain metadata are available, replace these by exact selections in PyMOL.
color yellow, receptor and resi 55-75
color palegreen, receptor and resi 110-140
color forest, receptor and resi 130-150
color limegreen, receptor and resi 150-175
color lightblue, receptor and resi 190-230
color marine, receptor and resi 230-270
color blue, receptor and resi 270-310
color purple, receptor and resi 310-340
color red, receptor and resi 330-360

# Partner broad highlights; domains differ between G proteins and arrestins.
color slate, partner
color grey40, partner and resi 330-380
color cyan, partner and resi 250-330
color yelloworange, partner and resi 180-250

select interface, byres ((receptor within 4.5 of partner) or (partner within 4.5 of receptor))
show sticks, interface
color red, interface and receptor
color black, interface and partner

orient complex
ray 2200, 1600
png {out_png}, dpi=300
save {out_pse}
quit
'''


def find_complex(candidates):
    dirs = [ROOT/'structural_complexes', ROOT]
    for d in dirs:
        for c in candidates:
            p = d / c
            if p.exists():
                return p
    return None


def run_pymol(pdb: Path, label: str) -> Path:
    outdir = ROOT/'images'/'figure3_panels'
    outdir.mkdir(parents=True, exist_ok=True)
    script = outdir / f'{label}.pml'
    png = outdir / f'figure3_panel_{label}.png'
    pse = outdir / f'figure3_panel_{label}.pse'
    script.write_text(PYMOL_SCRIPT_TEMPLATE.format(pdb=str(pdb), out_png=str(png), out_pse=str(pse)))
    exe = shutil.which('pymol') or shutil.which('pymol-open-source')
    if not exe:
        raise SystemExit('PyMOL not found. Install with: conda install -c conda-forge pymol-open-source')
    subprocess.check_call([exe, '-cq', str(script)])
    return png


def montage(panels):
    try:
        from PIL import Image, ImageOps, ImageDraw
    except Exception:
        print('Pillow not installed; panels were generated but montage skipped.')
        return None
    imgs = []
    for label, p in panels:
        im = Image.open(p).convert('RGB')
        im.thumbnail((900, 650))
        canvas = Image.new('RGB', (900, 700), 'white')
        canvas.paste(im, ((900-im.width)//2, 40))
        ImageDraw.Draw(canvas).text((20,20), label.upper(), fill='black')
        imgs.append(canvas)
    if not imgs:
        return None
    w = sum(i.width for i in imgs)
    h = max(i.height for i in imgs)
    out = Image.new('RGB', (w,h), 'white')
    x = 0
    for im in imgs:
        out.paste(im, (x,0)); x += im.width
    outpath = ROOT/'images/figure3_paper_style.png'
    out.save(outpath, dpi=(300,300))
    return outpath


def run_all():
    panels = []
    for label, candidates in DEFAULT_COMPLEXES.items():
        pdb = find_complex(candidates)
        if not pdb:
            print(f'Skipping {label}: no PDB found among {candidates}')
            continue
        panels.append((label, run_pymol(pdb, label)))
    return montage(panels)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--label', choices=list(DEFAULT_COMPLEXES)+['all'], default='all')
    args = ap.parse_args()
    if args.label == 'all':
        print(run_all())
    else:
        pdb = find_complex(DEFAULT_COMPLEXES[args.label])
        if not pdb:
            raise SystemExit(f'No PDB found for {args.label}')
        print(run_pymol(pdb, args.label))

if __name__ == '__main__':
    main()
