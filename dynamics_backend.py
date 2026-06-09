###################################################################################################
############################ MOREIRA LAB - DATA DRIVEN MOLECULAR DESIGN ###########################
###################################################################################################

from __future__ import annotations
import argparse
import shutil
import subprocess
from pathlib import Path
from typing import Dict

from opioid_metadata import ROOT, ensure_dirs
import dynamics_from_scratch
import legacy_dynamics_figures


def run_cmd(cmd: list[str | Path]) -> int:
    print('+', ' '.join(map(str, cmd)))
    return subprocess.call([str(c) for c in cmd])


def run_bio3d(limit: int | None = None) -> int:
    if not shutil.which('Rscript'):
        print('Rscript not found. Install the conda environment or use --backend prody.')
        return 127
    cmd = ['Rscript', ROOT / 'dynamics_bio3d.R', '--root', ROOT]
    if limit:
        cmd += ['--limit', str(limit)]
    return run_cmd(cmd)


def run_prody_like(cutoff: float = 10.0, limit: int | None = None) -> Dict[str, Path]:
    return dynamics_from_scratch.run_all(cutoff=cutoff, limit=limit)


def main() -> None:
    ap = argparse.ArgumentParser(description='Run dynamics backend: Bio3D reference, Python/ProDy-like fallback, or legacy RData plotting.')
    ap.add_argument('--backend', choices=['bio3d', 'prody', 'both', 'legacy'], default='bio3d')
    ap.add_argument('--cutoff', type=float, default=10.0)
    ap.add_argument('--limit', type=int, default=None)
    args = ap.parse_args()
    ensure_dirs(ROOT/'processed_results', ROOT/'summary', ROOT/'images')

    if args.backend in {'bio3d', 'both'}:
        code = run_bio3d(limit=args.limit)
        if code != 0 and args.backend == 'bio3d':
            raise SystemExit(code)
    if args.backend in {'prody', 'both'}:
        run_prody_like(cutoff=args.cutoff, limit=args.limit)
    if args.backend == 'legacy':
        legacy_dynamics_figures.figure_s3_from_legacy()


if __name__ == '__main__':
    main()
