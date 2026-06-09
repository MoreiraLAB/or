###################################################################################################
############################ MOREIRA LAB - DATA DRIVEN MOLECULAR DESIGN ###########################
###################################################################################################

from __future__ import annotations
import argparse, subprocess, shutil
from pathlib import Path
from opioid_metadata import ROOT, ensure_dirs
import pdb_analysis
import plot_paper_figures
import figure_s2_groups
import figure2_paper_exact
import consurf_pipeline
import consurf_auto
import consurf_results_parser
import dynamics_from_scratch
import dynamics_backend


def run(cmd):
    print('+', ' '.join(map(str, cmd)))
    return subprocess.call(list(map(str, cmd)))

def main():
    p=argparse.ArgumentParser(description='Rebuilt OR pipeline: Python-first, Bio3D only for interhelical/dynamics legacy where requested.')
    p.add_argument('--skip-structural', action='store_true')
    p.add_argument('--skip-interhelical', action='store_true')
    p.add_argument('--use-r-bio3d', action='store_true', help='Use interhelical_bio3d.R when Rscript+bio3d are available; otherwise Python fallback is used.')
    p.add_argument('--export-rdata', action='store_true')
    p.add_argument('--skip-figures', action='store_true')
    p.add_argument('--paper-exact-fig2', action='store_true', help='Generate Figure 2-style common OR-G-protein interaction signature.')
    p.add_argument('--paper-s2-groups', action='store_true', help='Generate Figure S2 chemical-group interface percentage heatmaps.')
    p.add_argument('--paper-figure3', action='store_true', help='Generate Figure 3-style PyMOL structural panels if PyMOL is installed.')
    p.add_argument('--paper-extra-figures', action='store_true', help='Generate all extra paper-matching figures: exact Fig2, S2 groups and Figure3.')
    p.add_argument('--consurf-export-fasta', action='store_true', help='Export FASTA for manual sequence-mode ConSurf use.')
    p.add_argument('--consurf-export-pdb', action='store_true', help='Create PDB manifest for current ConSurf structure-mode submission.')
    p.add_argument('--consurf-all-complexes', action='store_true', help='ConSurf: export/submit every complex. Default is only one representative structure per receptor, which is usually sufficient.')
    p.add_argument('--consurf-submit', action='store_true', help='Open/fill current ConSurf form using PDB files, job name, chain and email.')
    p.add_argument('--consurf-email')
    p.add_argument('--consurf-chain', default='A')
    p.add_argument('--consurf-limit', type=int, default=None)
    p.add_argument('--consurf-click-submit', action='store_true', help='Actually click Submit on ConSurf. Omit for browser-review mode.')
    p.add_argument('--consurf-collect', action='store_true')
    p.add_argument('--consurf-plot', action='store_true')
    p.add_argument('--consurf-fetch', action='store_true', help='Legacy alias: watch/download recorded ConSurf jobs from submission_log.json.')
    p.add_argument('--consurf-watch-submission-log', action='store_true', help='Poll submitted ConSurf jobs until completed, then download result files automatically.')
    p.add_argument('--consurf-watch-job', help='Poll one ConSurf job number or progress/results URL until completed, then download result files.')
    p.add_argument('--consurf-wait-minutes', type=int, default=180)
    p.add_argument('--consurf-poll-seconds', type=int, default=60)
    p.add_argument('--consurf-download-url', help='Completed ConSurf result URL to download all linked result files from.')
    p.add_argument('--consurf-job-name', default='manual_job', help='Folder name used with --consurf-download-url.')
    p.add_argument('--dynamics-from-scratch', action='store_true', help='Compute NMA-like CA elastic-network dynamics directly in Python and make Figure S3 outputs.')
    p.add_argument('--dynamics-backend', choices=['bio3d','prody','both','legacy'], default=None, help='Recalculate dynamics with Bio3D reference backend, Python/ProDy-like backend, both, or legacy RData plotting.')
    p.add_argument('--dynamics-cutoff', type=float, default=10.0, help='CA cutoff in Å for Python elastic-network dynamics.')
    p.add_argument('--inspect-rdata', action='store_true', help='Inspect legacy .RData objects and write object summaries/str output.')
    args=p.parse_args()
    ensure_dirs(ROOT/'processed_results', ROOT/'summary', ROOT/'images')

    if not args.skip_structural:
        pdb_analysis.run_structural()

    if not args.skip_interhelical:
        if args.use_r_bio3d and shutil.which('Rscript'):
            code=run(['Rscript', ROOT/'interhelical_bio3d.R'])
            if code != 0:
                print('R/Bio3D interhelical step failed; using Python fallback.')
                run(['python', ROOT/'interhelical_python_fallback.py'])
        else:
            run(['python', ROOT/'interhelical_python_fallback.py'])

    if args.inspect_rdata:
        if shutil.which('Rscript'):
            for rdata in [ROOT/'legacy_rdata/gprot-dynamics.RData', ROOT/'legacy_rdata/arrestin-dynamics.RData']:
                run(['Rscript', ROOT/'scripts/inspect_rdata.R', rdata, ROOT/'processed_results/rdata_inspection'])
        else:
            print('Rscript not available; cannot inspect .RData. Use the conda environment with r-base.')

    if args.export_rdata:
        if shutil.which('Rscript'):
            for rdata in [ROOT/'legacy_rdata/gprot-dynamics.RData', ROOT/'legacy_rdata/arrestin-dynamics.RData']:
                run(['Rscript', ROOT/'scripts/export_rdata_to_csv.R', rdata, ROOT/'processed_results/dynamics'])
        else:
            print('Rscript not available; cannot export .RData. Use the conda environment with r-base.')

    if args.dynamics_backend:
        if args.dynamics_backend in {'bio3d', 'both'}:
            code = dynamics_backend.run_bio3d()
            if code != 0 and args.dynamics_backend == 'bio3d':
                print('Bio3D backend failed; use --dynamics-backend prody for Python fallback or check R/Bio3D installation.')
        if args.dynamics_backend in {'prody', 'both'}:
            dynamics_backend.run_prody_like(cutoff=args.dynamics_cutoff)
        if args.dynamics_backend == 'legacy':
            legacy_script = ROOT/'legacy_dynamics_figures.py'
            run(['python', legacy_script])

    if args.dynamics_from_scratch:
        dynamics_from_scratch.run_all(cutoff=args.dynamics_cutoff)

    if args.consurf_export_fasta:
        consurf_pipeline.export_fasta(chain=args.consurf_chain)
    if args.consurf_export_pdb:
        consurf_pipeline.export_pdb(chain=args.consurf_chain, all_complexes=args.consurf_all_complexes)
    if args.consurf_submit:
        # Current ConSurf structure mode requires .pdb/.cif/.ent/.mmcif plus job name, chain and e-mail.
        if not args.consurf_email:
            raise SystemExit('--consurf-email is required for ConSurf submission')
        manifest = ROOT/'consurf/input/consurf_pdb_manifest.csv'
        if not manifest.exists():
            consurf_pipeline.export_pdb(chain=args.consurf_chain, all_complexes=args.consurf_all_complexes)
        import json
        logs = consurf_auto.submit_pdb_manifest(
            manifest,
            email=args.consurf_email,
            chain=args.consurf_chain,
            headless=False,
            limit=args.consurf_limit,
            click_submit=args.consurf_click_submit,
        )
        (ROOT/'consurf/submission_log.json').write_text(json.dumps(logs, indent=2))
    if args.consurf_fetch or args.consurf_watch_submission_log:
        consurf_auto.watch_submission_log(
            wait_minutes=args.consurf_wait_minutes,
            poll_seconds=args.consurf_poll_seconds,
            collect=args.consurf_collect,
            plot=args.consurf_plot,
        )
        # collection/plot already handled if requested; avoid duplicate below
        args.consurf_collect = False
        args.consurf_plot = False
    if args.consurf_watch_job:
        consurf_auto.watch_job(
            args.consurf_watch_job,
            job_name=args.consurf_job_name,
            wait_minutes=args.consurf_wait_minutes,
            poll_seconds=args.consurf_poll_seconds,
            collect=args.consurf_collect,
            plot=args.consurf_plot,
        )
        args.consurf_collect = False
        args.consurf_plot = False
    if args.consurf_download_url:
        consurf_auto.download_current_url(args.consurf_download_url, job_name=args.consurf_job_name)
    if args.consurf_collect:
        consurf_results_parser.collect_results(ROOT/'consurf/results_manual')
    if args.consurf_plot:
        consurf_results_parser.plot_consurf()

    if not args.skip_figures:
        plot_paper_figures.all_figures()

    if args.paper_exact_fig2 or args.paper_extra_figures:
        try:
            print('Generating paper-exact Figure 2...')
            print(figure2_paper_exact.run_all())
        except Exception as e:
            print(f'Figure 2 paper-exact generation failed: {e}')

    if args.paper_s2_groups or args.paper_extra_figures:
        try:
            print('Generating Figure S2 chemical-group panels...')
            print(figure_s2_groups.run_all())
        except Exception as e:
            print(f'Figure S2 chemical-group generation failed: {e}')

    if args.paper_figure3 or args.paper_extra_figures:
        try:
            print('Generating Figure 3 PyMOL panels...')
            import generate_figure3_pymol
            print(generate_figure3_pymol.run_all())
        except Exception as e:
            print(f'Figure 3 PyMOL generation failed: {e}')

    print('Done. See processed_results/, summary/, images/, consurf/.')

if __name__ == '__main__':
    main()
