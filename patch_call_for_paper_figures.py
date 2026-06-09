###################################################################################################
############################ MOREIRA LAB - DATA DRIVEN MOLECULAR DESIGN ###########################
###################################################################################################

#!/usr/bin/env python3
from pathlib import Path
p = Path('call.py')
txt = p.read_text()

# Add imports if missing
for imp in ['import figure2_paper_exact', 'import figure_s2_groups']:
    if imp not in txt:
        txt = txt.replace('import plot_paper_figures\n', 'import plot_paper_figures\n' + imp + '\n')
# Figure3 PyMOL optional: import inside execution to avoid failing if PyMOL deps absent.

# Add argparse flags after skip-figures
needle = "p.add_argument('--skip-figures', action='store_true')"
insert = """p.add_argument('--skip-figures', action='store_true')
    p.add_argument('--paper-exact-fig2', action='store_true', help='Generate Figure 2-style common OR-G-protein interaction signature.')
    p.add_argument('--paper-s2-groups', action='store_true', help='Generate Figure S2 chemical-group interface percentage heatmaps.')
    p.add_argument('--paper-figure3', action='store_true', help='Generate Figure 3-style PyMOL structural panels if PyMOL is installed.')
    p.add_argument('--paper-extra-figures', action='store_true', help='Generate all extra paper-matching figures: exact Fig2, S2 groups and Figure3.')"""
if needle in txt and '--paper-exact-fig2' not in txt:
    txt = txt.replace(needle, insert)

# Add execution block before final Done print
marker = "    print('Done. See processed_results/, summary/, images/, consurf/.')"
block = """    if args.paper_exact_fig2 or args.paper_extra_figures:
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

"""
if marker in txt and 'Generating paper-exact Figure 2' not in txt:
    txt = txt.replace(marker, block + marker)

p.write_text(txt)
print('call.py patched with --paper-exact-fig2, --paper-s2-groups, --paper-figure3, --paper-extra-figures')
