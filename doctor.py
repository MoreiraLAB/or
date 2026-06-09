###################################################################################################
############################ MOREIRA LAB - DATA DRIVEN MOLECULAR DESIGN ###########################
###################################################################################################

from pathlib import Path
import importlib.util
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
required = [
    'metadata/weinstein_numbering_opioids.csv',
    'metadata/weinstein_R_numbering_opioids.csv',
    'metadata/g_proteins_final_opioids.csv',
    'metadata/arrestins_final_opioids.csv',
    'metadata/alignment_vectors_opioids.csv',
    'metadata/partners_alternative_numbering.csv',
    'structural_complexes',
    'dynamics_complexes',
    'legacy_rdata/gprot-dynamics.RData',
    'legacy_rdata/arrestin-dynamics.RData',
    'dynamics_bio3d.R',
    'interhelical_bio3d.R',
]
issues = 0
for r in required:
    if not (ROOT / r).exists():
        print('MISSING', r)
        issues += 1

for mod in ['pandas', 'numpy', 'matplotlib', 'Bio', 'sklearn']:
    if importlib.util.find_spec(mod) is None:
        print('MISSING PYTHON MODULE', mod)
        issues += 1

for py in ROOT.glob('*.py'):
    rc = subprocess.call([sys.executable, '-m', 'py_compile', str(py)], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    if rc:
        print('SYNTAX ERROR', py.name)
        issues += 1

rscript = subprocess.call(['which', 'Rscript'], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL) == 0
print('Rscript:', 'available' if rscript else 'not found (needed for Bio3D/RData backend)')
if rscript:
    cmd = ['Rscript', '-e', "if(!requireNamespace('bio3d', quietly=TRUE)) quit(status=2)"]
    bio3d = subprocess.call(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL) == 0
    print('R package bio3d:', 'available' if bio3d else 'not found')
    if not bio3d:
        print('Install with conda env create -f environment.yml, or in R: install.packages("bio3d")')
else:
    print('Bio3D backend will be unavailable; Python/ProDy-like fallback still works.')

print('OK' if issues == 0 else f'{issues} issue(s)')
sys.exit(1 if issues else 0)
