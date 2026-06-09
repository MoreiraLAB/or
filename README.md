# Opioid Receptor Partners Pipeline — Functional 2026 Reconstruction

This repository is a functional 2026 reconstruction of the original `MoreiraLAB/or` opioid receptor partner-specificity analysis code associated with:

- Carlos A. V. Barreto*, Salete J. Baptista*, A. J. Preto, Daniel Silvério, Rita Melo, Irina S. Moreira. **Decoding partner specificity of opioid receptor family**. *Frontiers in Molecular Biosciences*, 2021. *Joint first authors.*
- A. J. Preto, Carlos A. V. Barreto, Salete J. Baptista, José Guilherme de Almeida, Agostinho Lemos, André Melo, M. Natália D. S. Cordeiro, Zeynep Kurkcuoglu, Rita Melo and Irina S. Moreira. **Understanding the Binding Specificity of G-Protein Coupled Receptors toward G-Proteins and Arrestins: Application to the Dopamine Receptor Family**. *Journal of Chemical Information and Modeling*, 2020.

The aim is to preserve the scientific logic of the original workflow while making the code runnable with current software and current web-server behaviour. The obsolete CoCoMaps and InterProSurf scraping steps have been replaced by local Python calculations. ConSurf remains external and is now handled through a receptor-level structure-upload workflow, because conservation is a property of the receptor sequence and does not need to be recalculated for every receptor-partner complex.

---

## Quick start

```bash
conda env create -f environment.yml
conda activate opioid_receptors_rebuilt
python doctor.py
python call.py
python call.py --export-rdata --skip-structural --skip-interhelical --skip-figures
python plot_paper_figures.py
```

Successful figure generation should report:

```text
{'s1': True, 's2': True, 'fig2': True, 's3': True}
```

Main generated outputs are written to:

```text
processed_results/
summary/
images/
consurf/
```

These folders are runtime outputs and should normally be ignored by Git.

---

## Main reconstruction status

| Component | Status | Notes |
|---|---:|---|
| Structural/interface tables | Working | Local Python replacement for obsolete server scraping. |
| S1-style interhelical plots | Working | Python fallback; optional R/Bio3D reference backend. |
| S2 amino-acid interface percentages | Working | Individual residue heatmaps. |
| S2 chemical-group interface percentages | Added | Groups: nonpolar aliphatic, aromatic, basic positive, acidic negative, polar uncharged. |
| S3 dynamics plots | Working | Legacy `.RData` validation and optional recalculation backends. |
| Figure 2-style interaction signature | Added | Paper-oriented BW residue versus partner-domain/common-interaction matrix. |
| Figure 3-style structural panels | Added | Requires PyMOL. Generates receptor/partner/interface panels. |
| ConSurf | Working with limitations | Default is 4 receptor-level jobs. Avoid all-complex automatic submission. |

---

## What changed relative to the original repository

### Replaced

- Old CoCoMaps scraping.
- Old InterProSurf scraping.
- Manual `chromedriver.exe` management.
- Fragile Selenium selectors based on obsolete server forms.
- R-only plotting for several downstream summaries.

These are replaced by local Python calculations for:

- Interface/contact tables;
- Receptor-side and partner-side amino-acid composition;
- Receptor-side and partner-side chemical-group composition;
- Hydrogen-bond and salt-bridge style summaries;
- Receptor/partner domain mapping;
- Figure 2-style common interaction signatures;
- ConSurf result parsing and integration.

### Retained and modernised

- Ballesteros-Weinstein numbering from the original CSV metadata.
- Receptor domains: TM1–TM7, ICLs and H8.
- Partner domains: H5, h4s6, finger loop, C-loop, lariat loop and related regions.
- Bio3D/R backend for dynamics/interhelical calculations where scientific equivalence with the original analysis is required.
- Python fallback backend for reproducibility on machines without R/Bio3D.
- Legacy `.RData` outputs as validation data for Figure S3-style dynamics plots.
- ConSurf as the external conservation server.

---

## Repository layout

```text
.
├── call.py                         # Main orchestration entry point
├── doctor.py                       # Environment/repository diagnostic
├── pdb_analysis.py                 # Local structural/interface analysis
├── plot_paper_figures.py           # S1/S2/Fig2/S3-style plotting
├── figure2_paper_exact.py          # Paper-oriented Figure 2 reconstruction
├── figure_s2_groups.py             # Chemical-group S2 plots
├── generate_figure3_pymol.py       # PyMOL-based Figure 3 panels
├── interhelical_bio3d.R            # R/Bio3D interhelical backend
├── interhelical_python_fallback.py # Python interhelical fallback
├── dynamics_bio3d.R                # R/Bio3D dynamics backend
├── dynamics_from_scratch.py        # Python Cα elastic-network fallback
├── dynamics_backend.py             # Dynamics dispatcher
├── compare_dynamics_backends.py    # Bio3D vs Python comparison
├── consurf_pipeline.py             # ConSurf PDB/FASTA export + parser wrapper
├── consurf_auto.py                 # ConSurf submit/watch/download automation
├── consurf_selenium_download_job.py# Browser-based downloader for finished jobs
├── consurf_results_parser.py       # Parser for ConSurf Download-all archives
├── metadata/                       # Original numbering/domain CSVs
├── structural_complexes/           # PDB complexes for interface/ConSurf input
├── dynamics_complexes/             # PDBs for dynamics calculations
├── legacy_rdata/                   # Original processed dynamic RData objects
├── processed_results/              # Generated tables; ignored by Git
├── summary/                        # Generated summaries; ignored by Git
├── images/                         # Generated figures; ignored by Git
└── consurf/                        # ConSurf runtime folders; ignored by Git
```

---

## Metadata files

The metadata CSVs are essential. They define the structural ontology used throughout the analysis:

```text
metadata/weinstein_numbering_opioids.csv
metadata/weinstein_R_numbering_opioids.csv
metadata/alignment_vectors_opioids.csv
metadata/template_comparisons_numbering_opioids.csv
metadata/g_proteins_final_opioids.csv
metadata/arrestins_final_opioids.csv
metadata/g_prot_specials.csv
metadata/arrestins_special.csv
metadata/partners_alternative_numbering.csv
metadata/receptors_special.csv
metadata/dimers_opioids.csv
```

They are used for:

- Ballesteros-Weinstein labels;
- TM/ICL/H8 definitions;
- G-protein and arrestin domain definitions;
- interhelical vectors/distances;
- receptor/partner motif mapping;
- paper-style interaction signatures.

---

## Environment

Recommended environment:

```yaml
name: opioid_receptors_rebuilt

channels:
  - conda-forge
  - bioconda
  - defaults

dependencies:
  - python=3.10
  - pip
  - pandas
  - numpy
  - scipy
  - scikit-learn
  - matplotlib
  - seaborn
  - pyyaml
  - openpyxl
  - biopython
  - prody
  - requests
  - beautifulsoup4
  - lxml
  - selenium>=4.25
  - r-base>=4.3
  - r-bio3d
  - r-tidyverse
  - r-ggplot2
  - r-ggrepel
  - r-cowplot
  - r-svglite
  - r-optparse
  - pip:
      - webdriver-manager>=4.0.2
```

Optional for Figure 3:

```bash
conda install -c conda-forge pymol-open-source pillow
```

---

## Structural/interface pipeline

Run the full local pipeline, excluding ConSurf:

```bash
python call.py --skip-figures
python call.py --export-rdata --skip-structural --skip-interhelical --skip-figures
python call.py --dynamics-backend legacy --skip-structural --skip-interhelical
python plot_paper_figures.py
```

This generates:

- S1-style interhelical distance plots;
- S2-style amino-acid interface percentage plots;
- Figure 2-style interaction/contact summaries where possible;
- S3-style dynamics plots;
- structural/contact/interface CSVs.

---

## Paper-oriented extra figures

These scripts improve agreement with the original article and supplementary figures.

### Figure 2-style interaction signature

```bash
python call.py --paper-exact-fig2 --skip-structural --skip-interhelical --skip-figures
```

Outputs:

```text
images/figure2_paper_exact.png
images/figure2_paper_exact.svg
summary/figure2_paper_exact_matrix.csv
```

This aims to reproduce the paper logic: receptor BW positions versus partner domains/residue motifs, focusing on common interaction signatures across receptor groups.

### S2 chemical-group plots

```bash
python call.py --paper-s2-groups --skip-structural --skip-interhelical --skip-figures
```

Outputs:

```text
images/figure_s2_groups_receptor.png
images/figure_s2_groups_partner.png
images/figure_s2_groups_receptor.svg
images/figure_s2_groups_partner.svg
summary/interface_group_percentages_receptor.csv
summary/interface_group_percentages_partner.csv
```

Chemical groups used:

```text
nonpolar_aliphatic: A, V, L, I, M
aromatic: F, Y, W
basic_positive: K, R, H
acidic_negative: D, E
polar_uncharged: S, T, N, Q, C, G, P
```

### Figure 3-style PyMOL structural panels

Requires PyMOL:

```bash
conda install -c conda-forge pymol-open-source pillow
python call.py --paper-figure3 --skip-structural --skip-interhelical --skip-figures
```

Outputs:

```text
images/figure3_panel_gio.png
images/figure3_panel_gs.png
images/figure3_panel_gq.png
images/figure3_panel_g12.png
images/figure3_paper.png
```

Run all paper-oriented extra figures:

```bash
python call.py --paper-extra-figures --skip-structural --skip-interhelical --skip-figures
```

---

## Dynamics and legacy RData

The original `.RData` files are retained as validation/reference outputs:

```text
legacy_rdata/gprot-dynamics.RData
legacy_rdata/arrestin-dynamics.RData
```

They contain final processed results, not raw trajectories:

```text
gprot_bc_scores             52 x 8
gprot_fluc_fold_change      52 x 11
arrestin_bc_scores          16 x 8
arrestin_fluc_fold_change   16 x 11
```

Domain order:

```text
TM1, TM2, TM3, TM4, TM5, TM6, TM7, H8
```

Commands:

```bash
python call.py --inspect-rdata --skip-structural --skip-interhelical --skip-figures
python call.py --export-rdata --skip-structural --skip-interhelical --skip-figures
python call.py --dynamics-backend legacy --skip-structural --skip-interhelical
python call.py --dynamics-backend bio3d --skip-structural --skip-interhelical --skip-figures
python call.py --dynamics-backend prody --skip-structural --skip-interhelical --skip-figures
python call.py --dynamics-backend both --skip-structural --skip-interhelical --skip-figures
python compare_dynamics_backends.py
```

---

# ConSurf workflow

## Important scientific point

ConSurf conservation scores are properties of the receptor sequence. They do not need to be recalculated for every receptor-partner complex. Therefore, for the opioid receptor dataset, the scientifically necessary ConSurf calculations are normally one representative structure per receptor:

```text
DOR
KOR
MOR
NOP
```

The pipeline therefore exports **4 representative ConSurf jobs by default** and maps receptor-level conservation scores back onto all receptor-partner interfaces by receptor and residue position.

Only use `--consurf-all-complexes` if you explicitly need a separate ConSurf calculation for every receptor-partner PDB complex. That is usually unnecessary and can trigger TAU/HPC anti-automation policies.

## 1. Export the recommended ConSurf PDB manifest

Default: one representative PDB per receptor.

```bash
python call.py --consurf-export-pdb --consurf-chain A --skip-structural --skip-interhelical --skip-figures
```

This creates:

```text
consurf/input/consurf_pdb_manifest.csv
```

Check how many jobs will be submitted:

```bash
wc -l consurf/input/consurf_pdb_manifest.csv
cat consurf/input/consurf_pdb_manifest.csv
```

Expected: header + up to 4 receptor-representative jobs.

If you really need all complexes:

```bash
python call.py --consurf-export-pdb --consurf-chain A --consurf-all-complexes --skip-structural --skip-interhelical --skip-figures
```

## 2. Browser-review submission, no click

Use this first. It opens Chrome, uploads the PDB, fills job name, email and chain, but does not click Submit.

```bash
python call.py \
  --consurf-submit \
  --consurf-email your@email.com \
  --consurf-chain A \
  --consurf-limit 1 \
  --skip-structural --skip-interhelical --skip-figures
```

The log is written to:

```text
consurf/submission_log.json
```

Each job should show:

```json
"uploaded": "True",
"job_name_filled": "True",
"email_filled": "True",
"chain_filled": "True"
```

## 3. Submit ConSurf jobs

Submit the 4 receptor-representative jobs:

```bash
python call.py \
  --consurf-submit \
  --consurf-email your@email.com \
  --consurf-chain A \
  --consurf-click-submit \
  --skip-structural --skip-interhelical --skip-figures
```

If you want to be conservative, submit one at a time:

```bash
python call.py \
  --consurf-submit \
  --consurf-email your@email.com \
  --consurf-chain A \
  --consurf-limit 1 \
  --consurf-click-submit \
  --skip-structural --skip-interhelical --skip-figures
```

Avoid submitting all receptor-partner complexes automatically unless needed.

## 4. Watch job numbers and download results

After submission, `consurf/submission_log.json` stores `run_number`/`progress_url` when captured.

```bash
python call.py \
  --consurf-watch-submission-log \
  --consurf-wait-minutes 180 \
  --consurf-poll-seconds 60 \
  --consurf-collect \
  --consurf-plot \
  --skip-structural --skip-interhelical --skip-figures
```

The ConSurf server can show an intermediate state:

```text
The files are being zipped. Please wait.
```

If this happens, wait a few minutes and run the watcher/downloader again.

## 5. Browser-based download for a finished job

If the requests-based watcher does not see the result links, use the Selenium downloader. This is the most reliable method after the job has reached:

```text
https://consurf.tau.ac.il/final_output/?number=<RUN_NUMBER>
```

Example:

```bash
python consurf_selenium_download_job.py \
  --number 1780989561 \
  --job-name KOR-Arr3_6PWC \
  --wait 90
```

This saves the downloaded `*_ConSurf.tar.gz` under:

```text
consurf/results_manual/<job-name>/
```

## 6. Parse ConSurf results and integrate them

If the `.tar.gz` archive exists:

```bash
python consurf_results_parser.py --extract consurf/results_manual/KOR-Arr3_6PWC/1780989561_ConSurf.tar.gz
python consurf_results_parser.py --collect
python consurf_results_parser.py --plot
```

Outputs:

```text
processed_results/consurf_scores.csv
processed_results/interface_contacts_with_consurf.csv
summary/consurf_domain_summary.csv
images/consurf_profile_<job>.png
images/consurf_profile_<job>.svg
images/consurf_domain_summary.png
images/consurf_domain_summary.svg
```

---

## Final validation sequence excluding ConSurf

```bash
python doctor.py
python -m py_compile *.py
python call.py --skip-figures
python call.py --export-rdata --skip-structural --skip-interhelical --skip-figures
python call.py --dynamics-backend legacy --skip-structural --skip-interhelical
python call.py --paper-exact-fig2 --paper-s2-groups --skip-structural --skip-interhelical --skip-figures
python plot_paper_figures.py
```

Optional Figure 3:

```bash
conda install -c conda-forge pymol-open-source pillow
python call.py --paper-figure3 --skip-structural --skip-interhelical --skip-figures
```


## Troubleshooting

### `zsh: no matches found: https://...progress/?number=...`

Quote URLs containing `?`:

```bash
python consurf_auto.py --download-current-url "https://consurf.tau.ac.il/progress/?number=1780947241" --job-name DOR
```

### `Your request has violated HPC secure policy`

Do not submit all receptor-partner complexes automatically. Use the default 4 receptor-representative jobs only, or submit one job manually at a time.

### `The files are being zipped. Please wait.`

The job finished but the archive is not ready yet. Wait and rerun the watcher or the Selenium downloader.

### `stand_alone_consurf_1.05.tar.gz`

This is the standalone software archive from the ConSurf menu, not a job result. Delete it with:

```bash
find consurf/results_manual -name "stand_alone_consurf_1.05.tar.gz" -delete
```
