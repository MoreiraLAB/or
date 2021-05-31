#!/usr/bin/env python

"""
Define path names and variables
"""

"""
Tables
"""
TEMPLATES_LOCATION = "templates"
WEINSTEIN_NUMBERING = TEMPLATES_LOCATION + "/weinstein_numbering_RMSD.csv"
GPCRDB_NUMBERING = TEMPLATES_LOCATION + "/weinstein_numbering_R_gproteins_opioids.csv"
GPCRDB_ALIGNMENT = TEMPLATES_LOCATION + "/g_prot_special_opioids.csv"
WEINSTEIN_NUMBERING_ORIGINAL = TEMPLATES_LOCATION + "/weinstein_numbering_opioids.csv"
TEMPLATES_NUMBERING = TEMPLATES_LOCATION + "/template_comparisons_numbering_opioids.csv"
DIMERS = TEMPLATES_LOCATION + "/dimers_opioids.csv"
ALIGN_VECTORS = TEMPLATES_LOCATION + "/alignment_vectors_opioids.csv"
FINAL_TABLE = "final"
CONSURF_RESULTS_FILE = "resources/consurf_results.csv"
POVME_RESULTS_FILE = "binding_pocket.csv"

"""
Pymol Sessions
"""
FULL_SESSION = "templates_models.pse"

"""
R associated files
"""
CIRCLE_GRAPH_R = "GPCR_circlize.R"
DYNAMICS_GRAPH_R = "NMA_script.R"
INTERHELICAL_GRAPH_R = "interhelical.R"
R_PATH = "C:/Users/marti/Anaconda3/pkgs/mro-base-3.4.3-0/Scripts/Rscript.exe"

"""
Executables
"""
CHROMEDRIVER = r'C:\Users\A.J. Preto\Documents\chromedriver_win32\chromedriver.exe'
PHANTOMJS = r'C:\Users\A.J. Preto\Documents\phantomjs-2.1.1-windows\phantomjs-2.1.1-windows\bin\phantomjs.exe'

"""
URLs
"""
COCOMAPS_SUBMIT = 'https://www.molnac.unisa.it/BioTools/cocomaps/'
COCOMAPS_OUTPUT = 'https://www.molnac.unisa.it/BioTools/cocomaps/view.psp'
INTERPROSURF = 'http://curie.utmb.edu/usercomplex.html'
CONSURF_URL = "http://consurf.tau.ac.il/index_proteins.php"

"""
Other variables
"""
TEMPLATES_LIST = ["3sn6_Gs","4zwj_Arr","6pwc_Arr","5dgy_Arr","5uz7_Gs","5vai_Gs","5w0p_Arr","6b3j_Gs","6cmo_Gi","6d9h_Gi2","6ddf_Gi","6g79_Go","6gdg_Gs"]
RECEPTORS_LIST = ["DOR","KOR","MOR","NOP"]
STRUCTURAL_FEATURES = ["POLAR area/energy", "APOLAR area/energy", "TOTAL area/energy", "Number of surface atoms", "Number of buried atoms"]
GPROTEIN_SUBSTRUCTURES = ["HN","hns1","S1","s1h1","H1","h1ha","HA", \
					"hahb","HB","hbhc","HC","hchd","HD","hdhe", \
					"HE","hehf","HF","hfs2","S2","s2s3","S3", \
					"s3h2","H2","h2s4","S4","s4h3","H3","h3s5", \
					"S5","s5hg","HG","hgh4","H4","h4s6","S6", \
					"s6h5","H5"]
INTERPROSURF_START = "INTERSURF"
INTERPROSURF_COMPLEX_VAR = "complex"
INTERPROSURF_RESIDUES_VAR = "residues"
COCOMAPS_MOL_VAR = "asa_mol1"
COCOMAPS_MOL_VAR_2 = "asa_mol2"
CONSURF_CSV_DEFAULT = "msa_aa_variety_percentage.csv"
CONSURF_START = "pos"
CONSURF_GRADE_COL = "ConSurf Grade"
INTRA_CHAIN_CA = "intra_chain"
INTER_CHAIN_CA = "inter_chain"
HBOND_VAR = "hbonds"
HBOND_VAR_GPCR = "toDR_hb"
HBOND_VAR_PROTEIN = "toprotein_hb"
TOTAL_STRING_VAR = "total"
COCOMAPS_START = "COCOMAPS"
NONPOLAR_ALIPHATIC = ["GLY","ALA","VAL","LEU","ILE","MET"]
POLAR_UNCHARGED = ["SER","THR","CYS","PRO","ASN","GLN"]
BASIC_CHARGED = ["LYS","ARG","HIS"]
ACIDIC_CHARGED = ["ASP","GLU"]
NONPOLAR_AROMATIC = ["PHE","TYR","TRP"]
AMINOACID_LIST = NONPOLAR_ALIPHATIC + POLAR_UNCHARGED + BASIC_CHARGED + ACIDIC_CHARGED + NONPOLAR_AROMATIC
GROUPS_DICT = {"Non polar aliphatic": NONPOLAR_ALIPHATIC, "Polar uncharged": POLAR_UNCHARGED,
                 "Basic positively charged": BASIC_CHARGED, "Acid negatively charged": ACIDIC_CHARGED,
                 "Nonpolar aromatic": NONPOLAR_AROMATIC}
POLAR_NAMES = ["Non polar aliphatic", "Polar uncharged", "Basic positively charged", "Acid negatively charged", "Nonpolar aromatic"]
HB_TOTAL_PREFIX = "toDR_hb_total"
SB_PREFIX = "novel_SB"
SUBSTRUCTURES_EVALUATED = ["ICL1", "ICL2", "ICL3", "H8"]
RMSD_SECTION = "section_RMSD"
RMSD_COMPLEX = "complex_RMSD"
RMSD_RECEPTOR = "receptor_RMSD"
POVME_LOC = "POVME"

"""
Folder paths
"""
DEFAULT_FOLDER = "C:/Users/marti/OneDrive/Desktop/silverio"
RESULTS_FOLDER = "results"
PROCESSED_RESULTS_FOLDER = "processed_results"
SUMMARY_FOLDER = "summary"
SUMMARY_TABLE = "summary_table"
RMSD_SECTION_SUMMARY = "RMSD_section_summary"
RMSD_COMPLEX_SUMMARY = "RMSD_complex_summary"
RMSD_RECEPTOR_SUMMARY = "RMSD_receptor_summary"
RMSIP_FILES_LOCATION = "DYNAMIC_ANALYSIS/RMSIP"
RMSIP_IMAGES_LOCATION = "images/RMSIP_heatmap"
RMSIP_BARIMAGES_LOCATION = "images/RMSIP_barplot"
FLUCT_FILES_LOCATION = "DYNAMIC_ANALYSIS/FLUCTUATIONS"
FLUCT_FILES_LOCATION_DIFF = "FLUCT_DIFF"
FLUCT_IMAGES_LOCATION = "images/fluctuations"
FLUCT_IMAGES_LOCATION_DIFF = "images/fluctuations_diff"
CENTRAL_FILES_LOCATION = "NODE_CENTRALITY"
CENTRAL_IMAGES_LOCATION = "images/node_centrality"
DELPHI_EXE = DEFAULT_FOLDER + "/Delphi/executable/delphicpp.exe"
DELPHI_RESULTS = DEFAULT_FOLDER + "/delphi_results"

"""
User data
"""
DEFAULT_EMAIL = "martinsgomes.jose@gmail.com"

"""
Numeric variables
"""
SB_DISTANCE = 4.0
HB_DISTANCE = 8.0
DECIMAL_HOUSES = 2