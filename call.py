#!/usr/bin/env python

"""
Apply analysis for all the .pdb files in the folder
"""

import os
import csv
import raw_parser
import gpcr_variables
from gpcr_variables import CIRCLE_GRAPH_R, R_PATH, DYNAMICS_SCRIPT_R, FINAL_TABLE, \
                            INTERHELICAL_GRAPH_R, RESULTS_FOLDER, PROCESSED_RESULTS_FOLDER, \
                            SUMMARY_FOLDER, IMAGES_FOLDER

__author__ = "A.J. Preto"
__email__ = "martinsgomes.jose@gmail.com"
__group__ = "Data-Driven Molecular Design"
__group_leader__ = "Irina S. Moreira"
__project__ = "GPCRs"

def call_all(circle_script, dynamics_script, interhelical_script, R_path = "Rscript", download_interface_data = False):

    """
    Import the Python modules to make them run.
    Run the R script separately 
    "submit_cocomaps_interprosurf" and the R script are the more 
    time consuming steps
    """
    create_folders = [RESULTS_FOLDER, PROCESSED_RESULTS_FOLDER, SUMMARY_FOLDER, IMAGES_FOLDER]
    to_create_folders = []
    for check_folder in create_folders:
        if not os.path.isdir(check_folder):
            to_create_folders.append(check_folder)
    for current_folder in to_create_folders:
        os.mkdir(current_folder)
    if download_interface_data == True:
        import submit_cocomaps_interprosurf
    import cocomaps_aligned
    import intersurf_aligned
    import interface
    import hbonds
    import salt_bridges
    import SB_aligned
    import interactions
    import generate_summary
    import interhelical
    import adapt_partner

    interhelical_graph_command = R_path + " " + interhelical_script
    os.system(interhelical_graph_command)

    circle_graph_command = R_path + " " + circle_script
    os.system(circle_graph_command)

    dynamics_graph_command = R_path + " " + dynamics_script
    os.system(dynamics_graph_command)

Rscript_circle = CIRCLE_GRAPH_R
R_executable = R_PATH
Rscript_interhelical = INTERHELICAL_GRAPH_R
Rscript_dynamics = DYNAMICS_SCRIPT_R
call_all(Rscript_circle, Rscript_dynamics, Rscript_interhelical, \
                R_path = R_executable, download_interface_data = True)

