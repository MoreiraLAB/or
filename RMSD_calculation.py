#!/usr/bin/env python

"""
Calculate total, receptor and section RMSD for a pymol session
"""

import pymol
from pymol import cmd
import os
import csv
import raw_parser
import sys
import gpcr_variables
from gpcr_variables import WEINSTEIN_NUMBERING, TEMPLATES_NUMBERING, TEMPLATES_LIST, \
                            DIMERS, FULL_SESSION, RMSD_SECTION, RMSD_COMPLEX, RMSD_RECEPTOR, \
                            SUMMARY_FOLDER, RMSD_SECTION_SUMMARY, RMSD_COMPLEX_SUMMARY, RMSD_RECEPTOR_SUMMARY


def dimer_gen(input_dict):

    """
    Generate the dimers from the receptor dictionary
    """
    output_dimers = []
    for receptor in input_dict.keys():
        for partner in input_dict[receptor]:
            output_dimers.append("-".join([receptor, partner]))
    return output_dimers

def total_rmsd(input_structure_1, input_structure_2):
    
    """
    Calculate total RMSD for two structures
    """
    output_file =  SUMMARY_FOLDER + "/" + RMSD_COMPLEX + "_" + input_structure_1 + ".csv"
    totalRMSD=cmd.super(input_structure_1, input_structure_2)
    print(input_structure_1, input_structure_2, totalRMSD[0])
    with open(output_file, "a") as csvfile:
        csv_writer = csv.writer(csvfile, delimiter=',')
        to_write_row = [input_structure_1, input_structure_2, str(totalRMSD[0])]
        csv_writer.writerow(to_write_row)

def section_rmsd(input_structure_1, input_structure_1_section, input_structure_2, input_structure_2_section):

    """
    Calculate section RMSD for two substructures
    """
    output_file = SUMMARY_FOLDER + "/" + RMSD_SECTION + "_" + input_structure_1 + ".csv" 
    section_1 = input_structure_1 + " and chain A and resi " + input_structure_1_section
    section_2 = input_structure_2 + " and chain A and resi " + input_structure_2_section
    sectionRMSD=cmd.super(section_1, section_2)
    print(section_1, section_2, sectionRMSD[0])
    with open(output_file, "a") as csvfile:
        csv_writer = csv.writer(csvfile, delimiter=',')
        to_write_row = [input_structure_1, input_structure_1_section, input_structure_2, input_structure_2_section, str(sectionRMSD[0])]
        csv_writer.writerow(to_write_row)

def chain_rmsd(input_structure_1, input_structure_2):

    """
    Calculate chain RMSD for two structures
    """
    output_file = SUMMARY_FOLDER + "/" + RMSD_RECEPTOR  + "_" + input_structure_1 + ".csv"
    chain_1 = input_structure_1 + " and chain A"
    chain_2 = input_structure_2 + " and chain A"
    chainRMSD=cmd.super(chain_1, chain_2)
    print(chain_1, chain_2, chainRMSD[0])
    with open(output_file, "a") as csvfile:
        csv_writer = csv.writer(csvfile, delimiter=',')
        to_write_row = [input_structure_1, input_structure_2, str(chainRMSD[0])]
        csv_writer.writerow(to_write_row)

def locate_partner(input_dict, mode = "arrestins"):

    """
    Generate partners lists
    """
    if mode == "arrestins":
        for receptor in input_dict.keys():
            return input_dict[receptor][0:2]
    if mode == "gproteins":
        output_list = []
        for receptor in input_dict.keys():
            for gprot in  input_dict[receptor][2:]:
                if gprot not in output_list:
                    output_list.append(gprot)
        return output_list 

def perform_calculations(template_sections, template, template_parts, receptors, model_parts, model, total_written):
    
    """
    Deploy total and section RMSD calculations
    """
    for section_template, section_model in zip(template_sections[template_parts[0]],receptors[model_parts[0]]):
        if total_written == False:
            total_written = True
            total_rmsd(template, model)
            chain_rmsd(template, model)
            section_rmsd(template, section_template, model, section_model)
        else:
            section_rmsd(template, section_template, model, section_model)

def calculate_RMSD(receptors, template_sections, templates, input_models, arrestins, G_proteins):

    """
    Deploy the entire pipeline
    """
    for model in input_models:
        model_parts = model.split("-")
        if model_parts[1] != "":
            for template in templates:
                print("Currently evaluating:",model_parts, template)
                template_parts = template.split("_")
                total_written = False
                if template_parts[1] == "Arr" and model_parts[1] in arrestins:
                    perform_calculations(template_sections, template, template_parts, receptors, model_parts, model, total_written)
                if template_parts[1] != "Arr" and model_parts[1] in G_proteins:
                    perform_calculations(template_sections, template, template_parts, receptors, model_parts, model, total_written)

"""
Variable initialization

- Structures for reference and their TM definition
- Models to calculate the RMSD and their TM definitions

"""

DXR = raw_parser.retrieve_clean_horizontal(WEINSTEIN_NUMBERING, delimiter = ";")
template_dict = raw_parser.retrieve_clean_horizontal(TEMPLATES_NUMBERING, delimiter = ";")
raw_models = raw_parser.retrieve_clean(DIMERS)
models = dimer_gen(raw_models)
templates = TEMPLATES_LIST
arrestins = locate_partner(raw_models, mode = "arrestins")
G_proteins = locate_partner(raw_models, mode = "gproteins")

"""
PyMOL session file
"""
session = FULL_SESSION
cmd.load(session)
calculate_RMSD(DXR, template_dict, templates, models, arrestins, G_proteins)

"""
Join the individual tables
"""
raw_parser.join_tables_RMSD(output_name = RMSD_COMPLEX_SUMMARY, table_type = RMSD_COMPLEX.split("_")[0], col_position = 2, mode = "complex")
raw_parser.join_tables_RMSD(output_name = RMSD_SECTION_SUMMARY, table_type = RMSD_SECTION.split("_")[0], col_position = 4, mode = "section")
raw_parser.join_tables_RMSD(output_name = RMSD_RECEPTOR_SUMMARY, table_type = RMSD_RECEPTOR.split("_")[0], col_position = 2, mode = "complex")