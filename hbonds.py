#!/usr/bin/env python

"""
Analyse interface output from COCOMAPS to yield
residue type and group information
"""

import os
import csv
import raw_parser
import gpcr_variables
from gpcr_variables import DEFAULT_FOLDER, RESULTS_FOLDER, PROCESSED_RESULTS_FOLDER, \
                            ALIGN_VECTORS, DIMERS, COCOMAPS_START, HBOND_VAR, HBOND_VAR_GPCR, \
                            HBOND_VAR_PROTEIN, TOTAL_STRING_VAR

__author__ = "A.J. Preto"
__email__ = "martinsgomes.jose@gmail.com"
__group__ = "Data-Driven Molecular Design"
__group_leader__ = "Irina S. Moreira"
__project__ = "GPCRs"

def process_input(input_file):

    """
    Process a simple input file into list of lists.
    Note that this is a generator
    """
    opened_file = open(input_file, "r").readlines()
    for row in opened_file:
        row = row.replace("\n","")
        new_row = row.split()
        yield new_row

def analyse_input(input_table, target_sequence, current_chain = "A"):

    """
    Yield the output sequence containing only the interfacial residues
    """
    output_list = []
    corrected_number = 0
    for res_number in range(1,len(target_sequence) + 1):
        if target_sequence[res_number - 1] == '-':
            corrected_number = corrected_number + 1
        new_res_number = res_number - corrected_number
        count = 0
        output_list.append('-')
        for row in input_table:
            if (int(row[1]) == new_res_number) and (row[0] == current_chain):
                if output_list[-1] != '-':
                    output_list[-1] = output_list[-1] + ',' + row[8]
                else:
                    output_list[-1] = row[8]
            if (int(row[5]) == new_res_number) and (row[4] == current_chain):
                if output_list[-1] != '-':
                    output_list[-1] = output_list[-1] + ',' + row[8]
                else:
                    output_list[-1] = row[8]
    return output_list

def write_output_file(output_name, output_data):

    """
    Write the final output file
    """
    new_file = open(output_name,'w')
    csv_writer=csv.writer(new_file,delimiter = ';')
    csv_writer.writerow(output_data)
    new_file.close()

def calculate_total_hbonds(input_list):

    """
    Return the ordered list of Hydrogen bonds per residue
    """
    output_list = []
    for value in input_list:
        if value == "-":
            output_list.append("-")
        else:
            output_list.append(len(value.split(",")))
    return output_list

def build_csv(input_file,seq_DR,seq_partner,out_name_DR,out_name_protein):
    
    """
    Build the csv file containing the Hydrogen bonds analysis
    """
    processed_hb = list(process_input(input_file))
    output_list_DR = analyse_input(processed_hb, seq_DR)
    output_list_protein = analyse_input(processed_hb, seq_partner, current_chain = "B")
    write_output_file(str(out_name_DR + ".csv"), output_list_DR)
    write_output_file(str(out_name_protein + ".csv"), output_list_protein)
    write_output_file(str(out_name_DR + "_" + TOTAL_STRING_VAR + ".csv"), calculate_total_hbonds(output_list_DR))
    write_output_file(str(out_name_protein + "_" + TOTAL_STRING_VAR + ".csv"), calculate_total_hbonds(output_list_protein))

def apply_on_folder(input_alignments, input_dimers):

    """
    Deploy the pipeline for Hydrogen bond analysis.
    Keep in mind the folder structure as it is determinant for the correct
    usage of the pipeline
    """
    for receptor in input_dimers.keys():
        for partner in input_dimers[receptor]:
            hbonds_name_DR = results_folder  + "/" + COCOMAPS_START + "_" + receptor + "-" + partner + "_" + HBOND_VAR + ".txt"
            output_name_DR = processed_results_folder  + "/" + receptor + "_" + partner + "_" + HBOND_VAR_GPCR
            output_name_partner = processed_results_folder  + "/" + receptor + "_" + partner + "_" + HBOND_VAR_PROTEIN
            try:
                build_csv(hbonds_name_DR, input_alignments[receptor], input_alignments[receptor], output_name_DR, output_name_partner)
            except:
                print("Failed to perform action on:", receptor, partner)

"""
Variable initialization
"""
home = DEFAULT_FOLDER
results_folder = home + "/" + RESULTS_FOLDER
processed_results_folder = home + "/" + PROCESSED_RESULTS_FOLDER
aligned_proteins = raw_parser.retrieve_clean(ALIGN_VECTORS)
input_dimers = raw_parser.retrieve_clean(DIMERS)
apply_on_folder(aligned_proteins, input_dimers)

