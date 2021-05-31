#!/usr/bin/env python

"""
Process salt bridges output
"""

import os
import csv
import raw_parser
import gpcr_variables
from gpcr_variables import DEFAULT_FOLDER, RESULTS_FOLDER, PROCESSED_RESULTS_FOLDER, ALIGN_VECTORS, DIMERS

__author__ = "A.J. Preto"
__email__ = "martinsgomes.jose@gmail.com"
__group__ = "Data-Driven Molecular Design"
__group_leader__ = "Irina S. Moreira"
__project__ = "GPCRs"

def process_table(input_table, delimiter = ";"):

    """
    Process the input table, skipping the header
    """
    opened_table = open(input_table, "r").readlines()
    for row in opened_table:
        row = row.replace("\n","")
        new_row = row.split(delimiter)
        if len(new_row) > 1:
            yield new_row

def identify_unique(input_row, threshold = 4.0):

    """
    Identify the output values and process them
    """
    output_values = []
    for value in input_row:
        if float(value) != 0.0 and float(value) < threshold:
            output_values.append(float(value))
    if len(output_values) == 0:
        return None
    if len(output_values) == 1:
        return str(output_values[0])
    if len(output_values) > 1:
        return ",".join(output_values)

def column_into_row(input_table, input_col_number):

    """
    Convert a table column into row for further iterations
    """
    output_col = []
    for row in input_table:
        output_col.append(row[input_col_number])
    return output_col

def build_csv(input_file, target_alignment, output_name, receptor = True):

    """
    Construct the final output final with the aligned salt bridges
    """
    output_list = []
    corrected_number = 1
    if receptor == True:
        for res in target_alignment:
            try:
                if res != "-":
                    corrected_number += 1 
                    value = identify_unique(input_file[corrected_number][1:])
                    if value == None:
                        value = "-"
                    output_list.append(value)
                else:
                    output_list.append("-")
            except:
                output_list.append("-")

    if receptor == False:
        for res in target_alignment:
            try:
                if res != "-":
                    corrected_number += 1
                    target_col = column_into_row(input_file, corrected_number)
                    value = identify_unique(target_col[1:])
                    if value == None:
                        value = "-"
                    output_list.append(value)
                else:
                    output_list.append("-")
            except:
                output_list.append("-")
    return output_list

def apply_on_folder(input_alignments, input_dimers):

    """
    Deploy the processing pipeline, pay attention to the folder structure
    """
    for receptor in input_dimers.keys():
        for partner in input_dimers[receptor]:
            input_SB = processed_results_folder  + "/novel_sb_" + receptor + "-" + partner + ".csv"
            output_name_DR = processed_results_folder  + "/" + receptor + "_" + partner + "_toDR_sb.csv"
            output_name_partner = processed_results_folder  + "/" + receptor + "_" + partner + "_toprotein_sb.csv"
            try:
                processed_table = list(process_table(input_SB))
                build_csv(processed_table, input_alignments[receptor], output_name_DR, receptor = True)
                build_csv(processed_table, input_alignments[partner], output_name_partner, receptor = False)
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