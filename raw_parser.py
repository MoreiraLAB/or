#!/usr/bin/env python

"""
File parser to be used by other functions
"""

import os
import csv
import gpcr_variables
from gpcr_variables import SUMMARY_FOLDER, SUMMARY_TABLE, FINAL_TABLE, \
                            DEFAULT_FOLDER
import pandas as pd

__author__ = "A.J. Preto"
__email__ = "martinsgomes.jose@gmail.com"
__group__ = "Data-Driven Molecular Design"
__group_leader__ = "Irina S. Moreira"
__project__ = "GPCRs"

def retrieve_clean(input_csv, delimiter = ","):

    """
    Clean a .csv file in order to be used by other functions
    """
    alignments_dict = {}
    opened_file = open(input_csv, "r").readlines()
    for row in opened_file:
        row = row.replace("\n","")
        row = row.replace("\r","")
        new_row = row.split(delimiter)
        alignments_dict[new_row[0]] = [i for i in new_row[1:]] 
    return alignments_dict

def retrieve_clean_weinstein(input_csv, delimiter = ",", with_x = True):

    """
    Clean a .csv file in order to be used by other functions
    """
    weinstein_dict = {}
    opened_file = open(input_csv, "r").readlines()
    for row in opened_file:
        row = row.replace("\n","")
        row = row.replace("\r","")
        new_row = row.split(delimiter)
        if with_x == False:
            weinstein_dict["_".join([new_row[0],new_row[1]])] = [i for i in new_row[2:-1]]
        else:
            weinstein_dict["_".join([new_row[0],new_row[1]])] = [i for i in new_row[2:]]
    return weinstein_dict

def retrieve_clean_horizontal(input_csv, delimiter = ","):

    """
    Clean a .csv file in order to be used by other functions, horizontal wise only
    """
    alignments_dict = {}
    opened_file = open(input_csv, "r").readlines()
    for row in opened_file:
        row = row.replace("\n","")
        row = row.replace("\r","")
        new_row = row.split(delimiter)
        alignments_dict[new_row[0]] = new_row[1:]
    return alignments_dict

def intersurf_parser(input_file, delimiter = ","):

    """
    Prepare intersurf values for final summary table
    """
    import pandas as pd
    opened_file = pd.read_csv(input_file, sep = delimiter)
    transposed_table = opened_file.transpose()
    return transposed_table.values.tolist()

def hb_total_counter(input_file, delimiter = ",", default_empty = "-"):

    """
    Open Hydrogen-Bonds associated file and count the bonds
    """
    opened_file = open(input_file,"r").readlines()[0].replace("\n","").split(delimiter)
    total = 0
    for value in opened_file:
        if value != default_empty:
            total += int(value)
    return total

def sb_total_counter(input_file, delimiter = ","):

    """
    Open Salt-Bridges associated file and count the bonds, total
    """
    import pandas as pd
    from gpcr_variables import SB_DISTANCE

    opened_file = open(input_file, "r").readlines()
    total_sb = 0
    for row in opened_file[1:]:
        row = row.split(delimiter)
        for cell in row[1:]:
            if len(cell.split(",")) > 1:
                for bond in cell.split(","):
                    if float(bond) < SB_DISTANCE:
                        total_sb += 1
            elif float(cell) != 0 and (float(cell) < SB_DISTANCE):
                total_sb += 1
    return total_sb

def cocomaps_parser(input_file):

    """
    Process cocomaps raw receptor ASA file
    """
    import pandas as pd

    opened_file = open(input_file, "r").readlines()
    processed_file = []
    for row in opened_file:
        row = row.replace("\n","").split()
        output_row = [row[1], row[-1]]
        processed_file.append(output_row)
    return processed_file

def parse_hb_compare(input_hb_file, delimiter = ","):

    """
    Get the original sequence to further calculate hydrogen bonds on substructures
    """
    from gpcr_variables import ALIGN_VECTORS
    opened_file = open(input_hb_file,"r").readlines()[0].replace("\n","").split(delimiter)
    current_receptor = input_hb_file.split("/")[1].split("_")[0]
    alignments = retrieve_clean(ALIGN_VECTORS)
    receptor_alignment = alignments[current_receptor]
    original_sequence_hb = []
    for hb_value, receptor_aa in zip(opened_file, receptor_alignment):
        if receptor_aa != "-":
            original_sequence_hb.append(hb_value)
    return original_sequence_hb

def join_tables(folder = SUMMARY_FOLDER, output_name = SUMMARY_TABLE, delimiter = ",", table_type = FINAL_TABLE, col_position = 1):

    """
    Generate a single summary table from the individual ones
    """
    os.chdir(SUMMARY_FOLDER)
    output_file = SUMMARY_FOLDER + "/" + output_name + ".csv"
    first = True
    for files in os.listdir(os.getcwd()):
        if files.endswith(".csv") and first == False and files.split("_")[0] == table_type:
            current_table = pd.read_csv(files)
            processed_table = pd.concat([processed_table, current_table.iloc[:,col_position]],axis = 1)
        if files.endswith(".csv") and first == True and files.split("_")[0] == table_type:
            current_table = pd.read_csv(files)
            processed_table = pd.DataFrame(current_table)
            first = False
    os.chdir(DEFAULT_FOLDER)
    processed_table.to_csv(output_file, sep = delimiter)

def join_tables_RMSD(folder = SUMMARY_FOLDER, output_name = SUMMARY_TABLE, delimiter = ",", table_type = FINAL_TABLE, col_position = 1, mode = "other"):

    """
    Generate a single summary table from the individual ones. Specific for RMSD associated tables
    """
    os.chdir(SUMMARY_FOLDER)
    output_file = SUMMARY_FOLDER + "/" + output_name + ".csv"
    first = True
    for files in os.listdir(os.getcwd()):
        if files.endswith(".csv") and first == False and files.split("_")[0] == table_type:
            current_table = pd.read_csv(files, header = None)
            if mode == "section":
                proper_rows = []
                holder_complex = ""
                holder_row = []
                for row in current_table.values.tolist():
                    if row[2] != holder_complex:
                        if holder_complex != "":
                            proper_rows.append(holder_row)
                        holder_row = row
                        holder_complex = row[2]
                    else:
                        holder_row.append(row[-1])
                proper_rows.append(holder_row)
                processed_table = pd.concat([processed_table, pd.DataFrame(proper_rows)],axis = 0)
            if mode == "complex":
                processed_table = pd.concat([processed_table, current_table],axis = 0)
        if files.endswith(".csv") and first == True and files.split("_")[0] == table_type:
            current_table = pd.read_csv(files, header = None)
            if mode == "section":
                proper_rows = []
                holder_complex = ""
                holder_row = []
                for row in current_table.values.tolist():
                    if row[2] != holder_complex:
                        if holder_complex != "":
                            proper_rows.append(holder_row)
                        holder_row = row
                        holder_complex = row[2]
                    else:
                        holder_row.append(row[-1])
                proper_rows.append(holder_row)
                processed_table = pd.DataFrame(proper_rows)
            if mode == "complex":
                processed_table = pd.DataFrame(current_table)
            first = False
    os.chdir(DEFAULT_FOLDER)
    processed_table.to_csv(output_file, sep = delimiter)