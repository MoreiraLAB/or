#!/usr/bin/env python

"""
Process COCOMAPS output
"""

import os
import csv
import raw_parser
import gpcr_variables
from gpcr_variables import DEFAULT_FOLDER, RESULTS_FOLDER, \
                            PROCESSED_RESULTS_FOLDER, ALIGN_VECTORS, DIMERS

__author__ = "A.J. Preto"
__email__ = "martinsgomes.jose@gmail.com"
__group__ = "Data-Driven Molecular Design"
__group_leader__ = "Irina S. Moreira"
__project__ = "GPCRs"

def build_csv(input_file,seq,out_name):

    """
    Build the output .csv files
    """
    output_list = []
    corrected_number = 0
    for res_number in range(1,len(seq) + 1):
        if seq[res_number - 1] == '-':
            corrected_number = corrected_number + 1
        new_res_number = res_number - corrected_number
        with open(input_file, 'r') as csvfile:
            csv_reader = csv.reader(csvfile, delimiter=' ', quotechar='|')
            count = 0
            output_list.append('-')
            for row in csv_reader:
                if count != 0:
                    if (int(row[1])) == (new_res_number):
                        output_list[-1] = row[-1]
                else:
                    count = count + 1
    new_file = open(out_name,'w')
    writer = csv.writer(new_file,delimiter = ';')
    writer.writerow(output_list)
    new_file.close()   

def apply_on_folder(input_alignments, input_dimers):

    """
    Deploy the processing pipeline, pay attention to the folder structure
    """
    for receptor in input_dimers.keys():
        for partner in input_dimers[receptor]:
            cocomaps_name_DR = results_folder  + "/COCOMAPS_" + receptor + "-" + partner + "_asa_mol1.txt"
            output_name_DR = processed_results_folder  + "/" + receptor + "_" + partner + "_toDR_cocomaps.csv"
            cocomaps_name_partner = results_folder  + "/COCOMAPS_" + receptor + "-" + partner + "_asa_mol2.txt"
            output_name_partner = processed_results_folder  + "/" + receptor + "_" + partner + "_toprotein_cocomaps.csv"
            build_csv(cocomaps_name_DR, input_alignments[receptor], output_name_DR)
            build_csv(cocomaps_name_partner, input_alignments[partner], output_name_partner)
"""
Variable initialization
"""

home = DEFAULT_FOLDER
results_folder = home + "/" + RESULTS_FOLDER
processed_results_folder = home + "/" + PROCESSED_RESULTS_FOLDER
aligned_proteins = raw_parser.retrieve_clean(ALIGN_VECTORS)
input_dimers = raw_parser.retrieve_clean(DIMERS)
apply_on_folder(aligned_proteins, input_dimers)
