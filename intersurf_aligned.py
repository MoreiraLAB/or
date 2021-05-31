#!/usr/bin/env python

"""
Process Intersurf output
"""

import os
import csv
import raw_parser
import gpcr_variables
from gpcr_variables import DEFAULT_FOLDER, RESULTS_FOLDER, PROCESSED_RESULTS_FOLDER, ALIGN_VECTORS, DIMERS, INTERPROSURF_START

__author__ = "A.J. Preto"
__email__ = "martinsgomes.jose@gmail.com"
__group__ = "Data-Driven Molecular Design"
__group_leader__ = "Irina S. Moreira"
__project__ = "GPCRs"

def build_csv(input_file,seq,out_name,chain):

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
            spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
            count = 0
            output_list.append('-')
            for row in spamreader:
                if count != 0:
                    if (int(row[1])) == (new_res_number) and (row[5]==chain):
                        output_list[-1] = row[-1]
                else:
                    count = count + 1
    new_file = open(out_name,'w')
    writer=csv.writer(new_file,delimiter = ';')
    writer.writerow(output_list)
    new_file.close()   

def apply_on_folder(input_alignments, input_dimers):

    """
    Deploy the processing pipeline, pay attention to the folder structure
    """
    for receptor in input_dimers.keys():
        for partner in input_dimers[receptor]:
            intersurf_name = results_folder  + "/" + INTERPROSURF_START + "_" + receptor + "-" + partner + "_residues.csv"
            output_name_DR = processed_results_folder + "/" + receptor + "_" + partner + "_toDR_amino.csv"
            output_name_partner = processed_results_folder + "/" + receptor + "_" + partner + "_toprotein_amino.csv"
            try:
                build_csv(intersurf_name, input_alignments[receptor], output_name_DR, "A")
                build_csv(intersurf_name, input_alignments[partner], output_name_partner, "B")
            except:
                print("Failed to perform action on:", receptor, partner)

"""
Variable initialization
"""     
home = DEFAULT_FOLDER
results_folder = home + "/" + RESULTS_FOLDER
processed_results_folder = home + "/" + PROCESSED_RESULTS_FOLDER
aligned_proteins = raw_parser.retrieve_clean(ALIGN_VECTORS)
dimers = raw_parser.retrieve_clean(DIMERS)
apply_on_folder(aligned_proteins, dimers)