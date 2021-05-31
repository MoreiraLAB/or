#!/usr/bin/env python

"""
Correct the partner numbering for interaction plot construction
"""


import pandas as pd 
import os
import csv
import sys
opened_corrector = pd.read_csv("templates/partners_alternative_numbering.csv", sep = ";", header = 0)
CONVERTER_AMINO_ACIDS_THREE_TO_ONE = {"A":"ALA", "R":"ARG", "N": "ASN", "D": "ASP", \
                                        "C":"CYS", "B": "ASX", "E": "GLU", "Q": "GLN", \
                                        "Z": "GLX", "G": "GLY", "H": "HIS", "I": "ILE", \
                                        "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE", \
                                        "P": "PRO", "S": "SER", "T": "THR", "W": "TRP", \
                                        "Y": "TYR", "V": "VAL"}
old_numbering = opened_corrector.apply(lambda x: x['old_numbering'][:-1], axis = 1)
amino_type = opened_corrector.apply(lambda x: CONVERTER_AMINO_ACIDS_THREE_TO_ONE[x["old_numbering"][-1]], axis = 1)
opened_corrector["old_numbering"] = old_numbering
opened_corrector["amino_type"] = amino_type

for files in os.listdir("processed_results"):
    file_loc =  "processed_results/" + files
    split_name = files.split("_")
    if len(split_name) >= 3:
        if [split_name[0],split_name[1],split_name[2]] == ["weinstein","inter","chain"]:
            try:
                partner = split_name[-2].split("-")[1]
            except:
                partner = split_name[-1].split("-")[1].split(".")[0]
            opened_file = pd.read_csv(file_loc, sep = ";", header = 0)
            rename_dict = {}
            for current_col in list(opened_file):
                only_number = current_col[0:-3]
                try:
                    new_name = opened_corrector.loc[(opened_corrector["old_numbering"] == only_number) & \
                                                    (opened_corrector["partner"] == partner)]["new_numbering"]
                    if not new_name.empty:
                        rename_dict[current_col] = new_name.values[0] + opened_corrector.loc[(opened_corrector["old_numbering"] == only_number) & \
                                                    (opened_corrector["partner"] == partner)]["amino_type"].values[0]
                except:
                    continue
            output_df = opened_file.rename(columns = rename_dict)
            output_loc = "processed_results/adapted_partner_" + files
            with open(output_loc, 'w') as csvfile:
                csv_writer = csv.writer(csvfile, delimiter=';')
                csv_writer.writerow(list(output_df))
                for index, row in output_df.iterrows():
                    csv_writer.writerow(list(row))

            #output_df.to_csv(output_loc, index = False, sep = ";", line_terminator='\n\n')
            print("Writing:", output_loc)
