#!/usr/bin/env python

"""
Write summary tables and overall table
"""

import os
import csv
import raw_parser
import pandas as pd
import gpcr_variables
from gpcr_variables import HBOND_VAR_GPCR, COCOMAPS_MOL_VAR, FINAL_TABLE, \
                            RESULTS_FOLDER, PROCESSED_RESULTS_FOLDER, DEFAULT_FOLDER, \
                            STRUCTURAL_FEATURES, STRUCTURAL_FEATURES, INTERPROSURF_COMPLEX_VAR, \
                            INTERPROSURF_START, INTERPROSURF_RESIDUES_VAR, DECIMAL_HOUSES, \
                            HB_TOTAL_PREFIX, SB_PREFIX, WEINSTEIN_NUMBERING_ORIGINAL, \
                            SUBSTRUCTURES_EVALUATED, COCOMAPS_START, SUMMARY_FOLDER, \
                            SUMMARY_TABLE, COCOMAPS_MOL_VAR_2, STRUCTURAL_PDBS_FOLDER               
__author__ = "A.J. Preto"
__email__ = "martinsgomes.jose@gmail.com"
__group__ = "Data-Driven Molecular Design"
__group_leader__ = "Irina S. Moreira"
__project__ = "GPCRs"

def calculate_interface(input_row):

    """
    Calculate interface values based on receptor, partner and total values
    """
    difference = float(input_row[0]) + float(input_row[1]) - float(input_row[2])
    return round(abs(difference),DECIMAL_HOUSES)

def generate_dict(input_list):

    """
    Generate a dictionary take list values as keys
    and starting at 0
    """
    output_dict = {}
    for residue in input_list:
        output_dict[residue] = 0
    return output_dict

def calculate_percentage(input_dict, input_total):

    percentages_dict = {}
    for key in input_dict:
        percentages_dict[key] = round(float(float(input_dict[key]) * 100.0) / float(input_total), DECIMAL_HOUSES)
    return percentages_dict

def evaluate_content(input_list):

    """
    Count amino acid in polar associated groups
    """
    from gpcr_variables import NONPOLAR_ALIPHATIC, POLAR_UNCHARGED, BASIC_CHARGED, \
                                ACIDIC_CHARGED, NONPOLAR_AROMATIC, AMINOACID_LIST, \
                                GROUPS_DICT, POLAR_NAMES

    current_groups_dict = generate_dict(POLAR_NAMES)
    current_amino_acid_dict = generate_dict(AMINOACID_LIST)
    total = len(input_list)
    for residue in input_list:
        for key in GROUPS_DICT.keys():
            for group_residue in GROUPS_DICT[key]:
                if residue == group_residue:
                    current_groups_dict[key] += 1
                    current_amino_acid_dict[residue] += 1
    return current_groups_dict, current_amino_acid_dict, total

def generate_interface_dict(input_table, chain = "A"):

    """
    Generate the processed dictionaries of amino acid
    interface group percentages
    """
    receptor_interface_list = input_table.loc[chain]["Three Letter Code"].values.tolist()
    chemical_groups, individual_frequencies, interface_residues_number = evaluate_content(receptor_interface_list)
    percentage_chemical_groups, percentage_individual_aa = calculate_percentage(chemical_groups, interface_residues_number), \
                                                            calculate_percentage(individual_frequencies, interface_residues_number)
    return percentage_chemical_groups, percentage_individual_aa

def single_substructure_interprosurf(input_column, substructure_values):

    """
    Evaluate one substructure and yield the sum of the ASA difference
    """
    residues = [x for x in range(int(substructure_values[0]), int(substructure_values[1]))]
    output_values = []
    for current_residue in residues:
        try:
            difference_value = input_column.loc[current_residue][0]
            output_values.append(difference_value)
        except:
            output_values.append(0)
    return round(sum(output_values), DECIMAL_HOUSES)

def evaluate_substructures_interprosurf(input_table, input_receptor):

    """
    Fetch the Weinstein intervals and yield the ASA difference sum for each substructure
    """
    weinstein_table = raw_parser.retrieve_clean_weinstein(WEINSTEIN_NUMBERING_ORIGINAL, delimiter = ";", with_x = False)
    ICL1, ICL2, ICL3, HX8 = [weinstein_table[input_receptor + "_TM1"][1],weinstein_table[input_receptor + "_TM2"][0]], \
                            [weinstein_table[input_receptor + "_TM3"][1],weinstein_table[input_receptor + "_TM4"][0]], \
                            [weinstein_table[input_receptor + "_TM5"][1],weinstein_table[input_receptor + "_TM7"][0]], \
                            weinstein_table[input_receptor + "_H8"]
    usable_table = input_table.loc["A"]
    res_number_col = usable_table[["Residue Number","Difference"]].set_index("Residue Number")
    ASA_ICL1, ASA_ICL2, ASA_ICL3, ASA_HX8 = single_substructure_interprosurf(res_number_col, ICL1), \
                                            single_substructure_interprosurf(res_number_col, ICL2), \
                                            single_substructure_interprosurf(res_number_col, ICL3), \
                                            single_substructure_interprosurf(res_number_col, HX8)
    return ASA_ICL1, ASA_ICL2, ASA_ICL3, ASA_HX8

def single_substructure_cocomaps(input_columns, substructure_values):

    """
    Evaluate one substructure and yield the sum of the ASA difference with cocomaps
    """
    residues = [x for x in range(int(substructure_values[0]), int(substructure_values[1]))]
    output_values = []
    for current_residue in residues:
        for row in input_columns:
            if int(row[0]) == int(current_residue):
                output_values.append(float(row[1]))
            else:
                output_values.append(0)
    return round(sum(output_values), DECIMAL_HOUSES)

def evaluate_substructures_cocomaps(input_table, input_receptor):

    """
    Fetch the Weinstein intervals and yield the ASA difference sum for each substructure using cocomaps
    """
    weinstein_table = raw_parser.retrieve_clean_weinstein(WEINSTEIN_NUMBERING_ORIGINAL, delimiter = ";", with_x = False)
    ICL1, ICL2, ICL3, HX8 = [weinstein_table[input_receptor + "_TM1"][1],weinstein_table[input_receptor + "_TM2"][0]], \
                            [weinstein_table[input_receptor + "_TM3"][1],weinstein_table[input_receptor + "_TM4"][0]], \
                            [weinstein_table[input_receptor + "_TM5"][1],weinstein_table[input_receptor + "_TM7"][0]], \
                            weinstein_table[input_receptor + "_H8"]
    ASA_ICL1, ASA_ICL2, ASA_ICL3, ASA_HX8 = single_substructure_cocomaps(input_table, ICL1), \
                                            single_substructure_cocomaps(input_table, ICL2), \
                                            single_substructure_cocomaps(input_table, ICL3), \
                                            single_substructure_cocomaps(input_table, HX8)
    return ASA_ICL1, ASA_ICL2, ASA_ICL3, ASA_HX8

def single_substructure_hbonds(input_column, substructure_values):

    """
    Evaluate one substructure and yield the hydrogen bond count
    """
    output_value = 0
    values = input_column[int(substructure_values[0]):int(substructure_values[1])]
    for row in values:
        if row == "-":
            continue
        elif len(row.split(",")) > 1:
            output_value += len(row.split(","))
        else:
            output_value += 1
    return output_value

def evaluate_substructures_hbonds(input_row, input_receptor):

    """
    Fetch the Weinstein intervals and yield the hydrogen bond count per substructure
    """
    weinstein_table = raw_parser.retrieve_clean_weinstein(WEINSTEIN_NUMBERING_ORIGINAL, delimiter = ";", with_x = False)
    ICL1, ICL2, ICL3, HX8 = [weinstein_table[input_receptor + "_TM1"][1],weinstein_table[input_receptor + "_TM2"][0]], \
                            [weinstein_table[input_receptor + "_TM3"][1],weinstein_table[input_receptor + "_TM4"][0]], \
                            [weinstein_table[input_receptor + "_TM5"][1],weinstein_table[input_receptor + "_TM7"][0]], \
                            weinstein_table[input_receptor + "_H8"]
    HB_ICL1, HB_ICL2, HB_ICL3, HB_HX8 = single_substructure_hbonds(input_row, ICL1), \
                                            single_substructure_hbonds(input_row, ICL2), \
                                            single_substructure_hbonds(input_row, ICL3), \
                                            single_substructure_hbonds(input_row, HX8)
    return HB_ICL1, HB_ICL2, HB_ICL3, HB_HX8

def single_substructure_sbridges(input_table, substructure_values, delimiter = ";"):

    """
    Evaluate one substructure and yield the salt bridges count
    """
    from gpcr_variables import SB_DISTANCE
    output_value = 0
    values = input_table[int(substructure_values[0]):int(substructure_values[1])]
    for row in values:
        row = row.split(delimiter)
        for cell in row[1:]:
            if len(cell.split(",")) > 1:
                for bond in cell.split(","):
                    if float(bond) < SB_DISTANCE:
                        output_value += 1
            elif int(float(cell)) != 0 and float(cell < SB_DISTANCE):
                output_value += 1
    return output_value

def evaluate_substructures_sbridges(input_file, input_receptor, delimiter = ","):

    """
    Open Salt-Bridges associated file and count the bonds associated to each substructure
    """
    import pandas as pd

    opened_file = open(input_file, "r").readlines()
    weinstein_table = raw_parser.retrieve_clean_weinstein(WEINSTEIN_NUMBERING_ORIGINAL, delimiter = ";", with_x = False)
    ICL1, ICL2, ICL3, HX8 = [weinstein_table[input_receptor + "_TM1"][1],weinstein_table[input_receptor + "_TM2"][0]], \
                            [weinstein_table[input_receptor + "_TM3"][1],weinstein_table[input_receptor + "_TM4"][0]], \
                            [weinstein_table[input_receptor + "_TM5"][1],weinstein_table[input_receptor + "_TM7"][0]], \
                            weinstein_table[input_receptor + "_H8"]
    
    SB_ICL1, SB_ICL2, SB_ICL3, SB_HX8 = single_substructure_sbridges(opened_file, ICL1), \
                                            single_substructure_sbridges(opened_file, ICL2), \
                                            single_substructure_sbridges(opened_file, ICL3), \
                                            single_substructure_sbridges(opened_file, HX8)
    return SB_ICL1, SB_ICL2, SB_ICL3, SB_HX8

def write_summary_table(receptor, partner, delimiter = ","):

    """
    Write a summary table for each complex
    """

    output_name = SUMMARY_FOLDER + "/" + FINAL_TABLE + "_" + receptor + "_" + partner + ".csv"
    with open(output_name, "w") as final_table:
        
        """
        Complex associated values
        """
        complex_header = "Complex" + delimiter + receptor + "-" + partner + "\n"
        final_table.write(complex_header)
        target_complex_name = RESULTS_FOLDER + "/" + INTERPROSURF_START + "_" + receptor + "-" + partner + "_" + INTERPROSURF_COMPLEX_VAR + ".csv"
        intersurf_table = raw_parser.intersurf_parser(target_complex_name)
        count = 1
        for row_index, row in zip(STRUCTURAL_FEATURES, intersurf_table):
            if count < 4:
                complex_row = [row_index] + [str(calculate_interface(row))]
            else:
                complex_row = [row_index] + [str(row[-1])]
            final_complex_row = delimiter.join(complex_row) + "\n"
            count += 1
            final_table.write(final_complex_row)

        """
        Partner associated values
        """
        cocomaps_partner = RESULTS_FOLDER + "/" + COCOMAPS_START + "_" + receptor + "-" + partner + "_" + COCOMAPS_MOL_VAR_2 + ".txt"
        cocomaps_table_partner = raw_parser.cocomaps_parser(cocomaps_partner)
        target_residues_name = RESULTS_FOLDER + "/" + INTERPROSURF_START + "_" + receptor + "-" + partner + "_" + INTERPROSURF_RESIDUES_VAR + ".csv"
        opened_table = pd.read_csv(target_residues_name, index_col = [5])

        partner_header = "Partner" + delimiter + partner + "\n"
        final_table.write(partner_header)

        interface_residues_partner = "total interface residues:" + delimiter + str(len(cocomaps_table_partner)) + '\n'
        final_table.write(interface_residues_partner)
        partner_group_percentage, partner_aa_percentage = generate_interface_dict(opened_table, chain = "B")

        for partner_key in partner_aa_percentage:
            partner_aa_row = partner_key + delimiter + str(partner_aa_percentage[partner_key]) + "\n"
            final_table.write(partner_aa_row)
        
        for partner_key_group in partner_group_percentage:
            partner_group_row = partner_key_group + delimiter + str(partner_group_percentage[partner_key_group]) + "\n"
            final_table.write(partner_group_row)   

        """
        Receptor associated values
        """
        cocomaps_receptor = RESULTS_FOLDER + "/" + COCOMAPS_START + "_" + receptor + "-" + partner + "_" + COCOMAPS_MOL_VAR + ".txt"
        cocomaps_table = raw_parser.cocomaps_parser(cocomaps_receptor)
        cocomaps_substructures_values = evaluate_substructures_cocomaps(cocomaps_table, receptor)
        receptor_header = "Receptor" + delimiter + receptor + "\n"
        final_table.write(receptor_header)
        interface_residues = "total interface residues:" + delimiter + str(len(cocomaps_table)) + '\n'
        final_table.write(interface_residues)
        receptor_group_percentage, receptor_aa_percentage = generate_interface_dict(opened_table)
        for receptor_key in receptor_aa_percentage:
            receptor_aa_row = receptor_key + delimiter + str(receptor_aa_percentage[receptor_key]) + "\n"
            final_table.write(receptor_aa_row)
        
        for receptor_key_group in receptor_group_percentage:
            receptor_group_row = receptor_key_group + delimiter + str(receptor_group_percentage[receptor_key_group]) + "\n"
            final_table.write(receptor_group_row)

        """
        Total hydrogen bonds and salt bridges
        """
        HB_name =  PROCESSED_RESULTS_FOLDER + "/" + receptor + "_" + partner + "_" + HB_TOTAL_PREFIX + ".csv"
        hb_total = raw_parser.hb_total_counter(HB_name, delimiter = ";")
        SB_name =  PROCESSED_RESULTS_FOLDER + "/" + SB_PREFIX + "_" + receptor + "-" + partner + ".csv"
        sb_total = raw_parser.sb_total_counter(SB_name, delimiter = ";")
        hb_sb = hb_total + sb_total
        hb_sb_row = "total HB/SB" + delimiter + str(hb_sb) + "\n"
        final_table.write(hb_sb_row)

        """
        Interprosurf ASA values
        """
        interprosurf_substructures_values = evaluate_substructures_interprosurf(opened_table, receptor)
        for substructure_ASA, inter_sub in zip(SUBSTRUCTURES_EVALUATED, interprosurf_substructures_values):
            substructure_row = "Intersurf ASA: " + substructure_ASA + delimiter + str(inter_sub) + "\n"
            final_table.write(substructure_row)

        """
        Cocomaps ASA values
        """
        for alternative_ASA, cocomaps_sub in zip(SUBSTRUCTURES_EVALUATED, cocomaps_substructures_values):
            substructure_row_cocomaps = "Cocomaps ASA: " + alternative_ASA + delimiter + str(cocomaps_sub) + "\n"
            final_table.write(substructure_row_cocomaps)

        """
        Hydrogen bonds and salt bridges by substructure
        """
        hbonds_name = PROCESSED_RESULTS_FOLDER + "/" + receptor + "_" + partner + "_" + HBOND_VAR_GPCR + ".csv"
        hbonds_var = raw_parser.parse_hb_compare(hbonds_name, delimiter = ";")
        hbonds_values = evaluate_substructures_hbonds(hbonds_var, receptor)
        sbridges_values = evaluate_substructures_sbridges(SB_name, receptor)
        hb_sb_values = [int(x) + int(y) for x, y in zip(hbonds_values, sbridges_values)]
        for substructure_hb_sb_name, hb_sb_val in zip(SUBSTRUCTURES_EVALUATED, hb_sb_values):
            substructure_row_hb_sb = "HB/SB: " + substructure_hb_sb_name + delimiter + str(hb_sb_val) + "\n"
            final_table.write(substructure_row_hb_sb)

"""
Apply over ".pdb" files on folder
"""
failed = []
for files in os.listdir(STRUCTURAL_PDBS_FOLDER):
    if files.endswith(".pdb"):
        file_loc = DEFAULT_FOLDER + "/" + STRUCTURAL_PDBS_FOLDER + "/" + files  
        split_name = files.split("-")
        receptor, partner = split_name[0], split_name[1][0:-4]
        write_summary_table(receptor, partner)

"""
Join all the individual summary tables
"""
raw_parser.join_tables(table_type = FINAL_TABLE)
print("The following PDBS failed:", failed)