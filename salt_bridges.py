#!/usr/bin/env python

"""
Measure inter-chain salt-bridges
"""

import numpy as np
import os
import csv
import gpcr_variables
from gpcr_variables import WEINSTEIN_NUMBERING_ORIGINAL, RECEPTORS_LIST, PROCESSED_RESULTS_FOLDER, SB_DISTANCE

__author__ = "A.J. Preto"
__email__ = "martinsgomes.jose@gmail.com"
__group__ = "Data-Driven Molecular Design"
__group_leader__ = "Irina S. Moreira"
__project__ = "GPCRs"

def csv_writer(dataset,name,res_names_A,res_names_B):

    """
    Write the output .csv file, in the pipeline, generates a different file for 
    each corresponding .pdb file
    """
    csv_name = name + '.csv'
    with open(csv_name, 'w') as csvfile:
        csv_writer = csv.writer(csvfile, delimiter=';')       
        res_names_B.insert(0,' ')
        csv_writer.writerow(res_names_B)
        count = 0
        for row in dataset:
            convert = []
            convert.append(res_names_A[count])
            for element in row:
                convert.append(element)
            csv_writer.writerow(convert)
            count = count + 1

def dist(input_coords):

    return np.sqrt(np.sum(input_coords * input_coords))

def calc_residue_dist(residue_one, residue_two) :

    """
    Returns the atomic distance between two residues oxygen or nitrogens
    LYS -> NZ
    ASP -> OD1, OD2
    GLU -> OE1, OE2
    ARG -> NH1, NH2
    """
    try:
        if (residue_one.resname == "LYS") and (residue_two.resname == "ASP"):
            diff_vector_1  = residue_one["NZ"].coord - residue_two["OD1"].coord
            diff_vector_2  = residue_one["NZ"].coord - residue_two["OD2"].coord
            diff_vector = [dist(diff_vector_1), dist(diff_vector_2)]
        if (residue_one.resname == "ARG") and (residue_two.resname == "ASP"):
            diff_vector_1  = residue_one["NH1"].coord - residue_two["OD1"].coord
            diff_vector_2  = residue_one["NH2"].coord - residue_two["OD1"].coord
            diff_vector_3  = residue_one["NH1"].coord - residue_two["OD2"].coord
            diff_vector_4  = residue_one["NH2"].coord - residue_two["OD2"].coord
            diff_vector = [dist(diff_vector_1), dist(diff_vector_2), dist(diff_vector_3), dist(diff_vector_4)]
        if (residue_one.resname == "LYS") and (residue_two.resname == "GLU"):
            diff_vector_1  = residue_one["NZ"].coord - residue_two["OE1"].coord
            diff_vector_2  = residue_one["NZ"].coord - residue_two["OE2"].coord
            diff_vector = [dist(diff_vector_1), dist(diff_vector_2)]
        if (residue_one.resname == "ARG") and (residue_two.resname == "GLU"):
            diff_vector_1  = residue_one["NH1"].coord - residue_two["OE1"].coord
            diff_vector_2  = residue_one["NH2"].coord - residue_two["OE1"].coord
            diff_vector_3  = residue_one["NH1"].coord - residue_two["OE2"].coord
            diff_vector_4  = residue_one["NH2"].coord - residue_two["OE2"].coord
            diff_vector = [dist(diff_vector_1), dist(diff_vector_2), dist(diff_vector_3), dist(diff_vector_4)]
        if (residue_one.resname == "GLU") and (residue_two.resname == "LYS"):
            diff_vector_1  = residue_one["OE1"].coord - residue_two["NZ"].coord
            diff_vector_2  = residue_one["OE2"].coord - residue_two["NZ"].coord
            diff_vector = [dist(diff_vector_1), dist(diff_vector_2)]
        if (residue_one.resname == "ASP") and (residue_two.resname == "LYS"):
            diff_vector_1  = residue_one["OD1"].coord - residue_two["NZ"].coord
            diff_vector_2  = residue_one["OD2"].coord - residue_two["NZ"].coord
            diff_vector = [dist(diff_vector_1), dist(diff_vector_2)]
        if (residue_one.resname == "GLU") and (residue_two.resname == "ARG"):
            diff_vector_1  = residue_one["OE1"].coord - residue_two["NH1"].coord
            diff_vector_2  = residue_one["OE1"].coord - residue_two["NH2"].coord
            diff_vector_3  = residue_one["OE2"].coord - residue_two["NH1"].coord
            diff_vector_4  = residue_one["OE2"].coord - residue_two["NH2"].coord
            diff_vector = [dist(diff_vector_1), dist(diff_vector_2), dist(diff_vector_3), dist(diff_vector_4)]
        if (residue_one.resname == "ASP") and (residue_two.resname == "ARG"):
            diff_vector_1  = residue_one["OD1"].coord - residue_two["NH1"].coord
            diff_vector_2  = residue_one["OD1"].coord - residue_two["NH2"].coord
            diff_vector_3  = residue_one["OD2"].coord - residue_two["NH1"].coord
            diff_vector_4  = residue_one["OD2"].coord - residue_two["NH2"].coord
            diff_vector = [dist(diff_vector_1), dist(diff_vector_2), dist(diff_vector_3), dist(diff_vector_4)]
        return diff_vector
    except:
        return 0.0

def calc_dist_matrix(chain_one, chain_two):

    """
    Returns a matrix of C-alpha distances between two chains
    """
    output_table = []
    for row, residue_one in enumerate(chain_one):
        current_row = []
        for col, residue_two in enumerate(chain_two):
            distance = calc_residue_dist(residue_one, residue_two)
            if type(distance) == list:
                distance = ",".join([str(value) for value in distance])
            current_row.append(distance)
        output_table.append(current_row)
    return output_table

def get_residues(model,chain,model_name):

    """
    Transform the input receptor sequence into
    weinstein numbered
    """
    RES=[]
    res_number = 0
    if chain == 'A':
        for residue in model[chain].get_residues():
            res_name =residue.get_resname()
            res_number = res_number + 1
            for dr in DXR_tags:
                try:
                    if dr in model_name:
                        weinstein_res_number = x_50(DXR_dict[dr],res_number)
                        res_info =str(weinstein_res_number) + res_name
                        RES.append(res_info)
                except:
                    res_info = str(res_number) + res_name
                    RES.append(res_info)
    else:
        for residue in model[chain].get_residues():
            res_name =residue.get_resname()
            res_number = res_number + 1
            res_info = str(res_number) + res_name
            RES.append(res_info)
    return RES

def process_csv(input_file, delimiter = ";"):

    """
    Process the input .csv weinstein numbering file
    """
    opened_file = open(input_file, "r").readlines()
    for row in opened_file:
        row = row.replace("\n","")
        new_row = row.split(delimiter)
        yield new_row

def process_template(input_csv, receptor_list):

    """
    Generate a dictionary in which each key is the receptor
    and the value is a list of lists of the substructure weinstein
     numberings
    """
    opened_file = list(process_csv(input_csv))
    output_dict = {}
    receptor_col = 0
    TM_col = 1
    receptor_holder = ""
    for row in opened_file[1:]:
        if row[receptor_col] != receptor_holder:
            receptor_holder = row[receptor_col] 
            output_dict[receptor_holder] = [row[1:]]
        else:
            output_dict[receptor_holder].append(row[1:])
    return output_dict

def x_50(DXR, res_number):

    """
    Renumber the substructures according to their x50 residues
    """
    count = 0
    for struct_element in DXR:
        count = count + 1 
        if res_number in range(int(struct_element[1]),int(struct_element[2])):
            weinstein_count = count
            if res_number == int(struct_element[3]):
                weinstein = 50
            elif res_number < int(struct_element[3]):
                weinstein = 50 - (int(struct_element[3]) - res_number)
            elif res_number > int(struct_element[3]):
                weinstein = 50 + (res_number - int(struct_element[3]))
            break
    full_weinstein = str(weinstein_count) + '.' + str(weinstein)
    return full_weinstein

def write_distance_table(contact_distance = 4.0):

    """
    Write the output table correspondant to each .pdb file in the folder.
    Change "contact_distance" to your desired distance (in angstroms)
    """    
    import Bio.PDB
    for files in os.listdir('.'):
        if files.endswith('.pdb'):
            name = PROCESSED_RESULTS_FOLDER + "/novel_sb_" + files[0:-4]
            structure = Bio.PDB.PDBParser().get_structure(name,files)
            model = structure[0]
            res_names_A = get_residues(model,"A",files)
            res_names_B = get_residues(model,"B",files)
            dist_matrix = calc_dist_matrix(model["A"], model["B"])
            csv_writer(dist_matrix,name,res_names_A,res_names_B)

"""
Variable Initialization
"""
weinstein_template = WEINSTEIN_NUMBERING_ORIGINAL
DXR_tags = RECEPTORS_LIST
DXR_dict = process_template(weinstein_template, DXR_tags)
write_distance_table(contact_distance = SB_DISTANCE)



