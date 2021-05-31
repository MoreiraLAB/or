#!/usr/bin/env python

"""
Measure alpha-Carbon interchain distance
"""

import numpy as np
import os
import csv
from gpcr_variables import WEINSTEIN_NUMBERING_ORIGINAL, RECEPTORS_LIST,\
                            PROCESSED_RESULTS_FOLDER, HB_DISTANCE, INTER_CHAIN_CA, \
                            GPCRDB_ALIGNMENT, GPROTEIN_SUBSTRUCTURES, \
                            STRUCTURAL_PDBS_FOLDER

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
            
def calc_residue_dist(residue_one, residue_two) :

    """
    Returns the C-alpha distance between two residues
    """
    diff_vector  = residue_one["CA"].coord - residue_two["CA"].coord
    return np.sqrt(np.sum(diff_vector * diff_vector))

def calc_dist_matrix(chain_one, chain_two):

    """
    Returns a matrix of C-alpha distances between two chains
    """
    answer = np.zeros((len(chain_one), len(chain_two)), np.float)
    for row, residue_one in enumerate(chain_one) :
        for col, residue_two in enumerate(chain_two) :
            answer[row, col] = calc_residue_dist(residue_one, residue_two)
    return answer

def process_gprot(input_gprotein, input_template = GPROTEIN_SUBSTRUCTURES):

    """
    Process the g-protein substructures numbering and match with original numbering
    """
    opened_file = open(GPCRDB_ALIGNMENT, "r").readlines()
    residue_count = 1
    output_dict = {}
    for row in opened_file:
        row = row.replace("\n","").split(".")
        if row[0] == input_gprotein:
            for sub_in_temp, sub in zip(row[1:], GPROTEIN_SUBSTRUCTURES):
                substructure_seq = list(filter(None, sub_in_temp.split(",")))
                substructure_count = 1
                for residue in substructure_seq:
                    if residue != "-":
                        final_name = sub + "." + str(substructure_count)
                        output_dict[residue_count] = final_name
                        residue_count += 1
                    substructure_count += 1
    return output_dict

def get_residues(model,chain,model_name, start_dict = {}):

    """
    Transform the input receptor sequence into
    weinstein numbered
    """
    RES = []
    partner_name = model_name.split("-")[1][0:-4]
    if partner_name not in start_dict.keys(): res_number = 0
    else: res_number = start_dict[partner_name]
    #if chain == "A": res_number = start_dict[model_name.split("-")[0]]

    if chain == 'A':
        for residue in model[chain].get_residues():
            res_name = residue.get_resname()
            res_number = res_number + 1
            for dr in DXR_tags:
                try:
                    if dr in model_name:
                        weinstein_res_number = x_50(DXR_dict[dr],res_number)
                        res_info = str(weinstein_res_number) + res_name
                        RES.append(res_info)
                except:
                    res_info = str(res_number + start_dict[model_name.split("-")[0]] - 1) + res_name
                    RES.append(res_info)
    else:
        for residue in model[chain].get_residues():
            res_number = res_number + 1
            res_name =residue.get_resname()
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

def write_distance_table(contact_distance = 5.0, target_folder = STRUCTURAL_PDBS_FOLDER):

    """
    Write the output table correspondant to each .pdb file in the folder.
    Change "contact_distance" to your desired distance (in angstroms)
    """    
    import Bio.PDB
    for files in os.listdir(target_folder):
        if files.endswith('.pdb'):
            name = PROCESSED_RESULTS_FOLDER + "/weinstein_" + INTER_CHAIN_CA + "_" + files[0:-4]
            file_loc = target_folder + "/" + files
            structure = Bio.PDB.PDBParser().get_structure(name, file_loc)
            model = structure[0]
            res_names_A = get_residues(model,"A", files, start_dict = {"DOR": 34, "KOR": 52, "MOR": 62, "NOP": 43})
            res_names_B = get_residues(model,"B", files, start_dict = {"Gq_6OIJ": 6, "Gq_6DDF": 6, "Gssh": 8, "Gslo": 8, \
                                                                        "Arr2_6PWC": 6, "Arr3_6PWC": 7,"Arr2_6U1N": 5, "Arr3_6U1N": 6, \
                                                                        "G11_6OIJ": 6, "G14_6OIJ": 4, "G14_6DDF": 4, \
                                                                        "G15_6OIJ": 9, "G12_6DDF": 25, "G12_6OIJ": 25, \
                                                                        "G11_6DDF": 6, "G15_6DDF": 9})
            dist_matrix = calc_dist_matrix(model["A"], model["B"])
            contact_map = dist_matrix < contact_distance
            csv_writer(dist_matrix,name, res_names_A, res_names_B)
            

"""
Variable Initialization
"""
weinstein_template = WEINSTEIN_NUMBERING_ORIGINAL
DXR_tags = RECEPTORS_LIST
DXR_dict = process_template(weinstein_template, DXR_tags)
write_distance_table(contact_distance = HB_DISTANCE)



