#!/usr/bin/env python

"""
Analyse interface output from COCOMAPS to yield
residue type and group information
"""

import os
import csv
import raw_parser
import gpcr_variables
from gpcr_variables import DEFAULT_FOLDER, ALIGN_VECTORS, DIMERS, PROCESSED_RESULTS_FOLDER

__author__ = "A.J. Preto"
__email__ = "martinsgomes.jose@gmail.com"
__group__ = "Data-Driven Molecular Design"
__group_leader__ = "Irina S. Moreira"
__project__ = "GPCRs"

def total_counter(input_seq_list,to_count):

    """
    Count all non-gaps in the alignment
    """
    count = 0
    for res in input_seq_list:
        if res != to_count:
            count = count + 1
    return count

def counter(input_seq_list,to_count):

    """
    Count a specific type of residue in a sequence
    """
    count = 0
    for res in input_seq_list:
        if res == to_count:
            count = count + 1
    return count

def res_type_counter(input_seq_list):

    """
    Count residues according to their residue type
    """
    res_name = {'A':'Alanine', 'R':'Arginine','D':'Aspartate',
                'N':'Asparagine','C':'Cysteine','E':'Glutamate','G':'Glycine',
                'Q':'Glutamine','H':'Histidine','I':'Isoleucine',
                'L':'Leucine','K':'Lysine','M':'Methionine',
                'F':'Phenylalanine','P':'Proline',
                'S':'Serine','T':'Threonine','W':'Tryptophan',
                'Y':'Tyrosine','V':'Valine'}
    res_type = {'A':0, 'R':0,'D':0,'N':0,'C':0,'E':0,'Q':0,'G':0,'H':0,'I':0,'L':0,'K':0,'M':0,'F':0,'P':0,'S':0,'T':0,'W':0,'Y':0,'V':0}
    res_list = ['G','A','V','L','I','M','S','T','C','P','N','Q','K','R','H','D','E','F','Y','W']
    new_dict = {}
    for target_res in res_list:
        res_type[target_res] = counter(input_seq_list,target_res)
    for rowing in res_list:
        new_dict[res_name[rowing]] = str((float(res_type[rowing])/total_counter(input_seq_list,'-'))*100)
    return new_dict

def res_group(input_dict):

    """
    Count residues according to their residue group
    - non polar aliphatic
    - polar uncharged
    - basic positively charged
    - acid negatively charged
    - non polar aromatic
    """
    non_pol_ali = ['Glycine','Alanine','Valine','Leucine','Isoleucine','Methionine']
    pol_un = ['Serine','Threonine','Cysteine','Proline','Asparagine','Glutamine']
    bas_pos = ['Lysine','Arginine','Histidine']
    acid_neg = ['Aspartate','Glutamate']
    non_pol_aro = ['Phenylalanine','Tyrosine','Tryptophan']
    list_of_groups = [non_pol_ali, pol_un, bas_pos, acid_neg, non_pol_aro]
    new_dict = {'Non polar aliphatic':0,'Polar uncharged':0,'Basic positively charged':0,'Acid negatively charged':0,'Non polar aromatic':0}
    for group in list_of_groups:
        for res in group:
            if group == non_pol_ali:
                new_dict['Non polar aliphatic'] = new_dict['Non polar aliphatic'] + float(input_dict[res])
            if group == pol_un:
                new_dict['Polar uncharged'] = new_dict['Polar uncharged'] + float(input_dict[res])
            if group == bas_pos:
                new_dict['Basic positively charged'] = new_dict['Basic positively charged'] + float(input_dict[res])
            if group == acid_neg:
                new_dict['Acid negatively charged'] = new_dict['Acid negatively charged'] + float(input_dict[res])
            if group == non_pol_aro:
                new_dict['Non polar aromatic'] = new_dict['Non polar aromatic'] + float(input_dict[res])
    return new_dict


def analyse_interface(input_file, seq, out_name, delimiter = ";"):

    """
    Run all interface analysis
    """
    interface_res = []
    values_row = open(input_file, "r").readlines()[0].replace("\n","").split(delimiter)
    for row_val, row_name in zip(values_row, seq):
        if row_val != "-":
            interface_res.append(row_name)
        else:
            interface_res.append("-")
    total_res = total_counter(interface_res,'-')
    interface_res_type = res_type_counter(interface_res)
    interface_group_type = res_group(interface_res_type)
    return interface_res, total_res, interface_res_type, interface_group_type

def apply_on_folder(input_alignments, input_dimers):

    """
    Run interface analysis for all the .pdb files in the folder
    """
    DR_interfaces = {}
    protein_interfaces = {}
    for receptor in input_dimers.keys():
        for partner in input_dimers[receptor]:
            cocomaps_name_DR = results_folder  + "/" + receptor + "_" + partner + "_toDR_cocomaps.csv"
            output_name_DR = results_folder  + "/" + receptor + "_" + partner + "_toDR_amino.csv"
            cocomaps_name_partner = results_folder  + "/" + receptor + "_" + partner + "_toprotein_cocomaps.csv"
            output_name_partner = results_folder  + "/" + receptor + "_" + partner + "_toprotein_amino.csv"
            try:
                interface_DR = analyse_interface(cocomaps_name_DR, input_alignments[receptor], output_name_DR)
                interface_prot = analyse_interface(cocomaps_name_partner, input_alignments[partner], output_name_partner)
            except:
                print("Failed to perform action on:", receptor, partner)

"""
Variable initialization
"""
home = DEFAULT_FOLDER
results_folder = home + "/" + PROCESSED_RESULTS_FOLDER
aligned_proteins = raw_parser.retrieve_clean(ALIGN_VECTORS)
dimers = raw_parser.retrieve_clean(DIMERS)
apply_on_folder(aligned_proteins, dimers)
