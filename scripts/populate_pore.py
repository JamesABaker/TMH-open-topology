from __future__ import division

import json
import os
import time
import urllib
from subprocess import check_output

import defusedxml.ElementTree as ET
import numpy as np
import pytz
import requests
from Bio import SeqIO
from Bio import SwissProt
from django.db import models
from django.utils import timezone
from requests import get

from scripts.populate_general_functions import *
# env LDFLAGS="-I/usr/local/opt/openssl/include -L/usr/local/opt/openssl/lib" pip install psycopg2


def check_porewalker():
    '''
    This checks the channel containing proteins against any already computed pore walker results and submits left overs to porewalker for next time.
    '''

    structures = fetch_channel_structures()
    for pdb_code in structures:

        pdb_code=clean_query(str(pdb_code))

        #porewalker_url = f"https://www.ebi.ac.uk/thornton-srv/software/PoreWalker/Results/vartmh{pdb_code}/vartmh{pdb_code}-marked-pdb.pdb"
        #porewalker_file = f"scripts/external_datasets/porewalker_results/{pdb_code}.pdb"
        ##try:
        #porewalker_url = f"https://www.ebi.ac.uk/thornton-srv/software/PoreWalker/Results/vartmh{pdb_code}/vartmh{pdb_code}-marked-pdb.pdb"
        #print(porewalker_url)
        #download(porewalker_url, porewalker_file)
        #porewalker_pdb_residue_list = open_porewalker_pdb(porewalker_file)
        #porewalker_to_database(pdb_code, porewalker_pdb_residue_list)
        #except:
        #    print(pdb_code, "Needs to be added to the laterbase.")



        porewalker_url = f"https://www.ebi.ac.uk/thornton-srv/software/PoreWalker/Results/vartmh{pdb_code}/vartmh{pdb_code}-aa_list.txt"
        porewalker_file = f"scripts/external_datasets/porewalker_results/{pdb_code}.txt"

        download(porewalker_url, porewalker_file)
        #porewalker_pdb_residue_list = open_porewalker_file(porewalker_file)
        porewalker_residue_list = open_porewalker_file(porewalker_file)
        porewalker_to_database(pdb_code, porewalker_residue_list)



def porewalker_to_database(this_pdb_id, residues):
    '''
    Goes through each line and adds it to the django databases
    '''
    #target_structure = Structure.objects.get(pdb_id=pdb_id)
    #completed_residues=[]
    for residue in residues:
        # pdb = pdb_residue_parse(residue)
        pdb = txt_residue_parse(residue) # Each residue is a line
        if pdb is not False:
            #if pdb["author_position"] not in completed_residues:
                #Structural_residue.objects.filter(structure__pdb_id=this_pdb_id, author_position=pdb["author_position"], pdb_chain=pdb["chain"]).update(porewalker_score=float(pdb['b_factor']))
            Structural_residue.objects.filter(structure__pdb_id=this_pdb_id, author_position=pdb["pdb_position"], pdb_chain=pdb["chain"]).update(pore_residue=True, porewalker_score=float(pdb['b_factor'])
            print(this_pdb_id, "residues containing pore information added to the database.")
            #else:
            #    pass

def pdb_residue_parse(pdb_line):
    pdb_line_list = pdb_line.split()
    print(pdb_line_list)
    if len(pdb_line_list) == 11:
        pdb_line_dictionary = {
            "atom_type": pdb_line_list[0],
            "atom_number": pdb_line_list[1],
            "element_type": pdb_line_list[2],
            "residue_type": pdb_line_list[3],
            "chain": pdb_line_list[4],
            "pdb_position": pdb_line_list[5],
            "x": pdb_line_list[6],
            "y": pdb_line_list[7],
            "z": pdb_line_list[8],
            "occupancy": pdb_line_list[9],
            "b_factor": pdb_line_list[10]
        }
        return(pdb_line_dictionary)
    else:
        return(False)



def txt_residue_parse(txt_line):
    txt_line_list = txt_line.split()
    print(txt_line_list)
    if len(txt_line_list) == 4:
        txt_line_dict = {
            "residue_type": txt_line_list[0],
            "pdb_position": txt_line_list[1],
            "chain": txt_line_list[2],
            "X-coordinate": txt_line_list[3],
        }
        return(txt_line_dict)

    elif txt_line_list[0] == "AA":
        return(False)

    else:
        return(False)



def open_porewalker_file(file_location):
    '''
    Opens the file and returns a list of lines.
    '''
    with open(file_location) as f:
        lineList = f. readlines()
    return(lineList)


def fetch_channel_structures():
    '''
    This uses Uniprot Keywords to fetch as many pore-containing proteins as possible and returns a list of corresponding pdb ids.
    '''
    list_of_channel_keyword= ["Ion channel", "Ion transport", "Translocase", "Translocation", "Sugar transport", "Nuclear pore complex", "Amino-acid transport", "Ammonia transport", "Copper transport", "Electron transport", "ER-Golgi transport", "Iron transport", "Lipid transport", "mRNA transport", "Neurotransmitter transport", "Peptide transport", "Phosphate transport", "Protein transport", "Sugar transport", "Transport", "Zinc transport"]
    pdb_ids=Structure.objects.filter(uniprot_protein_id__keywords__keyword__in=list_of_channel_keyword, uniprot_protein_id__total_tmh_number__gte=1).distinct('pk').values_list("pdb_id")
    # Just incase of duplicates we list a set of the list.
    return(list(set(pdb_ids)))


def run():
    check_porewalker()
    #check_porelogo()
