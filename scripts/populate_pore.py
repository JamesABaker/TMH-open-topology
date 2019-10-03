from __future__ import division
import requests
import urllib
from requests import get
import numpy as np
import os
import time
import json
from subprocess import check_output
import defusedxml.ElementTree as ET
from Bio import SeqIO
from Bio import SwissProt
# env LDFLAGS="-I/usr/local/opt/openssl/include -L/usr/local/opt/openssl/lib" pip install psycopg2
from django.db import models
from tmh_db.models import Database_Metadata, Subcellular_location, Flank, Flank_residue, Uniref, Go, Structure, Structural_residue, Funfam_residue, Funfamstatus, Protein, Residue, Tmh, Tmh_deltag, Tmh_hydrophobicity, Tmh_residue, Tmh_tmsoc, Variant, Keyword, Binding_residue
from datetime import datetime, timedelta
from django.utils import timezone
import pytz
from scripts.populate_general_functions import *


def check_porewalker():
    '''
    This checks the channel containing proteins against any already computed pore walker results and submits left overs to porewalker for next time.
    Spectrum range (0.00000 to 9.00000)
    '''

    structures = fetch_channel_structures()
    for pdb_code in structures:
        pdb_code=clean_query(str(pdb_code))
        porewalker_url = f"https://www.ebi.ac.uk/thornton-srv/software/PoreWalker/Results/vartmh{pdb_code}/vartmh{pdb_code}-marked-pdb.pdb"
        porewalker_file = f"scripts/external_datasets/porewalker_results/{pdb_code}.pdb"
        #try:
        download(porewalker_url, porewalker_file)
        porewalker_pdb_residue_list = open_porewalker_pdb(porewalker_file)
        porewalker_to_database(pdb_code, porewalker_pdb_residue_list)
        #except:
        #    print(pdb_code, "Needs to be added to the laterbase.")


def porewalker_to_database(pdb_id, residues):
    '''
    Goes through each line and adds it to the django databases
    '''
    #target_structure = Structure.objects.get(pdb_id=pdb_id)
    completed_residues=[]
    for residue in residues:
        pdb = pdb_residue_parse(residue)
        if pdb is not False:
            if pdb["pdb_position"] not in completed_residues:
                Structural_residue.objects.filter(pdb_position=pdb["pdb_position"], pdb_chain=pdb["chain"]).update(porewalker_score=pdb['b_factor'])
                print(pdb_id, "residues containing pore information added to the database.")
            else:
                pass


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



def open_porewalker_pdb(file_location):
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
    list_of_channel_keywords = ["Ion channel", "Ion transport", "Translocase", "Translocation", "Sugar transport", "Nuclear pore complex", "Amino-acid transport", "Ammonia transport", "Copper transport", "Electron transport",
                                "ER-Golgi transport", "Iron transport", "Lipid transport", "mRNA transport", "Neurotransmitter transport", "Peptide transport", "Phosphate transport", "Protein transport", "Sugar transport", "Transport", "Zinc transport"]
    pdb_ids = Structure.objects.filter(
        uniprot_protein__keywords__keyword__in=list_of_channel_keywords).values_list("pdb_id")
    # Just incase of duplicates we list a set of the list.
    return(list(set(pdb_ids)))


def run():
    check_porewalker()
    #check_porelogo()
