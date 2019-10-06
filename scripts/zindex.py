from __future__ import division
import requests
import urllib
import shutil
import collections
import json
from subprocess import check_output
import re
import defusedxml.ElementTree as ET
import Bio
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
# env LDFLAGS="-I/usr/local/opt/openssl/include -L/usr/local/opt/openssl/lib" pip install_vars psycopg2
from django.db import models
from tmh_db.models import Database_Metadata, Subcellular_location, Uniref, Go, Structure, Structural_residue, Funfam_residue, Funfamstatus, Protein, Residue, Tmh, Tmh_deltag, Tmh_hydrophobicity, Tmh_residue, Tmh_tmsoc, Variant, Keyword, Binding_residue
import pytz
from scripts.populate_general_functions import *
import matplotlib.pyplot as plt


def clean_positions(position_list):
    clean_list=[]
    for entry in position_list:
        for contents in entry:
            if contents == None:
                pass
            else:
                clean_list.append(contents)
    return(clean_list)

def plot_positions(all_positions):
    count_coordinates=[]
    position_coordinates=[]
    for number in range(min(all_positions), max(all_positions)):
        count=0
        for z in all_positions:
            if z == number:
                count = count+1
        count_coordinates.append(count)
        position_coordinates.append(number)
    plt.scatter(position_coordinates, count_coordinates)
    #plt.savefig('test.png')
    plt.show()
    plt.clf()


def get_variant_z_disease():
#print(Variant.objects.filter(residue__protein__uniprot_id=aline).values_list("residue__sequence_position"))
    return(list(Variant.objects.filter(disease_status='d', residue__tmh_residue__evidence="UniProt").values_list("residue__tmh_residue__amino_acid_location_in_to_out")))

def run():
    positions=get_variant_z_disease()
    clean_positions_list=clean_positions(positions)
    plot_positions(clean_positions_list)
