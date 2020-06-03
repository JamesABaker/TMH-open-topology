from __future__ import division
import json
import os
import re
import shutil
import subprocess
import time
import urllib
from datetime import date
from datetime import datetime
from datetime import timedelta
import numpy as np
import pytz
import requests
from Bio import AlignIO
from django.conf import settings
from django.db import models
from django.utils import timezone
from requests import get
from scripts.populate_general_functions import *
# env LDFLAGS="-I/usr/local/opt/openssl/include -L/usr/local/opt/openssl/lib" pip install psycopg2
import time


pause_time=5
#output("Usage:\npython manage.py runscript populate --traceback")

# How many days should be allowed to not enforce updates
today = date.today()
todaysdate = today.strftime("%d_%m_%Y")


def uniprot_to_funfams(a_query):
    '''
    Downloads the funfams associated with a uniprot id.
    It accepts a uniprot id.
    It returns a list of tuples for the superfamily and funfam id.
    '''
    funfam_url = f"http://www.cathdb.info/version/v4_2_0/api/rest/uniprot_to_funfam/{a_query}?content-type=application/json"
    funfam_file = f'scripts/external_datasets/funfam_bin/json/{a_query}.json'
    if check_local_file(funfam_file)==False:
        download(funfam_url, funfam_file, pause=pause_time)
    with open(funfam_file, "r") as json_file:
        funfam_json = json.load(json_file)
    funfams = []
    for i in funfam_json["data"]:
        funfams.append((i["superfamily_id"], i["funfam_number"]))
    return(funfams)


def funfam_to_stockholm(uniprot_id, superfamily_number, funfam_number):
    '''
    Takes the funfams and superfamilies and returns the stockholm alignment.
    '''
    #protein = Protein.objects.get(uniprot_id=a_query)
    stockholm_url = f"http://www.cathdb.info/version/v4_2_0/superfamily/{superfamily_number}/funfam/{funfam_number}/files/stockholm"
    stockholm_file = f'scripts/external_datasets/funfam_bin/stockholm/{superfamily_number}/{funfam_number}.sth'
    cath_superfamily_folder=f'scripts/external_datasets/funfam_bin/stockholm/{superfamily_number}'
    if check_local_file(stockholm_file)==False:
        if not os.path.exists(cath_superfamily_folder):
            os.makedirs(cath_superfamily_folder)
        download(stockholm_url, stockholm_file, pause=pause_time)
    return(stockholm_file)


def stockholm_to_database(a_query, stockholm_file_location):
    '''
    Parses a stockholm alignment and links residues and proteins to the funfam
    id in the VarTMH database.
    It adds the protein to the funfam, and the residue to the funfam_residue.
    '''
    print(stockholm_file_location)
    #align = AlignIO.read(stockholm_file_location, "stockholm")
    #for record in align:
    #    print(str(record.id))
    #    if str(record.id) in a_query:
    #        #print(f"{record.id}, {len(record)}")
    #        print(record.letter_annotations["scorecons_70"])
    #    pass
    return()


def run():
    '''
    This is what django runs.
    This is effectively the canonical script, even though django forces
    it to be in a function.
    '''

    ### Canonical script starts here ###

    input_query = input_query_get()

    os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'tmh_database.settings')

    ### Download uniprot files ###
    inputs = input_query_process(input_query)
    input_queries = inputs[0]
    input_query_set = inputs[1]

    for a_query in input_query:
        funfams_superfamilies = uniprot_to_funfams(str(clean_query(a_query)))
        if len(funfams_superfamilies) > 0:
            for a_funfam in funfams_superfamilies:
                superfamily_id = a_funfam[0]
                funfam_id = a_funfam[1]
                stockholm_file = funfam_to_stockholm(
                    a_query, superfamily_id, funfam_id)
                stockholm_to_database(a_query, stockholm_file)
