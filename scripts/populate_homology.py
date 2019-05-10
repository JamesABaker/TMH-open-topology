from __future__ import division
import requests
import urllib
from requests import get
import shutil
import numpy as np
import os
import collections
import time
import gzip
import subprocess
import json
from subprocess import check_output
import re
import sys
import defusedxml.ElementTree as ET
import Bio
from Bio import SeqIO
from Bio import SwissProt
from Bio.SeqUtils.ProtParam import ProteinAnalysis
# env LDFLAGS="-I/usr/local/opt/openssl/include -L/usr/local/opt/openssl/lib" pip install psycopg2
from django.conf import settings
from django.db import models
from tmh_db.models import Database_Metadata, Uniref, Go, Structure, Structural_residue, Funfam_residue, Funfamstatus, Protein, Residue, Tmh, Tmh_deltag, Tmh_hydrophobicity, Tmh_residue, Tmh_tmsoc, Variant, Keyword, Binding_residue
from datetime import datetime, timedelta
from django.utils import timezone
from datetime import date
import pytz
from lxml import etree
from scripts.populate_general_functions import *

print("Usage:\npython manage.py runscript populate --traceback")

# How many days should be allowed to not enforce updates
time_threshold = 7
today = date.today()
todaysdate = today.strftime("%d_%m_%Y")


def funfam_submit(a_query):
    print(a_query)
    protein = Protein.objects.get(uniprot_id=a_query)
    record_for_database, created = Funfamstatus.objects.update_or_create(
        protein=protein,
        defaults={
        }
    )

    # check our database for submission key to funfam
    funfam_key = Funfamstatus.objects.get(protein=protein).submission_key
    funfam_completed_date = Funfamstatus.objects.get(
        protein=protein).completed_date
    print("Funfam key for query", a_query, ":", funfam_key)

    # Convert the UniProt binned file to a fasta.
    fasta_file = f"./scripts/external_datasets/fasta_bin/{a_query}.fasta"
    # print(fasta_file)
    with open(str(f"./scripts/external_datasets/uniprot_bin/{a_query}.txt"), "rU") as input_handle:
        with open(str(fasta_file), "w") as output_handle:
            sequences = SeqIO.parse(input_handle, "swiss")
            count = SeqIO.write(sequences, output_handle, "fasta")

        # submit the query to CATH funfam
        base_url = 'http://www.cathdb.info/search/by_funfhmmer'

        with open(fasta_file) as x:
            fasta_contents = x.read()
            data = {'fasta': fasta_contents, "queue": "hmmscan_funvar"}
            headers = {'accept': 'application/json'}
            if len(str(fasta_contents)) > 2000:  # Quick FIX!!!
                pass
            else:

                r = requests.post(base_url, data=data, headers=headers)
                funfam_submission_code = r.json()
                print(funfam_submission_code)
                funfam_key = funfam_submission_code['task_id']

                print("submitted task: " + funfam_key)

            record_for_database, created = Funfamstatus.objects.update_or_create(
                protein=protein,
                defaults={
                    "submission_key": funfam_key,
                }
            )
        return(funfam_key)


def funfam_result(a_query, funfam_submission_code):
    base_url = f'http://www.cathdb.info/search/by_funfhmmer/check/{funfam_submission_code}'

    headers = {'accept': 'application/json'}

    r = requests.get(base_url, headers=headers)
    # Result is something like this: {'success': 0, 'data': {'date_completed': '', 'status': 'queued', 'worker_hostname': '', 'id': '715b00cba220424897cb09df7e81129f', 'date_started': ''}, 'message': 'queued'}
    funfam_status = r.json()
    print(funfam_status)

    # Results can take a while to complete. Best to just add those that have finished. A week in and everything should have settled down.
    if funfam_status['success'] == 1:
        headers = {'accept': 'application/json'}
        results_url = f'http://www.cathdb.info/search/by_funfhmmer/results/{funfam_submission_code}'
        r = requests.get(results_url, headers=headers)
        if str(r) == "<Response [204]>":
            print("No funfam hits for ", a_query)
        elif str(r) == "<Response [200]>":
            funfam_api_result = r.json()
            protein = Protein.objects.get(uniprot_id=a_query)
            funfam_to_update = Funfamstatus.objects.get(protein=protein)
            # record_for_database, created = Funfamstatus.objects.update(
            #    protein=protein,
            #    completed_date=timezone.now(),

            #    #"funfam":

            # )
            funfam_to_update = Funfamstatus.objects.get(protein=protein)
            funfam_to_update.completed_date = timezone.now()  # change field
            funfam_to_update.save()  # this will update only

            # This next bit isn't perfect. We assign each residue in the region the funfam score and id. There may be a way to elegantly put in another table, but given the queries we are asking, this will suffice.
            for key, value in funfam_api_result.items():
                print("Key:", key, "Value:", value)
                print("\n")

def uniref(a_query):
    url = 'https://www.uniprot.org/uploadlists/'

    params = {
        'from': 'ACC',
        'to': 'NF90',
        'format': 'tab',
        'query': a_query
    }

    data = urllib.parse.urlencode(params)
    request = urllib.request.Request(url, data.encode('utf-8'))
    # Please set a contact email address here to help us debug in case of problems (see https://www.uniprot.org/help/privacy).
    contact = ""
    request.add_header('User-Agent', 'Python %s' % contact)
    #request=bytes(request, 'utf-8')
    response = urllib.request.urlopen(request)
    page = response.read(200000)
    page=page.decode(encoding='utf-8', errors='strict')
    page=page.split('\n')
    for n, i in enumerate(page):
        i=i.split('\t')
        page[n]=i

    representative_id=page[1][1]

    print(a_query, "is represented in UniRef90 by", representative_id)
    protein = Protein.objects.get(uniprot_id=a_query)
    representative_id = clean_query(representative_id)
    record_for_database, created = Uniref.objects.update_or_create(
        protein=protein,
        representitive_id=representative_id,
    )

def run():
    '''
    This is what django runs. This is effectively the canonical script,
    even though django forces it to be in a function.
    This will go through several databases and extract the TMH boundaries from proteins,
    and then identify which variants are in those TMHs.
    $ python3 manage.py runscript populate --traceback
    '''

    ### Canonical script starts here ###

    # In full scale mode it will take a long time which may not be suitable for development.
    #input_query = get_uniprot()
    # Here we will just use a watered down list of tricky proteins. Uncomment this line for testing the whole list.
    input_query = ["P01850", "P22760", "Q5K4L6","Q7Z5H4", "P32897", "Q9NR77", "P31644", "Q9NS61"]

    # Also, parse the variant files which can be massive.
    # humsavar table
    print(input_query)
    print("Starting TMH database population script...")
    os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'tmh_database.settings')

    ### Download uniprot files ###
    inputs = input_query_process(input_query)
    input_queries = inputs[0]
    input_query_set = inputs[1]

    for a_query in input_query:
        a_query = clean_query(a_query)
        uniref(a_query)
    # The funfams need to be submitted, then checked for status and results.
    # This submits all the ids to the funfams and gets job ids.
    print("Finding closest funfams for records...")

    uniprotid_funfam_dict = {}
    for a_query in input_query:
        a_query = clean_query(a_query)
        this_funfam = funfam_submit(a_query)
        uniprotid_funfam_dict.update({a_query: this_funfam})

    # This uses the job id to wait until the job is complete and fetch the result.
    for a_query in input_query:
        a_query = clean_query(a_query)
        funfam = funfam_result(a_query, uniprotid_funfam_dict[a_query])
