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
import defusedxml.ElementTree as ET
import Bio
from Bio import SeqIO
from Bio import SwissProt
from Bio.SeqUtils.ProtParam import ProteinAnalysis

# env LDFLAGS="-I/usr/local/opt/openssl/include -L/usr/local/opt/openssl/lib" pip install psycopg2
from django.conf import settings
from django.db import models
from datetime import datetime, timedelta
from django.utils import timezone
from datetime import date
import pytz
import ast
from scripts.populate_general_functions import *

#output("Usage:\npython manage.py runscript populate --traceback")

# How many days should be allowed to not enforce updates
today = date.today()
todaysdate = today.strftime("%d_%m_%Y")


def funfam_submit(a_query):
    # output(a_query)
    protein = Protein.objects.get(uniprot_id=a_query)
    record_for_database, created = Funfamstatus.objects.update_or_create(
        protein=protein,
    )

    # check our database for submission key to funfam
    funfam_key = Funfamstatus.objects.get(protein=protein).submission_key
    #output(a_query, funfam_key)
    if len(funfam_key) < 1 or funfam_key == "NA":
        #output("Funfam key for query", a_query, ":", funfam_key)

        # Convert the UniProt binned file to a fasta.
        fasta_file = f"./scripts/external_datasets/fasta_bin/{a_query}.fasta"
        # #output(fasta_file)
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
                    #output(base_url, data, headers)

                    r = requests.post(base_url, data=data, headers=headers)
                    funfam_submission = r.json()
                    # output(funfam_submission)
                    funfam_key = funfam_submission['task_id']

                    #output("submitted task: " + funfam_key)

                record_for_database, created = Funfamstatus.objects.update_or_create(
                    protein=protein,
                    defaults={
                        "submission_key": funfam_key,
                    }
                )

    else:
        pass
    return(funfam_key)


def funfam_api_to_funfam_hits(result):
    output(result)
    funfam_hits = []
    if len(str(result)) > 2:  # Another lazy way to handle 204
        result = result["funfam_resolved_scan"]["results"][0]["hits"]

        for n, i in enumerate(result):
            # output(i)
            funfam_name = result[n]["match_name"]
            query_start_position = result[n]["hsps"][0]["query_start"]
            # output(funfam_name,query_start_position)
            funfam_hits.append([funfam_name, query_start_position])
    return(funfam_hits)


def query_to_funfam_result_json_file(a_query):
    file = "scripts/external_datasets/funfam_bin/" + a_query + ".ali"
    try:
        text_file = open(file, "r")
        #n = text_file.write('Welcome to pythonexamples.org')
        text_file.close()
    except(FileNotFoundError):
        funfam_submission_key = Funfamstatus.objects.filter(
            protein__uniprot_id=a_query).values_list("submission_key")
        funfam_submission_key = funfam_submission_key[0][0]
        time.sleep(0.1)  # Some sort of limit on CATH
        base_url = f'http://www.cathdb.info/search/by_funfhmmer/check/{funfam_submission_key}'
        headers = {'accept': 'application/json'}
        funfam_status = requests.get(base_url, headers=headers)
        if len(str(funfam_status))>1:
            funfam_status=funfam_status.json()
        # Result is something like this: {'success': 0, 'data': {'date_completed': '', 'status': 'queued', 'worker_hostname': '', 'id': '715b00cba220424897cb09df7e81129f', 'date_started': ''}, 'message': 'queued'}
        # output(funfam_status)
        # Results can take a while to complete. Best to just add those that have finished. A week in and everything should have settled down.
        if "error" in str(funfam_status):
            Funfamstatus.objects.filter(protein__uniprot_id=a_query).delete()
            pass
        elif funfam_status['success'] == 1:
            headers = {'accept': 'application/json'}
            results_url = f'http://www.cathdb.info/search/by_funfhmmer/results/{funfam_submission_key}'

            r = requests.get(results_url, headers=headers)
            output(r)
            output(r.text)

            json_data = r.text

            output(json_data)

            funfam_json_to_file(a_query, json_data)
    return()


def funfam_json_to_file(a_query, result_api):
    file = "scripts/external_datasets/funfam_bin/" + a_query + ".ali"

    funfam_json_file = "scripts/external_datasets/funfam_json_bin/" + a_query + ".json"
    json_file = open(funfam_json_file, "w")
    json_file.write(result_api)
    json_file.close()

    return()


def json_to_database(a_query, json_file):
    with open(json_file, 'r') as file:
        data = file.read().replace('\n', '')
        # output(data)
        # file=str(file).replace("'",'"')
        output(file)
        if len(data) > 1:
            funfam_api_result = json.loads(data)

        else:  # This is a really sloppy way to handle error 204
            funfam_api_result = "{}"
    funfam_match_ids = funfam_api_to_funfam_hits(funfam_api_result)
    output(a_query, "funfam_hits:", funfam_match_ids)

    for hits in funfam_match_ids:
        funfam = hits[0].replace("FF", "funfam")
        url = f'http://www.cathdb.info/version/v4_2_0/superfamily/{funfam}/files/stockholm'
        file = 'scripts/external_datasets/funfam_bin/' + \
            str(funfam.replace("/", "_")) + '.ali'
        if check_local_file(file) == False:
            download(url, file)

        funfam_ali_to_database(a_query, file)


def funfam_ali_to_database(a_query, file):
    protein = Protein.objects.get(uniprot_id=a_query)
    funfam_to_update = Funfamstatus.objects.get(protein=protein)
    # output(funfam_api_result)
    #output("results:", funfam_api_result)
    # record_for_database, created = Funfamstatus.objects.update(
    #    completed_date=timezone.now(),

    #    #"funfam":
    # )
    funfam_to_update = Funfamstatus.objects.get(protein=protein)
    funfam_to_update.completed_date = timezone.now()  # change field
    # funfam_to_update.save()  # this will update only
    with open(file) as json_file:
        for i in json_file:
            if a_query in i:
                # i is the alignment line
                output(i)

        # output(funfam_match_id)
        # This next bit isn't perfect. We assign each residue in the region the funfam score and id. There may be a way to elegantly put in another table, but given the queries we are asking, this will suffice.
        # for key, value in funfam_api_result.items():
        #    if key == "funfam_resolved_scan":
        #        #output("Key:", key, "Value:", value)
        #        for subkey in value.items():

        #            for subsubkey in subkey.items():
        #                output("This is:",subkey)
        #            #if subkey=="match_id":
        #            #    output(subkey)

        #    # output("\n")
        #    pass
        # return([key, value])
    return()


def run():
    '''
    This is what django runs. This is effectively the canonical script,
    even though django forces it to be in a function.
    This will go through several databases and extract the TMH boundaries from proteins,
    and then identify which variants are in those TMHs.
    $ python3 manage.py runscript populate --traceback
    '''

    ### Canonical script starts here ###

    input_query = input_query_get()

    os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'tmh_database.settings')

    ### Download uniprot files ###
    inputs = input_query_process(input_query)
    input_queries = inputs[0]
    input_query_set = inputs[1]

    for a_query in input_query:
        a_query = clean_query(a_query)
        #output("Submitting", a_query, "to FunFam in CATH...")
        funfam_file = "scripts/external_datasets/funfam_json_bin/" + a_query + ".json"
        check_local = check_local_file(funfam_file)

        if check_local == False:
            output(a_query, "is yet to be submitted as a json")
            this_funfam = funfam_submit(a_query)
            #output("id:", a_query, ", key:", this_funfam)
            # This uses the job id to wait until the job is complete and fetch the result.
            #output("Checking results for", a_query, "in FunFam in CATH...")
            funfam = query_to_funfam_result_json_file(a_query)

        # redo this
        check_local = check_local_file(funfam_file)
        if check_local == True:
            output(a_query, "already exists as a json.")
            json_to_database(a_query, funfam_file)
