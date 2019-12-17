from __future__ import division
import requests
import urllib
from requests import get
import shutil
import numpy as np
import os
import time
import gzip
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
from django.db import models
from datetime import datetime, timedelta
from django.utils import timezone
from datetime import date
import pytz
from scripts.populate_general_functions import *

def uniprot_table(query_id):
    filename = str(f"scripts/external_datasets/uniprot_bin/{query_id}.txt")
    input_format = "swiss"
    feature_type = "SIGNAL"

    target_protein = Protein.objects.get(uniprot_id=query_id)

    #print("Checking UniProt for TM annotation in", query_id, ".")
    for record in SeqIO.parse(filename, input_format):
        # features locations is a bit annoying as the start location needs +1 to match the sequence IO, but end is the correct sequence value.
        for i, f in enumerate(record.features):
            if f.type == feature_type:
                if "UnknownPosition" in str(f.location.start) or "UnknownPosition" in str(f.location.end):
                    #print(record.id, "Unknown position for TMH in record")
                    # This doesn't mean that there are coordinates, just that there is a TM somewhere in the protein.
                    pass
                else:
                    #+1 for human number
                    start_position=int(f.location.start) + 1
                    end_position=int(f.location.end)
                    signal_sequence=str(record.seq[(f.location.start):(f.location.end)])

                    record_for_database, created = Signal_peptide.objects.update_or_create(
                        protein=target_protein,
                        defaults={
                            "signal_sequence": signal_sequence,
                            "signal_start": start_position,
                            "signal_stop": end_position,
                        }
                    )

                    this_signal_peptide = Signal_peptide.objects.get(protein=target_protein, signal_start=start_position, signal_stop=end_position)
                    for signal_residue_number, a_residue in enumerate(signal_sequence):
                        sequence_position = int(start_position) + signal_residue_number
                        specific_residue= Residue.objects.get(protein=target_protein, sequence_position=int(sequence_position))
                        record_for_database, created = Signal_residue.objects.update_or_create(
                            residue=specific_residue,
                            the_signal_peptide= this_signal_peptide,
                            amino_acid_type=a_residue,
                        )

    record_for_database, created = Protein.objects.update_or_create(
        uniprot_id=query_id)

    #print("TM annotation found in", query_id, ".")

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

    # Also, parse the variant files which can be massive.
    # humsavar table
    #print(input_query)
    #print("Starting TMH database population script...")
    os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'tmh_database.settings')

    ### Download uniprot files ###
    inputs = input_query_process(input_query)
    input_queries = inputs[0]
    input_query_set = inputs[1]

    ### If UniProt says it is a TMH, add it to the protein table ###

    for query_number, a_query in enumerate(input_query):
        #print("Checking UniProt bin for", a_query)
        a_query = clean_query(a_query)
        uniprot_bin(a_query)
        #print("Adding UniProt record", a_query, " to table,", query_number + 1, "of", len(input_query), "records...")
        uniprot_table(a_query)
