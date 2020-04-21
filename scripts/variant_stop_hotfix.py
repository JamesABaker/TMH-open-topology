from __future__ import division

import os
import time
from datetime import date
from datetime import datetime
from datetime import timedelta

import pytz
from django.conf import settings
from django.db import models
from django.utils import timezone

from scripts.populate_general_functions import *
# env LDFLAGS="-I/usr/local/opt/openssl/include -L/usr/local/opt/openssl/lib" pip install psycopg2




varmap_files = {
    "ClinVar": "scripts/external_datasets/clinvar_varmap2019.tsv",
    "gnomAD3": "scripts/external_datasets/gnomad_coding_regions3.tsv",
    "gnomAD3": "scripts/external_datasets/gnomad_coding_regions2.tsv"
}



def varmap_columns_and_keys(column_headers):
    column_headers = column_headers.split('\t')
    varmap_col_dictionary = {}
    for column_number, column_title in enumerate(column_headers):
        varmap_col_dictionary[column_title] = column_number
    return(varmap_col_dictionary)

def varmap_header_dict(varmap_file):
    with open(varmap_file, encoding="ISO-8859-1") as inputfile:
        for line_number, var_database_entry in enumerate(inputfile):
            # This is the header line
            if line_number == 0:
                varmap_headers = varmap_columns_and_keys(var_database_entry)
            else:
                break
    return(varmap_headers)

def fix(varmap_file):
    varmap_file_location=varmap_files[varmap_file]
    varmap_index=varmap_header_dict(varmap_file_location)
    myfile=open(varmap_file_location, "r",  encoding="ISO-8859-1")

    while 1:
        lines = myfile.readlines(10000)
        if not lines:
            break
        for line in lines:
            varmap_item=line.rstrip('\n')
            varmap_item = varmap_item.strip().split('\t')
            var_record_location = varmap_item[varmap_index["SEQ_NO"]]
            var_record_id = varmap_item[varmap_index["USER_ID"]]
            uniprot_record = varmap_item[varmap_index["UNIPROT_ACCESSION"]]
            variant_review = "Unknown"
            aa_wt = varmap_item[varmap_index["UNIPROT_AA"]]
            user_id = varmap_item[varmap_index["USER_ID"]]
            if len(varmap_item[varmap_index["AA_CHANGE"]]) == 3 and "/" in varmap_item[varmap_index["AA_CHANGE"]]:
                aa_mut = varmap_item[varmap_index["AA_CHANGE"]].split("/")[1]
                if "*" in aa_mut:
                    print("stop codon (*) found in", user_id, varmap_file )
                    #Variant.objects.filter(variant_source_id=user_id).count()
                    offenses=Variant.objects.filter(variant_source_id=str(user_id), variant_source=str(varmap_file)).count()
                    if offenses > 0:
                        print("Bug identified in", user_id, offenses,"times." )
                        Variant.objects.get(variant_source_id=str(user_id), variant_source=str(varmap_file)).delete()



def run():



    for dataset, location in varmap_files.items():
        print(dataset)
        fix(dataset)
