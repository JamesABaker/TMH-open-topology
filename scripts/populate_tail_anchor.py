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
from tmh_db.models import Database_Metadata, Funfam, Funfam_residue, Funfamstatus, Go, Keyword, Protein, Residue, Structural_residue, Structure, Subcellular_location, Tmh, Tmh_deltag, Tmh_hydrophobicity, Tmh_residue, Tmh_tmsoc, Uniref, Variant

# env LDFLAGS="-I/usr/local/opt/openssl/include -L/usr/local/opt/openssl/lib" pip install psycopg2
from django.db import models
from datetime import datetime, timedelta
from django.utils import timezone
from datetime import date
import pytz
from scripts.populate_general_functions import *


def get_ta():
    with open("scripts/external_datasets/ta_reference_list.csv") as f:
        ta_protein_list = f.readlines()
    clean_ta_protein_list=[]
    for i in ta_protein_list:
        clean_ta_protein_list.append(clean_query(i))
    #print(clean_ta_protein_list)
    return(clean_ta_protein_list)


def cross_ref_ta_db():
    ta_list=get_ta()
    for i in ta_list:
        protein=Protein.objects.filter(uniprot_id=i)
        if protein.count()==1:
            protein.update(tail_anchor=True, tail_anchor_evidence="Thesis")

def run():
    cross_ref_ta_db()
