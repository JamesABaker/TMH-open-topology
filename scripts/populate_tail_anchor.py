from __future__ import division

import gzip
import json
import os
import re
import shutil
import sys
import time
import urllib
from datetime import date
from datetime import datetime
from datetime import timedelta
from subprocess import check_output

import Bio
import defusedxml.ElementTree as ET
import numpy as np
import pytz
import requests
from Bio import SeqIO
from Bio import SwissProt
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from django.db import models
from django.utils import timezone
from requests import get

from scripts.populate_general_functions import *
from tmh_db.models import Database_Metadata
from tmh_db.models import Funfam
from tmh_db.models import Funfam_residue
from tmh_db.models import Funfamstatus
from tmh_db.models import Go
from tmh_db.models import Keyword
from tmh_db.models import Protein
from tmh_db.models import Residue
from tmh_db.models import Structural_residue
from tmh_db.models import Structure
from tmh_db.models import Subcellular_location
from tmh_db.models import Tmh
from tmh_db.models import Tmh_deltag
from tmh_db.models import Tmh_hydrophobicity
from tmh_db.models import Tmh_residue
from tmh_db.models import Tmh_tmsoc
from tmh_db.models import Uniref
from tmh_db.models import Variant
# env LDFLAGS="-I/usr/local/opt/openssl/include -L/usr/local/opt/openssl/lib" pip install psycopg2


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
