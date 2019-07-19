
from __future__ import division
import requests
import urllib
from requests import get
import shutil
import os
import collections
import time
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
# env LDFLAGS="-I/usr/local/opt/openssl/include -L/usr/local/opt/openssl/lib" pip install_vars psycopg2
from django.conf import settings
from django.db import models
from tmh_db.models import Database_Metadata, Subcellular_location, Uniref, Go, Structure, Structural_residue, Funfam_residue, Funfamstatus, Protein, Residue, Tmh, Tmh_deltag, Tmh_hydrophobicity, Tmh_residue, Tmh_tmsoc, Variant, Keyword, Binding_residue
from datetime import datetime, timedelta
from datetime import date
import pytz
from scripts.populate_general_functions import *


file = open("lista_uniprots.txt", "r")

all_vars=0
tmh=0
non_tmh=0
flanks=0

print("Total, TMH, Flanks, Non-TMH")
for aline in file:
    aline=clean_query(aline)
    print(aline)
    #print(Variant.objects.filter(residue__protein__uniprot_id=aline).values_list("residue__sequence_position"))
    print(Variant.objects.filter(variant_source="ClinVar", disease_status='d', residue__protein__uniprot_id=aline).values_list( "residue__protein__uniprot_id", "residue__sequence_position"))
    all_vars_id_vars = Variant.objects.filter(variant_source="ClinVar", disease_status='d', residue__protein__uniprot_id=aline).count()
    tmh_id_vars = Variant.objects.filter(variant_source="ClinVar", disease_status='d', residue__protein__uniprot_id=aline, residue__tmh_residue__feature_location="TMH").count()
    flank_id_vars = Variant.objects.filter(variant_source="ClinVar", disease_status='d', residue__protein__uniprot_id=aline, residue__tmh_residue=True).exclude(residue__tmh_residue__feature_location="TMH").count()
    non_tmh_id_vars = Variant.objects.filter(variant_source="ClinVar", disease_status='d', residue__protein__uniprot_id=aline, residue__tmh_residue=None).count()

    print(all_vars_id_vars, tmh_id_vars, flank_id_vars, non_tmh_id_vars)
    all_vars = all_vars + all_vars_id_vars
    tmh = tmh + tmh_id_vars
    flanks = flanks + flank_id_vars
    non_tmh = non_tmh + non_tmh_id_vars

print("Total Total, TMH, Flanks, Non-TMH")
print(all_vars, tmh, flanks, non_tmh)
