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

a=Structure.objects.all().distinct('pk')

print("pdb, gnomAD, disease, pdb_residues, uniprot_residues")
for i in a.values_list("pdb_id"):
    id=clean_query(str(i))
    gnomad=Variant.objects.filter(residue__structural_residue__structure__pdb_id=id, variant_source__contains="gnomAD").distinct('pk').count()
    disease=Variant.objects.filter(residue__structural_residue__structure__pdb_id=id, disease_status="d").distinct('pk').count()
    pdb_residues=Structural_residue.objects.filter(structure__pdb_id=id).distinct('pk').count()
    uniprot_residues=Residue.objects.filter(structural_residue__structure__pdb_id=id).distinct('pk').count()
    print(id, gnomad, disease, pdb_residues, uniprot_residues)
