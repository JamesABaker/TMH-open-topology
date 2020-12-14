from __future__ import division

import ast
import collections
import gzip
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
from subprocess import check_output

import Bio
import defusedxml.ElementTree as ET
import numpy as np
import pytz
import requests
from Bio import SeqIO
from Bio import SwissProt
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from django.conf import settings
from django.db import models
from django.utils import timezone
from requests import get

from scripts.populate_general_functions import *

# env LDFLAGS="-I/usr/local/opt/openssl/include -L/usr/local/opt/openssl/lib" pip install psycopg2

a = Structure.objects.all().distinct("pk")

print("pdb, gnomAD, disease, pdb_residues, uniprot_residues, tmh_residues")
for i in a.values_list("pdb_id"):
    id = clean_query(str(i))
    gnomad = (
        Variant.objects.filter(
            residue__structural_residue__structure__pdb_id=id,
            variant_source__contains="gnomAD",
            residue__protein__total_tmh_number__gte=1,
        )
        .distinct("pk")
        .count()
    )
    disease = (
        Variant.objects.filter(
            residue__structural_residue__structure__pdb_id=id,
            disease_status="d",
            residue__protein__total_tmh_number__gte=1,
        )
        .distinct("pk")
        .count()
    )
    pdb_residues = (
        Structural_residue.objects.filter(
            structure__pdb_id=id, residue__protein__total_tmh_number__gte=1
        )
        .distinct("pk")
        .count()
    )
    uniprot_residues = (
        Residue.objects.filter(
            structural_residue__structure__pdb_id=id, protein__total_tmh_number__gte=1
        )
        .distinct("pk")
        .count()
    )
    tmh_residues = (
        Structural_residue.objects.filter(
            structure__pdb_id=id, residue__tmh_residue__tmh_id__meta_tmh=True
        )
        .distinct("pk")
        .count()
    )

    print(
        id,
        ",",
        gnomad,
        ",",
        disease,
        ",",
        pdb_residues,
        ",",
        uniprot_residues,
        ",",
        tmh_residues,
    )
