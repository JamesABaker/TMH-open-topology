import os
# env LDFLAGS="-I/usr/local/opt/openssl/include -L/usr/local/opt/openssl/lib" pip install psycopg2
from django.conf import settings
from django.db import models
from tmh_db.models import Database_Metadata, Uniref, Go, Structure, Structural_residue, Funfam_residue, Funfamstatus, Protein, Residue, Tmh, Tmh_deltag, Tmh_hydrophobicity, Tmh_residue, Tmh_tmsoc, Variant, Keyword, Binding_residue
from datetime import datetime, timedelta
from django.utils import timezone
from scripts.populate_general_functions import *



for i in list(Structure.objects.filter(uniprot_protein__keywords__keyword="Ion channel").values_list("pdb_id")):
    print(clean_query(str(i)))


Residue.objects.filter(structural_residue__pdb_position=1121, protein__structure__pdb_id="6j8h").values_list("protein__uniprot_id")

Residue.objects.filter(structural_residue__pdb_position=1121, protein__structure__pdb_id="6j8h").values_list("tmh_residue__tmh_id__tmh_id")

Tmh.objects.filter(protein__uniprot_id="Q15858", tmh_number=24)
