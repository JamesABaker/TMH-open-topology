# Shell Plus Model Imports
from django.contrib.admin.models import LogEntry
from django.contrib.auth.models import Group, Permission, User
from django.contrib.contenttypes.models import ContentType
from django.contrib.sessions.models import Session
from tmh_db.models import Database_Metadata, Flank, Flank_residue, Funfam, Funfam_residue, Funfamstatus, Go, Keyword, Non_tmh_helix, Non_tmh_helix_residue, Protein, Residue, Signal_peptide, Signal_residue, Structural_residue, Structure, Subcellular_location, Tmh, Tmh_deltag, Tmh_hydrophobicity, Tmh_residue, Tmh_tmsoc, Uniref, Variant
# Shell Plus Django Imports
from django.core.cache import cache
from django.conf import settings
from django.contrib.auth import get_user_model
from django.db import transaction
from django.db.models import Avg, Case, Count, F, Max, Min, Prefetch, Q, Sum, When, Exists, OuterRef, Subquery
from django.utils import timezone
from django.urls import reverse

def structure_residue_uniprot_residue_check():
    a=Structural_residue.objects.all().distinct('pk').values_list('structure__pdb_id','residue__protein__uniprot_id','structure_aa', 'residue__amino_acid_type')
    test_pass=True
    mismatches=[]
    for i in a:
        if i[2] != i[3]:
            test_pass=False
            mismatches.append(i)

    return(test_pass, mismatches)

def run():
    residue_structure_mismatches=structure_residue_uniprot_residue_check()
    if residue_structure_mismatches[0] == True:
        print("All structures match their uniprot residue aa types.")
    elif residue_structure_mismatches[0] == False:
        print(len(residue_structure_mismatches[1]), "structural residues did not match between their structure and their sequence.")
