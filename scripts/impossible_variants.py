# Shell Plus Model Imports
# diseasevariantsinbenignquery
import csv

from django.conf import settings
from django.contrib.admin.models import LogEntry
from django.contrib.auth import get_user_model
from django.contrib.auth.models import Group
from django.contrib.auth.models import Permission
from django.contrib.auth.models import User
from django.contrib.contenttypes.models import ContentType
from django.contrib.sessions.models import Session
from django.core.cache import cache
from django.db import transaction
from django.db.models import Avg
from django.db.models import Case
from django.db.models import Count
from django.db.models import Exists
from django.db.models import F
from django.db.models import Max
from django.db.models import Min
from django.db.models import OuterRef
from django.db.models import Prefetch
from django.db.models import Q
from django.db.models import Subquery
from django.db.models import Sum
from django.db.models import When
from django.urls import reverse
from django.utils import timezone

from scripts.populate_general_functions import *
from tmh_db.models import Database_Metadata
from tmh_db.models import Flank
from tmh_db.models import Flank_residue
from tmh_db.models import Funfam
from tmh_db.models import Funfam_residue
from tmh_db.models import Funfamstatus
from tmh_db.models import Go
from tmh_db.models import Keyword
from tmh_db.models import Non_tmh_helix
from tmh_db.models import Non_tmh_helix_residue
from tmh_db.models import Protein
from tmh_db.models import Residue
from tmh_db.models import Signal_peptide
from tmh_db.models import Signal_residue
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
# Shell Plus Django Imports


Variant.objects.all().prefetch_related("residue", "residue__protein")

def impossible_variants(variant_object):
    '''
    Checks if a variant substitution is allowed by 1 change in the codon tail.
    Accepts a variant object.
    Returns a True: The variant is impossible
    Returns a False: The variant is allowed by the codon table.
    '''

    #print(variant_object
    impossible_dict=impossible_subs()
    if str(variant_object.aa_mut) in impossible_dict[str(variant_object.aa_wt)]:
        #print(variant_object.residue.protein.uniprot_id, variant_object.residue.sequence_position, variant_object.aa_wt, variant_object.aa_mut, "is impossible")
        return True
    else:
        return False

def varmap_columns_and_keys(column_headers):
    #column_headers = column_headers.split('\t')
    varmap_col_dictionary = {}
    for column_number, column_title in enumerate(column_headers):
        varmap_col_dictionary[column_title] = column_number
    # print(varmap_col_dictionary)
    return varmap_col_dictionary

def varmap_tsv_to_list(path_to_varmap_tsv):
    list_of_lines=[]
    with open(path_to_varmap_tsv, 'r', encoding='cp437') as f:
        list=csv.reader(f,delimiter="\t")

        for n, i in enumerate(list):
            #print(i)

            if n==0:
                print(i)
                varmap_headers=varmap_columns_and_keys(i)
            list_of_lines.append(i)

    return((varmap_headers,list_of_lines))



variant_dict = {
    "disease": Variant.objects.filter(disease_status="d", variant_source="ClinVar").prefetch_related("residue", "residue__protein"),
    #"gnomad2": Variant.objects.filter(residue__protein__uniprot_id="Q9HCH0", variant_source="gnomAD2").prefetch_related("residue", "residue__protein"),
    #"gnomad3": Variant.objects.filter(residue__protein__uniprot_id="Q9HCH0", variant_source="gnomAD3").prefetch_related("residue", "residue__protein"),
}



path='scripts/external_datasets/clinvar_varmap2019.tsv'

clinvar_list=varmap_tsv_to_list(path)

impossible_variant_list=[]
for dataset in variant_dict:
    for variant in variant_dict[dataset]:
        #print(variant)
        if impossible_variants(variant) == True:
            impossible_variant_list.append([variant.residue.protein.uniprot_id, variant.residue.sequence_position, variant.aa_wt, variant.aa_mut])

for x in impossible_variant_list:
    #print(x)
    for i in clinvar_list[1]:
        uniprot_accession=i[clinvar_list[0]["UNIPROT_ACCESSION"]]
        uniprot_position=i[clinvar_list[0]["SEQ_NO"]]
        try:
            uniprot_position = int(uniprot_position)
        except(ValueError):
            pass
        reference_aa=i[clinvar_list[0]["UNIPROT_AA"]]
        variant_aa=i[clinvar_list[0]["AA_CHANGE"]]
        if len(variant_aa) == 3 and "/" in variant_aa:
            #if variant_aa.split("/")[0] != reference_aa:
            #    print("OH GOD THIS BROKE!")
            #    print([uniprot_accession, uniprot_position, reference_aa, variant_aa])
            variant_aa = variant_aa.split("/")[1]
            #print([uniprot_accession, uniprot_position, reference_aa, variant_aa], x)

        #print(variant_aa)
        if [uniprot_accession, uniprot_position, reference_aa, variant_aa] == x:
            print(i)
