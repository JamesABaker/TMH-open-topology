# Shell Plus Model Imports
from django.contrib.admin.models import LogEntry
from django.contrib.contenttypes.models import ContentType
from django.contrib.sessions.models import Session
from tmh_db.models import Binding_residue, Database_Metadata, Funfam_residue, Funfamstatus, Go, Keyword, Protein, Residue, Structural_residue, Structure, Tmh, Tmh_deltag, Tmh_hydrophobicity, Tmh_residue, Tmh_tmsoc, Variant
# Shell Plus Django Imports
from django.core.cache import cache
from django.conf import settings
from django.contrib.auth import get_user_model
from django.db import transaction
from django.db.models import Avg, Case, Count, F, Max, Min, Prefetch, Q, Sum, When, Exists, OuterRef, Subquery
from django.utils import timezone
from django.urls import reverse
from scripts.populate_general_functions import *
from scripts.graphs import *


#This will eventually act as a test file


def protein_query(args):
    print(args)
    return()


def raw_variant_query(*args):
    f = [i.strip('\n').split(',') for i in open(*args)]
    tmh_count = 0
    tmh_flank_count = 0
    non_tmh_count = 0

    total_flank = 0
    total_tmh = 0
    total_non_tmh= 0

    #tmh_normalised_per_protein = 0
    tmh_flank_count_normalised_per_protein = 0
    non_tmh_count_normalised_per_protein = 0
    tmp_ids_done=[]

    for paired_data in f:

        uniprot_ids = clean_query(paired_data[0])
        positions = clean_query(paired_data[1])

        total_residues = Residue.objects.filter(protein__uniprot_id=uniprot_ids).count()
        non_tmh_residues = Residue.objects.filter(protein__uniprot_id=uniprot_ids, tmh_residue=None).count()
        tmh_residues = Residue.objects.filter(protein__uniprot_id=uniprot_ids, tmh_residue__evidence="UniProt", tmh_residue__feature_location="TMH").count()
        flank_residues = Residue.objects.filter(protein__uniprot_id=uniprot_ids, tmh_residue__evidence="UniProt").exclude(tmh_residue__feature_location="TMH").count()

        if uniprot_ids not in tmp_ids_done:
            print("\n",uniprot_ids,"\nTotal residues\t",total_residues ,"\nTMH residues\t", tmh_residues, "\nNon-TMH residues\t", non_tmh_residues, "\nFlank Residues\t", flank_residues)
            total_flank = total_flank + flank_residues
            total_tmh = total_tmh + tmh_residues
            total_non_tmh= total_non_tmh + non_tmh_residues
            tmp_ids_done.append(uniprot_ids)


        tmh_residue_type="non-TMH"
        is_it_in_a_tmh = Residue.objects.filter(protein__uniprot_id=uniprot_ids, sequence_position=positions, tmh_residue__evidence="UniProt")
        try:
            if is_it_in_a_tmh.count() == 0:
                non_tmh_count=non_tmh_count+1
                #non_tmh_count_normalised_per_protein=non_tmh_count_normalised_per_protein+1/(non_tmh_residues/total_residues)


            elif is_it_in_a_tmh.count() == 1:
                is_it_in_a_tmh_core = Residue.objects.filter(protein__uniprot_id=uniprot_ids, sequence_position=positions, tmh_residue__evidence="UniProt", tmh_residue__feature_location="TMH")
                if is_it_in_a_tmh_core.count() == 1:
                    tmh_residue_type="TMH"
                    tmh_count=tmh_count+1
                    #tmh_normalised_per_protein = tmh_normalised_per_protein+1/(tmh_residues/total_residues)

                is_it_in_a_tmh_flank = Residue.objects.filter(protein__uniprot_id=uniprot_ids, sequence_position=positions, tmh_residue__evidence="UniProt").exclude(tmh_residue__feature_location="TMH")
                if is_it_in_a_tmh_flank.count() == 1:
                    tmh_residue_type="Flank"
                    tmh_flank_count=tmh_flank_count+1
                    #tmh_flank_count_normalised_per_protein=tmh_flank_count_normalised_per_protein+1/(flank_residues/total_residues)

        except(ZeroDivisionError):
            pass

        print(uniprot_ids,"\t",positions, "\t", tmh_residue_type)


    objects = ["TMH", "±5 flanking residues", "Not in TMH nor ±5 residues"]
    performance = [tmh_count, tmh_flank_count, non_tmh_count]
    source="Cys loop receptor disease variants"
    state = "d"
    print(source, "\n" ,objects, "\n", performance)
    barchart(objects, performance, source, state, "Variant location", "Variant count")

    #objects = ["TMH", "±5 flanking residues", "Not in TMH nor ±5 residues"]
    #performance = [tmh_normalised_per_protein, tmh_flank_count_normalised_per_protein, non_tmh_count_normalised_per_protein]
    #source="Cys loop receptor disease variants normalised 1 div per protein tm to non-tm ratio"
    #state = "d"
    #print(source, "\n" ,objects, "\n", performance)
    #barchart(objects, performance, source, state, "Variant location", "Variant count")

    objects = ["TMH", "±5 flanking residues", "Not in TMH nor ±5 residues"]
    performance = [tmh_count/total_tmh, tmh_flank_count/total_flank, non_tmh_count/total_non_tmh]
    source="Cys loop receptor disease variants normalised at global tm to non-tm ratio"
    state = "d"
    print(source, "\n" ,objects, "\n", performance)
    barchart(objects, performance, source, state, "Variant location", "Variant count")


    return()


def run(*args):
    'python manage.py runscript query --script-args filename'
    #if len(args) == 1:
    #    protein_query(*args)

    raw_variant_query(*args)
