# Shell Plus Model Imports
from django.contrib.admin.models import LogEntry
from django.contrib.auth.models import Group, Permission, User
from django.contrib.contenttypes.models import ContentType
from tmh_db.models import Binding_residue, Database_Metadata, Funfam, Funfam_residue, Funfamstatus, Go, Keyword, Pfam, Pfam_residue, Protein, Residue, Structural_residue, Structure, Subcellular_location, Tmh, Tmh_deltag, Tmh_hydrophobicity, Tmh_residue, Tmh_tmsoc, Uniref, Variant
# Shell Plus Django Imports
from django.db import transaction
from django.db.models import Avg, Case, Count, F, Max, Min, Prefetch, Q, Sum, When, Exists, OuterRef, Subquery
from django.utils import timezone
from django.urls import reverse
# Charts
import matplotlib.pyplot as plt
import numpy as np
import collections
from scripts.graphs import *
from scripts.populate_general_functions import *

aa_list_baezo_order=['K', 'R', 'E', 'D', 'Q', 'H', 'N', 'P', 'Y', 'W', 'C', 'M', 'T', 'S', 'G', 'V', 'F', 'A', 'I', 'L']

def heatmap_normalised_by_heatmap(title, heatmap_one, heatmap_two):
    new_heatmap = []
    for row_number, row in enumerate(heatmap_one):
        new_heatmap.append([])
        for column_number, column in enumerate(heatmap_one):
            new_heatmap[row_number].append('')
            try:
                value=heatmap_one[row_number][column_number]/heatmap_two[row_number][column_number]
            except:
                value=0
            new_heatmap[row_number][column_number]=value
    heatmap(np.array(new_heatmap), title, aa_list_baezo_order, "d", None)

def remove_duplicate_variants(list_of_variants):
    remove_duplicate_list_of_variants = set(list_of_variants)
    remove_duplicate_list_of_variants = list(remove_duplicate_list_of_variants)
    truncate_list=[]
    for variant in remove_duplicate_list_of_variants:
        truncate_list.append((variant[0], variant[1]))
    return(truncate_list)

outside_flank_disease_variants = heatmap_array(remove_duplicate_variants(list(Variant.objects.filter(residue__tmh_residue__feature_location="Outside flank", residue__tmh_residue__evidence="TOPDB", disease_status='d', variant_source="ClinVar").values_list("aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id"))), aa_list_baezo_order)
outside_flank_gnomad_variants = heatmap_array(remove_duplicate_variants(list(Variant.objects.filter(residue__tmh_residue__feature_location="Outside flank", residue__tmh_residue__evidence="TOPDB",  variant_source="gnomAD").exclude(aa_mut=F("aa_wt")).values_list("aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id"))), aa_list_baezo_order)

inside_flank_disease_variants = heatmap_array(remove_duplicate_variants(list(Variant.objects.filter(residue__tmh_residue__feature_location="Inside flank", residue__tmh_residue__evidence="TOPDB", disease_status='d', variant_source="ClinVar").values_list("aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id"))), aa_list_baezo_order)
inside_flank_gnomad_variants = heatmap_array(remove_duplicate_variants(list(Variant.objects.filter(residue__tmh_residue__feature_location="Inside flank", residue__tmh_residue__evidence="TOPDB", variant_source="gnomAD").exclude(aa_mut=F("aa_wt")).values_list("aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id"))), aa_list_baezo_order)

tmh_disease_variants = heatmap_array(remove_duplicate_variants(list(Variant.objects.filter(residue__tmh_residue__feature_location="TMH", residue__tmh_residue__evidence="TOPDB", disease_status='d', variant_source="ClinVar").values_list("aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id"))), aa_list_baezo_order)
tmh_gnomad_variants = heatmap_array(remove_duplicate_variants(list(Variant.objects.filter(residue__tmh_residue__feature_location="TMH", residue__tmh_residue__evidence="TOPDB", variant_source="gnomAD").exclude(aa_mut=F("aa_wt")).values_list("aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id"))), aa_list_baezo_order)

heatmap_normalised_by_heatmap("ClinVar normalised by gnomad TOPDB Outside flank", outside_flank_disease_variants, outside_flank_gnomad_variants)
heatmap_normalised_by_heatmap("ClinVar normalised by gnomad TOPDB Inside flank", inside_flank_disease_variants, inside_flank_gnomad_variants)
heatmap_normalised_by_heatmap("ClinVar normalised by gnomad TOPDB TMHnormalised by gnomad", tmh_disease_variants, tmh_gnomad_variants)
