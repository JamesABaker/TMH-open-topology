# Shell Plus Model Imports
from django.contrib.admin.models import LogEntry
from django.contrib.auth.models import Group, Permission, User
from django.contrib.contenttypes.models import ContentType
from tmh_db.models import Database_Metadata, Funfam, Funfam_residue, Funfamstatus, Go, Keyword, Protein, Residue, Structural_residue, Structure, Subcellular_location, Tmh, Tmh_deltag, Tmh_hydrophobicity, Tmh_residue, Tmh_tmsoc, Uniref, Variant
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
from scripts.populate_general_functions import impossible_subs, aa_baezo_order, heatmap_array
from matplotlib.colors import LogNorm


aa_list_baezo_order=aa_baezo_order()
impossible_subs_dict=impossible_subs()

def heatmap_normalised_by_heatmap(title, heatmap_one, heatmap_two):
    new_heatmap = []

    for row_number, row in enumerate(heatmap_one):
        new_heatmap.append([])
        for column_number, column in enumerate(heatmap_one):
            new_heatmap[row_number].append('')
            #print(aa_list_baezo_order[row_number], impossible_subs_dict[aa_list_baezo_order[column_number]])
            if heatmap_two[row_number][column_number] == 0 or heatmap_one[row_number][column_number] == 0:
                value=0
            elif aa_list_baezo_order[row_number] in impossible_subs_dict[aa_list_baezo_order[column_number]]:
                value=0
            else:
                value=int(heatmap_one[row_number][column_number])/int(heatmap_two[row_number][column_number])
            new_heatmap[row_number][column_number]=value
    heatmap(np.array(new_heatmap), title, aa_list_baezo_order, "PuRd", None)

def remove_duplicate_variants(list_of_variants):
    remove_duplicate_list_of_variants = set(list_of_variants)
    remove_duplicate_list_of_variants = list(remove_duplicate_list_of_variants)
    truncate_list=[]
    for variant in remove_duplicate_list_of_variants:
        truncate_list.append((variant[0], variant[1]))
    return(truncate_list)

outside_disease_query=Variant.objects.filter(residue__flank_residue__feature_location="Outside flank", residue__flank_residue__flank__tmh__meta_tmh=True, disease_status='d', variant_source="ClinVar").values_list("aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id")
print(len(outside_disease_query), "disease variants in the outside flank")
outside_flank_disease_variants = heatmap_array(remove_duplicate_variants(list(outside_disease_query)), aa_list_baezo_order)

outside_gnomad_query=Variant.objects.filter(residue__flank_residue__feature_location="Outside flank", residue__flank_residue__flank__tmh__meta_tmh=True,  variant_source="gnomAD3").exclude(aa_mut=F("aa_wt")).values_list("aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id")
print(len(outside_gnomad_query), "gnomad variants in the outside flank")
outside_flank_gnomad_variants = heatmap_array(remove_duplicate_variants(list(outside_gnomad_query)), aa_list_baezo_order)



inside_disease_query=Variant.objects.filter(residue__flank_residue__feature_location="Inside flank", residue__flank_residue__flank__tmh__meta_tmh=True, disease_status='d', variant_source="ClinVar").values_list("aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id")
print(len(inside_disease_query), "disease variants in the inside flank")
inside_flank_disease_variants = heatmap_array(remove_duplicate_variants(list(inside_disease_query)), aa_list_baezo_order)

inside_gnomad_query=Variant.objects.filter(residue__flank_residue__feature_location="Inside flank", residue__flank_residue__flank__tmh__meta_tmh=True, variant_source="gnomAD3").exclude(aa_mut=F("aa_wt")).values_list("aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id")
print(len(inside_gnomad_query), "gnomad variants in the inside flank")
inside_flank_gnomad_variants = heatmap_array(remove_duplicate_variants(list(inside_gnomad_query)), aa_list_baezo_order)



tmh_disease_query=Variant.objects.filter(residue__tmh_residue__feature_location="TMH", residue__tmh_residue__tmh_id__meta_tmh=True, disease_status='d', variant_source="ClinVar").values_list("aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id")
print(len(tmh_disease_query), "disease variants in the tmh")
tmh_disease_variants = heatmap_array(remove_duplicate_variants(list(tmh_disease_query)), aa_list_baezo_order)
heatmap(np.array(tmh_disease_variants), "test", aa_list_baezo_order, "PuRd", None)

tmh_gnomad_query=Variant.objects.filter(residue__tmh_residue__feature_location="TMH", residue__tmh_residue__tmh_id__meta_tmh=True, variant_source="gnomAD3").exclude(aa_mut=F("aa_wt")).values_list("aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id")
print(len(tmh_gnomad_query), "gnomad variants in the tmh")
tmh_gnomad_variants = heatmap_array(remove_duplicate_variants(list(tmh_gnomad_query)), aa_list_baezo_order)

# I came across some TMHs labelled incorrectly as helices. They mayhave snuck into the database, so to ensure they are not counted as variants, meta-tmh exclude is needed.
helix_disease_query=Variant.objects.filter(residue__non_tmh_helix_residue__nont_tmh_helix_id__helix_start__gte=0, disease_status='d', variant_source="ClinVar").exclude(residue__tmh_residue__tmh_id__meta_tmh=True).values_list("aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id")
print(len(helix_disease_query), "disease variants in the helix")
helix_disease_variants = heatmap_array(remove_duplicate_variants(list(helix_disease_query)), aa_list_baezo_order)

helix_gnomad_query=Variant.objects.filter(residue__non_tmh_helix_residue__nont_tmh_helix_id__helix_start__gte=0, variant_source="gnomAD3").exclude(aa_mut=F("aa_wt")).exclude(residue__tmh_residue__tmh_id__meta_tmh=True).values_list("aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id")
print(len(helix_gnomad_query), "gnomad variants in the helix")
helix_gnomad_variants = heatmap_array(remove_duplicate_variants(list(helix_gnomad_query)), aa_list_baezo_order)



sp_disease_query=Variant.objects.filter(residue__signal_residue__the_signal_peptide__signal_start__gte=0, disease_status='d', variant_source="ClinVar").values_list("aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id")
print(len(sp_disease_query), "disease variants in the signal peptides")
sp_disease_variants = heatmap_array(remove_duplicate_variants(list(sp_disease_query)), aa_list_baezo_order)

sp_gnomad_query=Variant.objects.filter(residue__signal_residue__the_signal_peptide__signal_start__gte=0, variant_source="gnomAD3").exclude(aa_mut=F("aa_wt")).values_list("aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id")
print(len(sp_gnomad_query), "gnomad variants in the signal peptides")
sp_gnomad_variants = heatmap_array(remove_duplicate_variants(list(sp_gnomad_query)), aa_list_baezo_order)




heatmap_normalised_by_heatmap("ClinVar normalised by gnomad v3 meta-tmh Outside flank", outside_flank_disease_variants, outside_flank_gnomad_variants)
heatmap_normalised_by_heatmap("ClinVar normalised by gnomad v3 meta-tmh Inside flank", inside_flank_disease_variants, inside_flank_gnomad_variants)
heatmap_normalised_by_heatmap("ClinVar normalised by gnomad v3 meta-tmh TMH", tmh_disease_variants, tmh_gnomad_variants)
heatmap_normalised_by_heatmap("ClinVar normalised by gnomad v3 non-TMH Helix", helix_disease_variants, helix_gnomad_variants)
heatmap_normalised_by_heatmap("ClinVar normalised by gnomad v3 Signal Peptides", sp_disease_variants, sp_gnomad_variants)
