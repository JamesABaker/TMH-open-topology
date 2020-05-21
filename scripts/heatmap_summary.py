# Shell Plus Model Imports
import collections

import matplotlib.pyplot as plt
import numpy as np
from django.contrib.admin.models import LogEntry
from django.contrib.auth.models import Group
from django.contrib.auth.models import Permission
from django.contrib.auth.models import User
from django.contrib.contenttypes.models import ContentType
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
from matplotlib.colors import LogNorm

# Shell Plus Django Imports

from scripts.graphs import *
from scripts.populate_general_functions import aa_baezo_order
from scripts.populate_general_functions import heatmap_array
from scripts.populate_general_functions import impossible_subs
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

date = datetime.now().strftime('%Y-%m-%d %H:%M:%S')

with open('bug_exclusion_list.txt') as f:
    buggy_uniprots = f.read().splitlines()

aa_list_baezo_order = aa_baezo_order()
impossible_subs_dict = impossible_subs()


def heatmap_normalised_by_heatmap(title, heatmap_one, heatmap_two):
    new_heatmap = []

    for row_number, row in enumerate(heatmap_one):
        new_heatmap.append([])
        for column_number, column in enumerate(heatmap_one):
            new_heatmap[row_number].append('')
            #print(aa_list_baezo_order[row_number], impossible_subs_dict[aa_list_baezo_order[column_number]])
            if heatmap_two[row_number][column_number] == 0 or heatmap_one[row_number][column_number] == 0:
                value = 0
            elif aa_list_baezo_order[row_number] in impossible_subs_dict[aa_list_baezo_order[column_number]]:
                value = 0
            else:
                value = int(heatmap_one[row_number][column_number]) / \
                    int(heatmap_two[row_number][column_number])
            new_heatmap[row_number][column_number] = value
    heatmap(np.array(new_heatmap), title, aa_list_baezo_order,
            "PuRd", None, annotation_format="", bars=False)


def remove_duplicate_variants(list_of_variants):
    remove_duplicate_list_of_variants = set(list_of_variants)
    remove_duplicate_list_of_variants = list(remove_duplicate_list_of_variants)
    truncate_list = []
    for variant in remove_duplicate_list_of_variants:
        if variant[3] not in buggy_uniprots:
            truncate_list.append((variant[0], variant[1]))
    return(truncate_list)


#    redundant_list=[]
#    for variant in list_of_variants:
#        redundant_list.append((variant[0], variant[1]))
#        print(variant)
#    return(redundant_list)
#


print("Feature, disease variants, benign variants, residues")


#Inside flanks
inside_flank_residues=Residue.objects.filter(flank_residue__feature_location="Inside flank", flank_residue__flank__tmh__meta_tmh=True).distinct('pk').count()
inside_disease_query=Variant.objects.filter(residue__flank_residue__feature_location="Inside flank", residue__flank_residue__flank__tmh__meta_tmh=True, disease_status='d', variant_source="ClinVar").distinct('pk').values_list("aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id")
inside_flank_disease_variants = heatmap_array(remove_duplicate_variants(list(inside_disease_query)), aa_list_baezo_order)

inside_benign_query=Variant.objects.filter(residue__flank_residue__feature_location="Inside flank", residue__flank_residue__flank__tmh__meta_tmh=True, disease_status='n', variant_source="ClinVar").distinct('pk').values_list("aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id")
inside_flank_benign_variants = heatmap_array(remove_duplicate_variants(list(inside_benign_query)), aa_list_baezo_order)
print(f"Inside flanks, {len(inside_disease_query)}, {len(inside_benign_query)}, {inside_flank_residues}")


heatmap_normalised_by_heatmap("ClinVar disease normalised by benign meta-tmh inside flank", inside_flank_disease_variants, inside_flank_benign_variants)


#Outside flanks
outside_flank_residues=Residue.objects.filter(flank_residue__feature_location="Outside flank", flank_residue__flank__tmh__meta_tmh=True).distinct('pk').count()
outside_disease_query=Variant.objects.filter(residue__flank_residue__feature_location="Outside flank", residue__flank_residue__flank__tmh__meta_tmh=True, disease_status='d', variant_source="ClinVar").distinct('pk').values_list("aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id")
outside_flank_disease_variants = heatmap_array(remove_duplicate_variants(list(outside_disease_query)), aa_list_baezo_order)

outside_benign_query=Variant.objects.filter(residue__flank_residue__feature_location="Outside flank", residue__flank_residue__flank__tmh__meta_tmh=True, disease_status='n', variant_source="ClinVar").distinct('pk').values_list("aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id")
outside_flank_benign_variants = heatmap_array(remove_duplicate_variants(list(outside_benign_query)), aa_list_baezo_order)

print(f"Outside flanks, {len(outside_disease_query)}, {len(outside_benign_query)}, {outside_flank_residues}")

heatmap_normalised_by_heatmap("ClinVar disease normalised by benign meta-tmh multipass outside flank", outside_flank_disease_variants, outside_flank_benign_variants)


### Single pass

singlepass_residues=Residue.objects.filter(protein__total_tmh_number=1, tmh_residue__tmh_id__meta_tmh=True).distinct('pk').count()

single_tmh_disease_query=Variant.objects.filter(residue__tmh_residue__feature_location="TMH", residue__protein__total_tmh_number=1, residue__tmh_residue__tmh_id__meta_tmh=True, disease_status='d', variant_source="ClinVar").distinct('pk').values_list("aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id")
single_tmh_disease_variants = heatmap_array(remove_duplicate_variants(list(single_tmh_disease_query)), aa_list_baezo_order)

single_tmh_benign_query=Variant.objects.filter(residue__tmh_residue__feature_location="TMH", residue__protein__total_tmh_number=1, residue__tmh_residue__tmh_id__meta_tmh=True, disease_status='n', variant_source="ClinVar").distinct('pk').values_list("aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id")
single_tmh_benign_variants = heatmap_array(remove_duplicate_variants(list(single_tmh_benign_query)), aa_list_baezo_order)

print(f"Singlepass TMHs, {len(single_tmh_disease_query)}, {len(single_tmh_benign_query)}, {singlepass_residues}")

heatmap_normalised_by_heatmap("ClinVar disease normalised by benign meta-tmh singlepass TMHs", single_tmh_disease_variants, single_tmh_benign_variants)


### multipass

multipass_residues=Residue.objects.filter(protein__total_tmh_number__gt=1, tmh_residue__tmh_id__meta_tmh=True).distinct('pk').count()

multi_tmh_disease_query=Variant.objects.filter(residue__tmh_residue__feature_location="TMH", residue__protein__total_tmh_number__gt=1, residue__tmh_residue__tmh_id__meta_tmh=True, disease_status='d', variant_source="ClinVar").distinct('pk').values_list("aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id")
multi_tmh_disease_variants = heatmap_array(remove_duplicate_variants(list(multi_tmh_disease_query)), aa_list_baezo_order)

multi_tmh_benign_query=Variant.objects.filter(residue__tmh_residue__feature_location="TMH", residue__protein__total_tmh_number__gt=1, residue__tmh_residue__tmh_id__meta_tmh=True, disease_status='n', variant_source="ClinVar").distinct('pk').values_list("aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id")
multi_tmh_benign_variants = heatmap_array(remove_duplicate_variants(list(multi_tmh_benign_query)), aa_list_baezo_order)

print(f"Multipass TMHs, {len(multi_tmh_disease_query)}, {len(multi_tmh_benign_query)}, {multipass_residues}")

heatmap_normalised_by_heatmap("ClinVar disease normalised by benign meta-tmh multipass TMHs", multi_tmh_disease_variants, multi_tmh_benign_variants)


### Helix
helix_residues=Residue.objects.filter(non_tmh_helix_residue__nont_tmh_helix_id__helix_start__gte=0).count()

helix_disease_query=Variant.objects.filter(residue__non_tmh_helix_residue__nont_tmh_helix_id__helix_start__gte=0, disease_status='d', variant_source="ClinVar").exclude(residue__tmh_residue__tmh_id__meta_tmh=True).distinct('pk').values_list("aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id")
helix_disease_variants = heatmap_array(remove_duplicate_variants(list(helix_disease_query)), aa_list_baezo_order)

helix_benign_query=Variant.objects.filter(residue__non_tmh_helix_residue__nont_tmh_helix_id__helix_start__gte=0, disease_status='n', variant_source="ClinVar").exclude(residue__tmh_residue__tmh_id__meta_tmh=True).distinct('pk').values_list("aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id")
helix_benign_variants = heatmap_array(remove_duplicate_variants(list(helix_benign_query)), aa_list_baezo_order)

print(f"Non membrane helices, {len(helix_disease_query)}, {len(helix_benign_query)}, {helix_residues}")

heatmap_normalised_by_heatmap("ClinVar disease normalised by benign meta-tmh non-TMHs", helix_disease_variants, helix_benign_variants)

### Pore
pore_residues = Residue.objects.filter(tmh_residue__feature_location="TMH", structural_residue__pore_residue=True, tmh_residue__tmh_id__meta_tmh=True).distinct('pk').count()

pore_disease_query = Variant.objects.filter(residue__tmh_residue__feature_location="TMH", residue__structural_residue__pore_residue=True, residue__tmh_residue__tmh_id__meta_tmh=True, disease_status='d').distinct("pk").values_list("aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id")
pore_disease_variants = heatmap_array(remove_duplicate_variants(list(pore_disease_query)), aa_list_baezo_order)

pore_benign_query = Variant.objects.filter(residue__tmh_residue__feature_location="TMH", residue__structural_residue__pore_residue=True, residue__tmh_residue__tmh_id__meta_tmh=True, disease_status='n').distinct("pk").values_list("aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id")
pore_benign_variants = heatmap_array(remove_duplicate_variants(list(pore_benign_query)), aa_list_baezo_order)

print(f"Pore variants, {len(pore_disease_query)}, {len(pore_benign_query)}, {pore_residues}")

heatmap_normalised_by_heatmap("ClinVar disease normalised by benign meta-tmh pore residues", pore_disease_variants, pore_benign_variants)












def run():
    print("complete")

### Past stuff:

'''
### Multipass starts here ###

#Outside flanks
outside_disease_query=Variant.objects.filter(residue__flank_residue__feature_location="Outside flank", residue__protein__total_tmh_number__gte=2, residue__flank_residue__flank__tmh__meta_tmh=True, disease_status='d', variant_source="ClinVar").values_list("aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id")
outside_disease_query_res=Residue.objects.filter(flank_residue__feature_location="Outside flank", protein__total_tmh_number__gte=2, flank_residue__flank__tmh__meta_tmh=True)
print(len(outside_disease_query), "disease variants in the multipass outside flank of", outside_disease_query_res.count(), "residues.")
outside_flank_disease_variants = heatmap_array(remove_duplicate_variants(list(outside_disease_query)), aa_list_baezo_order)
heatmap(np.array(outside_flank_disease_variants), "ClinVar disease variants in multipass outside flanks", aa_list_baezo_order, "Reds", None)


outside_gnomad3_query=Variant.objects.filter(residue__flank_residue__feature_location="Outside flank", residue__protein__total_tmh_number__gte=2, residue__flank_residue__flank__tmh__meta_tmh=True,  variant_source="gnomAD3").exclude(aa_mut=F("aa_wt")).values_list("aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id")
print(len(outside_gnomad3_query), "gnomad v3 variants in the multipass outside flank")
outside_flank_gnomad3_variants = heatmap_array(remove_duplicate_variants(list(outside_gnomad3_query)), aa_list_baezo_order)
heatmap(np.array(outside_flank_gnomad3_variants), "gnomAD v3 disease variants in multipass outside flanks", aa_list_baezo_order, "Greens", None)


outside_gnomad2_query=Variant.objects.filter(residue__flank_residue__feature_location="Outside flank", residue__protein__total_tmh_number__gte=2, residue__flank_residue__flank__tmh__meta_tmh=True,  variant_source="gnomAD2").exclude(aa_mut=F("aa_wt")).values_list("aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id")
print(len(outside_gnomad2_query), "gnomad v2 variants in the multipass outside flank")
outside_flank_gnomad3_variants = heatmap_array(remove_duplicate_variants(list(outside_gnomad2_query)), aa_list_baezo_order)
heatmap(np.array(outside_flank_gnomad3_variants), "gnomAD v2 disease variants in multipass outside flanks", aa_list_baezo_order, "Greens", None)

#Inside flanks
inside_disease_query=Variant.objects.filter(residue__flank_residue__feature_location="Inside flank", residue__protein__total_tmh_number__gte=2, residue__flank_residue__flank__tmh__meta_tmh=True, disease_status='d', variant_source="ClinVar").values_list("aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id")
print(len(inside_disease_query), "disease variants in the multipass inside flank")
inside_flank_disease_variants = heatmap_array(remove_duplicate_variants(list(inside_disease_query)), aa_list_baezo_order)
heatmap(np.array(inside_flank_disease_variants), "ClinVar disease variants in multipass inside flanks", aa_list_baezo_order, "Reds", None)


inside_gnomad3_query=Variant.objects.filter(residue__flank_residue__feature_location="Inside flank", residue__protein__total_tmh_number__gte=2, residue__flank_residue__flank__tmh__meta_tmh=True, variant_source="gnomAD3").exclude(aa_mut=F("aa_wt")).values_list("aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id")
print(len(inside_gnomad3_query), "gnomad v3 variants in the multipass inside flank")
inside_flank_gnomad3_variants = heatmap_array(remove_duplicate_variants(list(inside_gnomad3_query)), aa_list_baezo_order)
heatmap(np.array(inside_flank_gnomad3_variants), "gnomAD v3 disease variants in multipass inside flanks", aa_list_baezo_order, "Greens", None)


inside_gnomad2_query=Variant.objects.filter(residue__flank_residue__feature_location="Inside flank", residue__protein__total_tmh_number__gte=2, residue__flank_residue__flank__tmh__meta_tmh=True, variant_source="gnomAD2").exclude(aa_mut=F("aa_wt")).values_list("aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id")
print(len(inside_gnomad2_query), "gnomad v2 variants in the multipass inside flank")
inside_flank_gnomad2_variants = heatmap_array(remove_duplicate_variants(list(inside_gnomad2_query)), aa_list_baezo_order)
heatmap(np.array(inside_flank_gnomad2_variants), "gnomAD v2 disease variants in multipass inside flanks", aa_list_baezo_order, "Greens", None)


#TMH
tmh_disease_query=Variant.objects.filter(residue__tmh_residue__feature_location="TMH", residue__protein__total_tmh_number__gte=2, residue__tmh_residue__tmh_id__meta_tmh=True, disease_status='d', variant_source="ClinVar").values_list("aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id")
print(len(tmh_disease_query), "disease variants in the multipasstmh")
tmh_disease_variants = heatmap_array(remove_duplicate_variants(list(tmh_disease_query)), aa_list_baezo_order)
heatmap(np.array(tmh_disease_variants), "ClinVar disease variants in multipass TMHs", aa_list_baezo_order, "Reds", None)

tmh_benign_query=Variant.objects.filter(residue__tmh_residue__feature_location="TMH", residue__protein__total_tmh_number__gte=2, residue__tmh_residue__tmh_id__meta_tmh=True, disease_status='n', variant_source="ClinVar").exclude(aa_mut=F("aa_wt")).values_list("aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id")
print(len(tmh_benign_query), "clinvar benign variants in the multipasstmh")
tmh_benign_variants = heatmap_array(remove_duplicate_variants(list(tmh_benign_query)), aa_list_baezo_order)
heatmap(np.array(tmh_benign_variants), "ClinVar benign variants in multipass TMHs", aa_list_baezo_order, "Blues", None)

tmh_gnomad3_query=Variant.objects.filter(residue__tmh_residue__feature_location="TMH", residue__protein__total_tmh_number__gte=2, residue__tmh_residue__tmh_id__meta_tmh=True, variant_source="gnomAD3").exclude(aa_mut=F("aa_wt")).values_list("aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id")
print(len(tmh_gnomad3_query), "gnomad v3 variants in multipass tmh")
tmh_gnomad3_variants = heatmap_array(remove_duplicate_variants(list(tmh_gnomad3_query)), aa_list_baezo_order)
heatmap(np.array(tmh_gnomad3_variants), "gnomAD v3 disease variants in multipass TMHs", aa_list_baezo_order, "Greens", None)

tmh_gnomad2_query=Variant.objects.filter(residue__tmh_residue__feature_location="TMH", residue__protein__total_tmh_number__gte=2, residue__tmh_residue__tmh_id__meta_tmh=True, variant_source="gnomAD2").exclude(aa_mut=F("aa_wt")).values_list("aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id")
print(len(tmh_gnomad2_query), "gnomad v2 variants in multipass tmh")
tmh_gnomad2_variants = heatmap_array(remove_duplicate_variants(list(tmh_gnomad2_query)), aa_list_baezo_order)
heatmap(np.array(tmh_gnomad2_variants), "gnomAD v2 disease variants in multipass TMHs", aa_list_baezo_order, "Greens", None)




### SINGLEPASS STARTS HERE ###

#Outside flank
single_outside_disease_query=Variant.objects.filter(residue__flank_residue__feature_location="Outside flank", residue__protein__total_tmh_number=1, residue__flank_residue__flank__tmh__meta_tmh=True, disease_status='d', variant_source="ClinVar").values_list("aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id")
print(len(single_outside_disease_query), "disease variants in the singlepass outside flank")
single_outside_flank_disease_variants = heatmap_array(remove_duplicate_variants(list(single_outside_disease_query)), aa_list_baezo_order)
heatmap(np.array(single_outside_flank_disease_variants), "ClinVar disease variants in singlepass outside flanks", aa_list_baezo_order, "Reds", None)


single_outside_gnomad3_query=Variant.objects.filter(residue__flank_residue__feature_location="Outside flank", residue__protein__total_tmh_number=1, residue__flank_residue__flank__tmh__meta_tmh=True,  variant_source="gnomAD3").exclude(aa_mut=F("aa_wt")).values_list("aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id")
print(len(single_outside_gnomad3_query), "gnomad v3 variants in the singlepass outside flank")
single_outside_flank_gnomad3_variants = heatmap_array(remove_duplicate_variants(list(single_outside_gnomad3_query)), aa_list_baezo_order)
heatmap(np.array(single_outside_flank_gnomad3_variants), "gnomAD v3 disease variants in singlepass outside flanks", aa_list_baezo_order, "Greens", None)

single_outside_gnomad2_query=Variant.objects.filter(residue__flank_residue__feature_location="Outside flank", residue__protein__total_tmh_number=1, residue__flank_residue__flank__tmh__meta_tmh=True,  variant_source="gnomAD2").exclude(aa_mut=F("aa_wt")).values_list("aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id")
print(len(single_outside_gnomad2_query), "gnomad v2 variants in the singlepass outside flank")
single_outside_flank_gnomad2_variants = heatmap_array(remove_duplicate_variants(list(single_outside_gnomad2_query)), aa_list_baezo_order)
heatmap(np.array(single_outside_flank_gnomad2_variants), "gnomAD v2 disease variants in singlepass outside flanks", aa_list_baezo_order, "Greens", None)

#Inside flank
single_inside_disease_query=Variant.objects.filter(residue__flank_residue__feature_location="Inside flank", residue__protein__total_tmh_number=1, residue__flank_residue__flank__tmh__meta_tmh=True, disease_status='d', variant_source="ClinVar").values_list("aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id")
print(len(single_inside_disease_query), "disease variants in the singlepass inside flank")
single_inside_flank_disease_variants = heatmap_array(remove_duplicate_variants(list(single_inside_disease_query)), aa_list_baezo_order)
heatmap(np.array(single_inside_flank_disease_variants), "ClinVar disease variants in singlepass inside flanks", aa_list_baezo_order, "Reds", None)


single_inside_gnomad3_query=Variant.objects.filter(residue__flank_residue__feature_location="Inside flank", residue__protein__total_tmh_number=1, residue__flank_residue__flank__tmh__meta_tmh=True, variant_source="gnomAD3").exclude(aa_mut=F("aa_wt")).values_list("aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id")
print(len(single_inside_gnomad3_query), "gnomad v3 variants in the singlepass inside flank")
single_inside_flank_gnomad3_variants = heatmap_array(remove_duplicate_variants(list(single_inside_gnomad3_query)), aa_list_baezo_order)
heatmap(np.array(single_outside_flank_gnomad3_variants), "gnomAD v3 disease variants in singlepass inside flanks", aa_list_baezo_order, "Greens", None)

single_inside_gnomad2_query=Variant.objects.filter(residue__flank_residue__feature_location="Inside flank", residue__protein__total_tmh_number=1, residue__flank_residue__flank__tmh__meta_tmh=True, variant_source="gnomAD2").exclude(aa_mut=F("aa_wt")).values_list("aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id")
print(len(single_inside_gnomad2_query), "gnomad v2 variants in the singlepass inside flank")
single_inside_flank_gnomad2_variants = heatmap_array(remove_duplicate_variants(list(single_inside_gnomad2_query)), aa_list_baezo_order)
heatmap(np.array(single_outside_flank_gnomad2_variants), "gnomAD v2 disease variants in singlepass inside flanks", aa_list_baezo_order, "Greens", None)

#TMHs
single_tmh_disease_query=Variant.objects.filter(residue__tmh_residue__feature_location="TMH", residue__protein__total_tmh_number=1, residue__tmh_residue__tmh_id__meta_tmh=True, disease_status='d', variant_source="ClinVar").values_list("aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id")
print(len(single_tmh_disease_query), "disease variants in the singlepass TMHs")
single_tmh_disease_variants = heatmap_array(remove_duplicate_variants(list(single_tmh_disease_query)), aa_list_baezo_order)
heatmap(np.array(single_tmh_disease_variants), "ClinVar disease variants in singlepass TMHs", aa_list_baezo_order, "Reds", None)

single_tmh_gnomad3_query=Variant.objects.filter(residue__tmh_residue__feature_location="TMH", residue__protein__total_tmh_number=1, residue__tmh_residue__tmh_id__meta_tmh=True, variant_source="gnomAD3").exclude(aa_mut=F("aa_wt")).values_list("aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id")
print(len(single_tmh_gnomad3_query), "gnomad v3 variants in singlepass tmhs")
single_tmh_gnomad3_variants = heatmap_array(remove_duplicate_variants(list(single_tmh_gnomad3_query)), aa_list_baezo_order)
heatmap(np.array(single_tmh_gnomad3_variants), "gnomAD v3 disease variants in singlepass TMHs", aa_list_baezo_order, "Greens", None)

single_tmh_gnomad2_query=Variant.objects.filter(residue__tmh_residue__feature_location="TMH", residue__protein__total_tmh_number=1, residue__tmh_residue__tmh_id__meta_tmh=True, variant_source="gnomAD2").exclude(aa_mut=F("aa_wt")).values_list("aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id")
print(len(single_tmh_gnomad2_query), "gnomad v2 variants in singlepass tmhs")
single_tmh_gnomad2_variants = heatmap_array(remove_duplicate_variants(list(single_tmh_gnomad2_query)), aa_list_baezo_order)
heatmap(np.array(single_tmh_gnomad2_variants), "gnomAD v3 disease variants in singlepass TMHs", aa_list_baezo_order, "Greens", None)



### Alternative stuff starts here ###
# Non-TMH helices
# I came across some TMHs labelled incorrectly as helices. They mayhave snuck into the database, so to ensure they are not counted as variants, meta-tmh exclude is needed.
helix_disease_query=Variant.objects.filter(residue__non_tmh_helix_residue__nont_tmh_helix_id__helix_start__gte=0, disease_status='d', variant_source="ClinVar").exclude(residue__tmh_residue__tmh_id__meta_tmh=True).values_list("aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id")
print(len(helix_disease_query), "disease variants in the non-TMH helix")
helix_disease_variants = heatmap_array(remove_duplicate_variants(list(helix_disease_query)), aa_list_baezo_order)
heatmap(np.array(helix_disease_variants), "ClinVar disease variants in helices", aa_list_baezo_order, "Reds", None)

helix_gnomad3_query=Variant.objects.filter(residue__non_tmh_helix_residue__nont_tmh_helix_id__helix_start__gte=0, variant_source="gnomAD3").exclude(aa_mut=F("aa_wt")).exclude(residue__tmh_residue__tmh_id__meta_tmh=True).values_list("aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id")
print(len(helix_gnomad3_query), "gnomad v3 variants in the non-TMH helix")
helix_gnomad3_variants = heatmap_array(remove_duplicate_variants(list(helix_gnomad3_query)), aa_list_baezo_order)
heatmap(np.array(helix_gnomad3_variants), "gnomAD v3 disease variants in helices", aa_list_baezo_order, "Greens", None)

helix_gnomad2_query=Variant.objects.filter(residue__non_tmh_helix_residue__nont_tmh_helix_id__helix_start__gte=0, variant_source="gnomAD2").exclude(aa_mut=F("aa_wt")).exclude(residue__tmh_residue__tmh_id__meta_tmh=True).values_list("aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id")
print(len(helix_gnomad2_query), "gnomad v2 variants in the non-TMH helix")
helix_gnomad2_variants = heatmap_array(remove_duplicate_variants(list(helix_gnomad2_query)), aa_list_baezo_order)
heatmap(np.array(helix_gnomad2_variants), "gnomAD v2 disease variants in helices", aa_list_baezo_order, "Greens", None)


# Signal peptides
sp_disease_query=Variant.objects.filter(residue__signal_residue__the_signal_peptide__signal_start__gte=0, disease_status='d', variant_source="ClinVar").values_list("aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id")
print(len(sp_disease_query), "disease variants in the signal peptides")
sp_disease_variants = heatmap_array(remove_duplicate_variants(list(sp_disease_query)), aa_list_baezo_order)
heatmap(np.array(sp_disease_variants), "ClinVar disease variants in signal_peptides", aa_list_baezo_order, "Reds", None)

sp_benign_query=Variant.objects.filter(residue__signal_residue__the_signal_peptide__signal_start__gte=0, disease_status='n', variant_source="ClinVar").exclude(aa_mut=F("aa_wt")).values_list("aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id")
print(len(sp_benign_query), "beign variants in the signal peptides")
sp_benign_variants = heatmap_array(remove_duplicate_variants(list(sp_benign_query)), aa_list_baezo_order)
heatmap(np.array(sp_benign_variants), "ClinVar benign variants in signal_peptides", aa_list_baezo_order, "Blues", None)


sp_gnomad3_query=Variant.objects.filter(residue__signal_residue__the_signal_peptide__signal_start__gte=0, variant_source="gnomAD3").exclude(aa_mut=F("aa_wt")).values_list("aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id")
print(len(sp_gnomad3_query), "gnomad v3 variants in the signal peptides")
sp_gnomad3_variants = heatmap_array(remove_duplicate_variants(list(sp_gnomad3_query)), aa_list_baezo_order)
heatmap(np.array(sp_gnomad3_variants), "gnomAD v3 disease variants in signal_peptides", aa_list_baezo_order, "Reds", None)



### Families and more specific queries ###

# Pore residues
pore_disease_query=Variant.objects.filter(residue__tmh_residue__feature_location="TMH", residue__structural_residue__pore_residue=True, residue__tmh_residue__tmh_id__meta_tmh=True, disease_status='d', variant_source="ClinVar").values_list("aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id")
print(len(pore_disease_query), "disease variants in the tmh pore residues")
pore_disease_variants = heatmap_array(remove_duplicate_variants(list(pore_disease_query)), aa_list_baezo_order)
heatmap(np.array(pore_disease_variants), "ClinVar disease variants in pore residue TMHs", aa_list_baezo_order, "Reds", None)

pore_benign_query=Variant.objects.filter(residue__tmh_residue__feature_location="TMH", residue__structural_residue__pore_residue=True, residue__tmh_residue__tmh_id__meta_tmh=True, disease_status='n', variant_source="ClinVar").exclude(aa_mut=F("aa_wt")).values_list("aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id")
print(len(pore_benign_query), "benign variants in the tmh pore residues")
pore_benign_variants = heatmap_array(remove_duplicate_variants(list(pore_benign_query)), aa_list_baezo_order)
heatmap(np.array(pore_benign_variants), "ClinVar benign variants in pore residue TMHs", aa_list_baezo_order, "Blues", None)

pore_gnomad3_query=Variant.objects.filter(residue__tmh_residue__feature_location="TMH", residue__structural_residue__pore_residue=True, residue__tmh_residue__tmh_id__meta_tmh=True, variant_source="gnomAD3").exclude(aa_mut=F("aa_wt")).values_list("aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id")
print(len(pore_gnomad3_query), "gnomad v3 variants in the TMH pore residues")
pore_gnomad3_variants = heatmap_array(remove_duplicate_variants(list(pore_gnomad3_query)), aa_list_baezo_order)
heatmap(np.array(pore_gnomad3_variants), "gnomAD v3 disease variants in tmh pore residues", aa_list_baezo_order, "Greens", None)

pore_gnomad2_query=Variant.objects.filter(residue__tmh_residue__feature_location="TMH", residue__structural_residue__pore_residue=True, residue__tmh_residue__tmh_id__meta_tmh=True, variant_source="gnomAD2").exclude(aa_mut=F("aa_wt")).values_list("aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id")
print(len(pore_gnomad2_query), "gnomad v2 variants in the helix")
pore_gnomad2_variants = heatmap_array(remove_duplicate_variants(list(pore_gnomad2_query)), aa_list_baezo_order)
heatmap(np.array(pore_gnomad2_variants), "gnomAD v2 disease variants in tmh pore residues", aa_list_baezo_order, "Greens", None)


# GPCRs
#disease_mp_tmh_variants = list(Variant.objects.filter(residue__tmh_residue__feature_location="TMH", residue__protein__keywords__keyword="G-protein coupled receptor").filter(disease_status='d').exclude(residue__protein__total_tmh_number=1).values_list("aa_wt", "aa_mut"))
#disease_mp_tmh_variants_dict = substitution_dictionary(disease_mp_tmh_variants)
#heatmap(sub_dict_to_heatmap(disease_mp_tmh_variants_dict) ,title, aa_list_baezo_order, "coolwarm", None)
#
#title = "Disease propensity in GPCRs"
#gnomad_mp_tmh_variants = list(Variant.objects.exclude(aa_mut=F("aa_wt")).filter(residue__tmh_residue__feature_location="TMH", residue__protein__keywords__keyword="G-protein coupled receptor").filter(variant_source='gnomAD').exclude(residue__protein__total_tmh_number=1).values_list("aa_wt", "aa_mut"))
#gnomad_mp_tmh_variants_dict = substitution_dictionary(gnomad_mp_tmh_variants)
#print(title, "disease:", len(disease_mp_tmh_variants), "gnomAD:", len(gnomad_mp_tmh_variants))
#mp_disease_propensity = subs_normalise_by_dic(disease_mp_tmh_variants_dict, gnomad_mp_tmh_variants_dict)
#mp_disease_propensity_array=sub_dict_to_heatmap(mp_disease_propensity)


heatmap_normalised_by_heatmap("ClinVar disease normalised by gnomad v3 meta-tmh multipass Outside flank", outside_flank_disease_variants, outside_flank_gnomad3_variants)
heatmap_normalised_by_heatmap("ClinVar disease normalised by gnomad v3 meta-tmh multipass Inside flank", inside_flank_disease_variants, inside_flank_gnomad3_variants)
heatmap_normalised_by_heatmap("ClinVar disease normalised by gnomad v3 meta-tmh multipass TMH", tmh_disease_variants, tmh_gnomad3_variants)
heatmap_normalised_by_heatmap("ClinVar disease normalised by ClinVar benign meta-tmh multipass TMH", tmh_disease_variants, tmh_benign_variants)

heatmap_normalised_by_heatmap("ClinVar disease normalised by gnomad v3 meta-tmh singlepass Outside flank", single_outside_flank_disease_variants, single_outside_flank_gnomad3_variants)
heatmap_normalised_by_heatmap("ClinVar disease normalised by gnomad v3 meta-tmh singlepass Inside flank", single_inside_flank_disease_variants, single_inside_flank_gnomad3_variants)
heatmap_normalised_by_heatmap("ClinVar disease normalised by gnomad v3 meta-tmh singlepass TMH",single_tmh_disease_variants, single_tmh_gnomad3_variants)

heatmap_normalised_by_heatmap("ClinVar disease normalised by gnomad v3 non-TMH Helix", helix_disease_variants, helix_gnomad3_variants)
heatmap_normalised_by_heatmap("ClinVar disease normalised by gnomad v3 Signal Peptides", sp_disease_variants, sp_gnomad3_variants)
heatmap_normalised_by_heatmap("ClinVar disease normalised by ClinVar benign Signal Peptides", sp_disease_variants, sp_benign_variants)
heatmap_normalised_by_heatmap("ClinVar disease normalised by ClinVar benign Pore residues", pore_disease_variants, pore_benign_variants)


# QUICK!!!
pore_disease_query = Variant.objects.filter(residue__tmh_residue__feature_location="TMH", residue__structural_residue__pore_residue=True, residue__tmh_residue__tmh_id__meta_tmh=True,
                                            disease_status='d').distinct("pk").values_list("aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id")
print(len(pore_disease_query), "disease variants in the tmh pore residues")
pore_disease_variants = heatmap_array(remove_duplicate_variants(
    list(pore_disease_query)), aa_list_baezo_order)
heatmap(np.array(pore_disease_variants),
        "ClinVar disease variants in pore residue TMHs", aa_list_baezo_order, "Reds", None)

pore_gnomad2_query = Variant.objects.filter(residue__tmh_residue__feature_location="TMH", residue__structural_residue__pore_residue=True, residue__tmh_residue__tmh_id__meta_tmh=True,
                                            variant_source="gnomAD2").distinct("pk").exclude(aa_mut=F("aa_wt")).values_list("aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id")
print(len(pore_gnomad2_query), "gnomad v2 variants in the pore residue")
pore_gnomad2_variants = heatmap_array(remove_duplicate_variants(
    list(pore_gnomad2_query)), aa_list_baezo_order)
heatmap(np.array(pore_gnomad2_variants),
        "gnomAD v2 disease variants in tmh pore residues", aa_list_baezo_order, "Greens", None)

heatmap_normalised_by_heatmap(
    "Disease variants normalised by gnomAD version 2 residues in the pore", pore_disease_variants, pore_gnomad2_variants)


tmh_disease_query = Variant.objects.filter(residue__tmh_residue__feature_location="TMH", residue__protein__total_tmh_number__gte=1, residue__tmh_residue__tmh_id__meta_tmh=True,
                                           disease_status='d').distinct("pk").values_list("aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id")
print(len(tmh_disease_query), "disease variants in the tmhs")
tmh_disease_variants = heatmap_array(remove_duplicate_variants(
    list(tmh_disease_query)), aa_list_baezo_order)
heatmap(np.array(tmh_disease_variants),
        "ClinVar disease variants in multipass TMHs", aa_list_baezo_order, "Reds", None)

tmh_gnomad2_query = Variant.objects.filter(residue__tmh_residue__feature_location="TMH", residue__protein__total_tmh_number__gte=1, residue__tmh_residue__tmh_id__meta_tmh=True,
                                           variant_source='gnomAD2').distinct("pk").values_list("aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id")
print(len(tmh_gnomad2_query), "disease variants in the tmhs")
tmh_gnomad2_variants = heatmap_array(remove_duplicate_variants(
    list(tmh_gnomad2_query)), aa_list_baezo_order)
heatmap(np.array(tmh_gnomad2_variants),
        "ClinVar disease variants in multipass TMHs", aa_list_baezo_order, "Greens", None)

heatmap_normalised_by_heatmap(
    "Disease variants normalised by gnomAD version 2 residues in TMHs", tmh_disease_variants, tmh_gnomad2_variants)
'''
