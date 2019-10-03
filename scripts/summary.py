# Shell Plus Model Imports
from django.contrib.admin.models import LogEntry
from django.contrib.auth.models import Group, Permission, User
from django.contrib.contenttypes.models import ContentType
from django.contrib.sessions.models import Session
from tmh_db.models import Binding_residue, Database_Metadata, Flank, Funfam, Funfam_residue, Funfamstatus, Go, Keyword, Pfam, Pfam_residue, Protein, Residue, Structural_residue, Structure, Subcellular_location, Tmh, Tmh_deltag, Tmh_hydrophobicity, Tmh_residue, Tmh_tmsoc, Uniref, Variant
# Shell Plus Django Imports
from django.core.cache import cache
from django.conf import settings
from django.contrib.auth import get_user_model
from django.db import transaction
from django.db.models import Avg, Case, Count, F, Max, Min, Prefetch, Q, Sum, When, Exists, OuterRef, Subquery
from django.utils import timezone
from django.urls import reverse
# Charts
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
import collections
from matplotlib.colors import LogNorm
from scripts.graphs import *

from scripts.populate_general_functions import *


def normalise_tmh_resid_array(var_freqs_list, aa_list):
    large_array = []
    freq_dict={}

    for order_number, aa in enumerate(aa_list):
        residue_count=Residue.objects.exclude(tmh_residue=None).filter(amino_acid_type=aa).count()
        freq_dict[aa]=residue_count

    for aa_mut_order_number, aa_mut_row in enumerate(var_freqs_list):
        aa_mut_row_normalised = []
        for aa_wt_order_number, aa_wt_count in enumerate(aa_mut_row):
            aa_mut_row_normalised.append(aa_wt_count/freq_dict[aa_list[aa_wt_order_number]])
        large_array.append(aa_mut_row_normalised)
    #print(np.array(large_array))
    return(np.array(large_array))


def heatmap_run():
    print("HEATMAPS")
    ### HEATMAPS ###
    aa_list_alpha=['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', 'X']

    aa_list_baezo_order=['K', 'R', 'E', 'D', 'Q', 'H', 'N', 'P', 'Y', 'W', 'C', 'M', 'T', 'S', 'G', 'V', 'F', 'A', 'I', 'L']
    title = "TMH±5_disease"

    Residue.objects.all().prefetch_related("variant")

    var_freq=list(Variant.objects.exclude(residue__tmh_residue=None).filter(disease_status='d').values_list("aa_wt", "aa_mut"))
    var_freqs_list=heatmap_array(var_freq, aa_list_baezo_order)
    heatmap(var_freqs_list, title, aa_list_baezo_order, "coolwarm", None)
    title = "TMH±5_disease_normalised_wt_aa_tmh_freq"
    heatmap(normalise_tmh_resid_array(var_freqs_list, aa_list_baezo_order), title, aa_list_baezo_order, "reds", None)

    title = "non_TMH±5_disease"
    var_freq=list(Variant.objects.filter(residue__tmh_residue=None).filter(disease_status='d').values_list("aa_wt", "aa_mut"))
    heatmap(heatmap_array(var_freq, aa_list_baezo_order), title, aa_list_baezo_order, "reds", None)

    title = "TMH_disease_clinvar"
    clinvar_tmh_var_freq=Variant.objects.filter(residue__tmh_residue__feature_location="TMH", variant_source="ClinVar", disease_status='d').values_list("aa_wt", "aa_mut")
    print(title, clinvar_tmh_var_freq.count())
    clinvar_tmh_disease_var_array = heatmap_array(list(clinvar_tmh_var_freq), aa_list_baezo_order)
    heatmap(clinvar_tmh_disease_var_array, title, aa_list_baezo_order , "reds", None)

    title = "Flank_disease_clinvar"
    clinvar_flank_var_freq=Variant.objects.exclude(residue__tmh_residue=None).exclude(residue__tmh_residue__feature_location="TMH").filter(variant_source="ClinVar", disease_status='d').values_list("aa_wt", "aa_mut")
    print(title, clinvar_flank_var_freq.count())
    clinvar_flank_disease_var_array = heatmap_array(list(clinvar_flank_var_freq), aa_list_baezo_order)
    heatmap(clinvar_flank_disease_var_array, title, aa_list_baezo_order , "reds", None)

    title = "non_TMH±5_disease_clinvar"
    clinvar_non_tmh_var_freq=Variant.objects.filter(residue__tmh_residue=None).filter(variant_source="ClinVar", disease_status='d').values_list("aa_wt", "aa_mut")
    print(title, clinvar_non_tmh_var_freq.count())
    clinvar_non_tmh_disease_var_array = heatmap_array(list(clinvar_non_tmh_var_freq), aa_list_baezo_order)
    heatmap(clinvar_non_tmh_disease_var_array, title, aa_list_baezo_order, "reds", None)

    title = "TMH±5_disease_humsavar"
    var_freq=list(Variant.objects.exclude(residue__tmh_residue=None).filter(variant_source="Humsavar", disease_status='d').values_list("aa_wt", "aa_mut"))
    heatmap(heatmap_array(var_freq, aa_list_baezo_order), title, aa_list_baezo_order, "reds", None)

    title = "non_TMH±5_disease_humsavar"
    var_freq=list(Variant.objects.filter(residue__tmh_residue=None).filter(variant_source="Humsavar", disease_status='d').values_list("aa_wt", "aa_mut"))
    heatmap(heatmap_array(var_freq, aa_list_baezo_order), title, aa_list_baezo_order, "reds", None)

    title = "TMH_gnomAD"
    gnomad_tmh_var_freq=Variant.objects.exclude(aa_mut=F("aa_wt")).filter(residue__tmh_residue__feature_location="TMH").filter(variant_source='gnomAD').values_list("aa_wt", "aa_mut")
    print(title, gnomad_tmh_var_freq.count())
    gnomad_tmh_var_array = heatmap_array(gnomad_tmh_var_freq, aa_list_baezo_order)
    heatmap(gnomad_tmh_var_array, title, aa_list_baezo_order, "n", None)

    title = "Flank_gnomAD"
    gnomad_flank_var_freq=Variant.objects.exclude(aa_mut=F("aa_wt")).exclude(residue__tmh_residue=None).exclude(residue__tmh_residue__feature_location="TMH").filter(variant_source='gnomAD').values_list("aa_wt", "aa_mut")
    print(title, gnomad_flank_var_freq.count())
    gnomad_flank_var_array = heatmap_array(gnomad_flank_var_freq, aa_list_baezo_order)
    heatmap(gnomad_flank_var_array, title, aa_list_baezo_order, "n", None)

    title = "non_TMH±5_gnomAD"
    gnomad_non_tmh_var_freq=Variant.objects.exclude(aa_mut=F("aa_wt")).filter(residue__tmh_residue=None).filter(variant_source='gnomAD').values_list("aa_wt", "aa_mut")
    print(title,gnomad_non_tmh_var_freq.count())
    gnomad_non_tmh_var_array=heatmap_array(gnomad_non_tmh_var_freq, aa_list_baezo_order)
    heatmap(gnomad_non_tmh_var_array, title, aa_list_baezo_order, "n", None)

    # bonkers contingency table idea #
    p_matrix = []
    for var_aa_num, aa_var in enumerate(aa_list_baezo_order):
        p_list = []
        for wt_aa_num, aa_wt in enumerate(aa_list_baezo_order):
            fisher_oddsratio, fisher_pvalue = stats.fisher_exact([[clinvar_tmh_disease_var_array[var_aa_num][wt_aa_num], clinvar_non_tmh_disease_var_array[var_aa_num][wt_aa_num]], [gnomad_tmh_var_array[var_aa_num][wt_aa_num], gnomad_non_tmh_var_array[var_aa_num][wt_aa_num]]], alternative='two-sided')
            p_list.append(fisher_pvalue)
        p_matrix.append(p_list)
    heatmap(np.array(p_matrix), "fisher_pvalue_[clin_tmh±5,clin_nontmh±5][gnom_tmh±5,gnom_nontmh±5]", aa_list_baezo_order, "s", LogNorm()) #needs a better scale

    # normalised to residue count presence #
    residue_count_dict={}
    print("Residue counts that fall outside TMH±5 in TMPs")
    for aa in aa_list_baezo_order:
        aa_count = Residue.objects.filter(tmh_residue=None, amino_acid_type=aa).count()
        print(aa, aa_count)
        residue_count_dict[aa]=aa_count

    tmh_residue_count_dict={}
    print("Residue counts that fall within the TMH in TMPs")
    for aa in aa_list_baezo_order:
        aa_count = Residue.objects.filter(tmh_residue__feature_location="TMH", amino_acid_type=aa).count()
        print(aa, aa_count)
        tmh_residue_count_dict[aa]=aa_count

    color_dict_list = {
    "gnomAD" : "n",
    "ClinVar": "reds"
    }

    #Normalise according to residue populations

    tmh_datasets= [gnomad_tmh_var_array, clinvar_tmh_disease_var_array]
    for n, dataset in enumerate(tmh_datasets):
        residue_normalised_count_matrix = []
        for var_aa_num, aa_var in enumerate(aa_list_baezo_order):
            residue_normalised_list = []
            for wt_aa_num, aa_wt in enumerate(aa_list_baezo_order):
                #print(dataset[var_aa_num][wt_aa_num], tmh_residue_count_dict[aa_wt])
                residue_normalised_value = dataset[var_aa_num][wt_aa_num]/tmh_residue_count_dict[aa_wt]
                residue_normalised_list.append(residue_normalised_value)
            residue_normalised_count_matrix.append(residue_normalised_list)
        if n == 1:
            source="ClinVar"
        elif n == 0:
            source="gnomAD"

        heatmap(np.array(residue_normalised_count_matrix), f"Residue normalised according to WT residue in TMH residue population in {color_dict_list[source]} state", aa_list_baezo_order, color_dict_list[source], None) #needs a better scale


    non_tmh_datasets= [gnomad_non_tmh_var_array, clinvar_non_tmh_disease_var_array]
    for n, dataset in enumerate(non_tmh_datasets):
        residue_normalised_count_matrix = []
        for var_aa_num, aa_var in enumerate(aa_list_baezo_order):
            residue_normalised_list = []
            for wt_aa_num, aa_wt in enumerate(aa_list_baezo_order):
                residue_normalised_value = dataset[var_aa_num][wt_aa_num]/residue_count_dict[aa_wt]
                residue_normalised_list.append(residue_normalised_value)
            residue_normalised_count_matrix.append(residue_normalised_list)
        if n == 1:
            source="ClinVar"
        elif n == 0:
            source="gnomAD"

        heatmap(np.array(residue_normalised_count_matrix), f"Residue normalised according to WT residue in non_TMH residue population in {color_dict_list[source]} state", aa_list_baezo_order, color_dict_list[source], None) #needs a better scale


    #Normalise ClinVar acording to gnomAD in TMHs
    residue_normalised_count_matrix = []
    for var_aa_num, aa_var in enumerate(aa_list_baezo_order):
        residue_normalised_list = []
        for wt_aa_num, aa_wt in enumerate(aa_list_baezo_order):

            #print(dataset[var_aa_num][wt_aa_num], tmh_residue_count_dict[aa_wt])
            normalised_gnomad_value=gnomad_tmh_var_array[var_aa_num][wt_aa_num]

            if var_aa_num==wt_aa_num:
                residue_normalised_value = 0

            elif normalised_gnomad_value == 0:
                residue_normalised_value = clinvar_tmh_disease_var_array[var_aa_num][wt_aa_num]

            else:
                residue_normalised_value = clinvar_tmh_disease_var_array[var_aa_num][wt_aa_num]/gnomad_tmh_var_array[var_aa_num][wt_aa_num]


            residue_normalised_list.append(residue_normalised_value)
        residue_normalised_count_matrix.append(residue_normalised_list)
    source="ClinVar"
    heatmap(np.array(residue_normalised_count_matrix), f"Residue normalised according to ClinVar div gnomAd in TMH residue population in {color_dict_list[source]} state", aa_list_baezo_order, color_dict_list[source], None)




    #Normalise ClinVar acording to gnomAD in flanks
    residue_normalised_count_matrix = []
    for var_aa_num, aa_var in enumerate(aa_list_baezo_order):
        residue_normalised_list = []
        for wt_aa_num, aa_wt in enumerate(aa_list_baezo_order):

            #print(dataset[var_aa_num][wt_aa_num], tmh_residue_count_dict[aa_wt])
            normalised_gnomad_value=gnomad_flank_var_array[var_aa_num][wt_aa_num]

            if var_aa_num==wt_aa_num:
                residue_normalised_value = 0

            elif normalised_gnomad_value == 0:
                residue_normalised_value = clinvar_flank_disease_var_array[var_aa_num][wt_aa_num]

            else:
                residue_normalised_value = clinvar_flank_disease_var_array[var_aa_num][wt_aa_num]/gnomad_flank_var_array[var_aa_num][wt_aa_num]


            residue_normalised_list.append(residue_normalised_value)
        residue_normalised_count_matrix.append(residue_normalised_list)
    source="ClinVar"
    heatmap(np.array(residue_normalised_count_matrix), f"Residue normalised according to ClinVar div gnomAd in flank residue population in {color_dict_list[source]} state", aa_list_baezo_order, color_dict_list[source], None)




    #Normalise ClinVar acording to gnomAD in non-TMHs
    residue_normalised_count_matrix = []
    for var_aa_num, aa_var in enumerate(aa_list_baezo_order):
        residue_normalised_list = []
        for wt_aa_num, aa_wt in enumerate(aa_list_baezo_order):

            #print(dataset[var_aa_num][wt_aa_num], tmh_residue_count_dict[aa_wt])
            normalised_gnomad_value=gnomad_non_tmh_var_array[var_aa_num][wt_aa_num]

            if var_aa_num==wt_aa_num:
                residue_normalised_value = 0

            elif normalised_gnomad_value == 0:
                residue_normalised_value = clinvar_non_tmh_disease_var_array[var_aa_num][wt_aa_num]

            else:
                residue_normalised_value = clinvar_non_tmh_disease_var_array[var_aa_num][wt_aa_num]/gnomad_non_tmh_var_array[var_aa_num][wt_aa_num]


            residue_normalised_list.append(residue_normalised_value)
        residue_normalised_count_matrix.append(residue_normalised_list)
    source="ClinVar"
    heatmap(np.array(residue_normalised_count_matrix), f"Residue normalised according to ClinVar div gnomAd in non_TMH residue population in {color_dict_list[source]} state", aa_list_baezo_order, color_dict_list[source], None)


    ##Normalise ClinVar acording to gnomAD and get a p-value
    #residue_normalised_count_matrix = []
    #for var_aa_num, aa_var in enumerate(aa_list_baezo_order):
    #    residue_normalised_list = []
    #    for wt_aa_num, aa_wt in enumerate(aa_list_baezo_order):

    #        #print(dataset[var_aa_num][wt_aa_num], tmh_residue_count_dict[aa_wt])
    #        normalised_gnomad_value=gnomad_tmh_var_array[var_aa_num][wt_aa_num]

    #        if var_aa_num==wt_aa_num:
    #            residue_normalised_value = 1 #fake test data

    #        elif normalised_gnomad_value == 0:
    #            residue_normalised_value = 1

    #        else:
    #            stats_oddsratio, stats_pvalue = stats.binom_test(clinvar_tmh_disease_var_array[var_aa_num][wt_aa_num], n=gnomad_tmh_var_array[var_aa_num][wt_aa_num], p=sum(sum(np.array(clinvar_tmh_disease_var_array)))/sum(sum(np.array(gnomad_tmh_var_array))), alternative='two-sided')
    #            residue_normalised_value = stats_pvalue
    #        print(residue_normalised_value)

    #        residue_normalised_list.append(residue_normalised_value)

    #    residue_normalised_count_matrix.append(residue_normalised_list)

    #source="ClinVar"
    #heatmap(np.array(residue_normalised_count_matrix), f"Residue normalised according to ClinVar gnomAdchi sqaure in TMH residue population in {color_dict_list[source]} state", aa_list_baezo_order, color_dict_list[source], LogNorm())


def advanced_heatmaps():
    #This could be any order
    aa_list_baezo_order=['K', 'R', 'E', 'D', 'Q', 'H', 'N', 'P', 'Y', 'W', 'C', 'M', 'T', 'S', 'G', 'V', 'F', 'A', 'I', 'L']

    #prefetch the residues with variants
    Residue.objects.all().prefetch_related("variant")


    variants_for_heatmap=list(Variant.objects.exclude(residue__tmh_residue__feature_location="TMH").filter(disease_status='d').values_list("aa_wt", "aa_mut"))
    subsitutions_in_tmhs=substitution_dictionary(variants_for_heatmap)

def substitution_dictionary(variant_list):
    aa_list_alpha=['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', 'X']
    subs={}
    for reference_amino_acid in aa_list_alpha:
        for mutant_amino_acid in aa_list_alpha:
            #print(variant_list)
            instances=variant_list.count((reference_amino_acid, mutant_amino_acid))
            subs[(reference_amino_acid,mutant_amino_acid)]=instances
    return(subs)





def basic_num():
    print("STATS OVERVIEW\n",)
    # Tmh.objects.exclude(residue__tmh_residue=None).filter(disease_status='d').count()
    print("\n\nResidues\n")
    protein_num = Protein.objects.count()
    print("UniProt IDs,", protein_num)
    residue_num = Residue.objects.count()
    print("Residues,", residue_num)
    tmh_residue_num = Residue.objects.filter(tmh_residue__feature_location="TMH").count()
    print("TMH residues,", tmh_residue_num)
    flank_residue_num = Residue.objects.filter(Q(flank_residue__feature_location="Inside flank") | Q(flank_residue__feature_location="Outside flank")).count()
    print("flank residues,", flank_residue_num)

    non_tmh_residue_num = Residue.objects.filter(tmh_residue=None).count()
    print("Non-TMH residues,", non_tmh_residue_num)
    structural_residue_num = Residue.objects.exclude(structural_residue=None).count()
    print("Residues with a map to at least one residue in a structure,", structural_residue_num)
    non_structural_residue_num = Residue.objects.filter(structural_residue=None).count()
    print("Residues with no map to a residue in a structure,",non_structural_residue_num)
    tmh_structural_residue_num = Residue.objects.exclude(tmh_residue=None, structural_residue=None).count()
    print("TMH residues with a map to at least one residue in a structure,",tmh_structural_residue_num)
    tmh_non_structural_residue_num = Residue.objects.exclude(tmh_residue=None).filter(structural_residue=None).count()
    print("TMH residues with no map to a residue in a structure,",tmh_non_structural_residue_num)

    print("\n\nTMH boundaries\n")
    tmh_boundary_num = Tmh.objects.count()
    print("TMP TMH boundaries,", tmh_boundary_num)
    uniprot_tmh_boundary = Tmh.objects.filter(tmh_evidence='UniProt').count()
    print("UniProt TMH boundaries,", uniprot_tmh_boundary)
    topdb_tmh_boundary = Tmh.objects.filter(tmh_evidence='TOPDB').count()
    print("TopDB TMH boundaries,", topdb_tmh_boundary)
    mptopo_tmh_boundary = Tmh.objects.filter(tmh_evidence='MPTOPO').count()
    print("MPTOPO TMH boundaries,", mptopo_tmh_boundary)

    print("\n\nVariants\n")
    # Complex query example. How many variants are in the TMH?
    g_variants_num = Variant.objects.filter(variant_source='gnomAD').count()
    print("gnomAD variants,", g_variants_num)
    tmh_g_variants_num = Variant.objects.filter(residue__tmh_residue__feature_location = "TMH", variant_source='gnomAD').count()
    print("TMH gnomAD variants,", tmh_g_variants_num)
    flank_g_variants_num = Variant.objects.exclude(residue__tmh_residue=None, residue__tmh_residue__feature_location = "TMH").filter(variant_source='gnomAD').count()
    print("Flank gnomAD variants,", flank_g_variants_num)
    non_tmh_g_variants_num = Variant.objects.filter(residue__tmh_residue=None, variant_source='gnomAD').count()
    print("Non-TMH gnomAD variants,", non_tmh_g_variants_num)

    d_variants_num = Variant.objects.filter(disease_status='d').count()
    print("Disease variants,", d_variants_num)
    tmh_d_variants_num = Variant.objects.filter(residue__tmh_residue__feature_location = "TMH", disease_status='d').count()
    print("TMH disease variants,", tmh_d_variants_num)
    flank_d_variants_num = Variant.objects.exclude(residue__tmh_residue=None, residue__tmh_residue__feature_location = "TMH").filter(disease_status='d').count()
    print("Flank disease variants,", flank_d_variants_num)
    d_variants_clinvar_num = Variant.objects.filter(disease_status='d', variant_source="ClinVar").count()
    print("ClinVar disease variants,", d_variants_clinvar_num)
    tmh_d_variants_clinvar_num = Variant.objects.filter(residue__tmh_residue__feature_location = "TMH", disease_status='d', variant_source="ClinVar").count()
    print("TMH ClinVar disease variants,", tmh_d_variants_clinvar_num)
    flank_d_variants_clinvar_num = Variant.objects.exclude(residue__tmh_residue=None, residue__tmh_residue__feature_location = "TMH").filter(disease_status='d', variant_source="ClinVar").count()
    print("Flank ClinVar disease variants,", flank_d_variants_num)

    flank_d_variants_clinvar_num = Variant.objects.exclude(residue__tmh_residue=None, residue__tmh_residue__feature_location = "TMH").filter(disease_status='d', variant_source="Humsavar").count()
    print("Flank Humsavar disease variants,", flank_d_variants_num)

    tmh_d_variants_humsavar_num = Variant.objects.exclude(residue__tmh_residue=None).filter(disease_status='d', variant_source="Humsavar").count()

    d_variants_humsavar_num = Variant.objects.filter(disease_status='d', variant_source="Humsavar").count()
    print("Humsavar disease variants,", d_variants_humsavar_num)
    print("TMH Humsavar disease variants,", tmh_d_variants_humsavar_num)
    non_tmh_d_variants_num = Variant.objects.filter(residue__tmh_residue=None).filter(disease_status='d').count()
    print("Non-TMH disease variants,", non_tmh_d_variants_num)
    non_tmh_d_variants_clinvar_num = Variant.objects.filter(residue__tmh_residue=None, disease_status='d', variant_source="ClinVar").count()
    print("Non-TMH ClinVar disease variants,", non_tmh_d_variants_clinvar_num)
    non_tmh_d_variants_humsavar_num = Variant.objects.filter(residue__tmh_residue=None, disease_status='d', variant_source="Humsavar").count()
    print("Non-TMH Humsavar disease variants,", non_tmh_d_variants_humsavar_num)



    disease_status = {"reds":"disease", "n": "benign"}
    for state in ["reds", "n"]:
        residues = Residue.objects.filter(variant__disease_status=state, tmh_residue__feature_location="TMH").prefetch_related("protein")
        sites = Residue.objects.filter(variant__disease_status=state, tmh_residue__feature_location="TMH").distinct().count()
        proteins = Protein.objects.filter(residue__in=residues).distinct().count()



        print("Proteins with",disease_status[state],"variants in the TMH", proteins)
        print("Residues in the TMH that have",disease_status[state],"variants", sites)
        id_tmh_var_dict=collections.Counter([i.protein.uniprot_id for i in residues])

        performance = list(id_tmh_var_dict.values())
        #print(objects)
        #print(performance)
        histogram(performance, f"Frequency of proteins with TMH {disease_status[state]} variants", state,"Number of TMH disease variants", "Number of proteins")


        #proteins =  Protein.objects.filter(residue__in=residues).distinct()
        #location_tmh_var_dict=[]
        #for i in proteins:
        #    multi_maps=i.subcellular_locations
        #    for a_location in multi_maps:
        #        location_tmh_var_dict.append(a_location.location)

        #collections.Counter(location_tmh_var_dict)
        #print(location_tmh_var_dict)
        #performance = list(location_tmh_var_dict.values())
        #objects = list(location_tmh_var_dict.keys())
        #barchart(objects, performance, f"Frequency of TMH {disease_status[state]} variants in subceullar locations", state, "Subcellular Location", "Number of variants")

    #Variant.objects.filter(residue__tmh_residue__feature_location="TMH", disease_status='d', variant_source="ClinVar").count()

    print("\n\nVariant enrichment\n")
    print("Disease variants per residue,", d_variants_num / residue_num)
    print("Disease variants per TMH residue,", tmh_d_variants_num / tmh_residue_num)
    print("Disease variants per non-TMH residue,", non_tmh_d_variants_num / non_tmh_residue_num)

    fisher_oddsratio, fisher_pvalue = stats.fisher_exact([[tmh_d_variants_num, non_tmh_d_variants_num], [tmh_g_variants_num, non_tmh_g_variants_num]], alternative='two-sided')


    print("Fisher test Disease TMH non-TMH versus gnomAD TMH non-TMH,", fisher_pvalue)

    obs = np.array([[tmh_d_variants_num, non_tmh_d_variants_num], [tmh_g_variants_num, non_tmh_g_variants_num]])
    chi_test = stats.chi2_contingency(obs)

    print("Chi 2 test Disease TMH non-TMH versus gnomAD TMH non-TMH,", chi_test)

    objects = ("Residues", "TMH",  "±5 residues", "Non-TMH residues")
    performance = [d_variants_num/residue_num, tmh_d_variants_num/tmh_residue_num, flank_d_variants_num/flank_residue_num, non_tmh_d_variants_num / non_tmh_residue_num]
    barchart(objects, performance, "TMP_disease_variants", "reds", "Residue type", "Variants per residue")

    objects = ("Residues", "TMH",  "±5 residues", "Non-TMH residues")
    performance = [d_variants_clinvar_num/residue_num, tmh_d_variants_clinvar_num/tmh_residue_num, flank_d_variants_clinvar_num/flank_residue_num, non_tmh_d_variants_clinvar_num / non_tmh_residue_num]
    barchart(objects, performance, "ClinVar_disease_variants", "reds", "Residue type", "Variants per residue")

    #objects = ("Residues", "TMH",  "±5 residues", "Non-TMH residues")
    #performance = [d_variants_humsavar_num / residue_num, tmh_d_variants_humsavar_num/tmh_residue_num, non_tmh_d_variants_humsavar_num / non_tmh_residue_num]
    #barchart(objects, performance, "Humsavar_disease_variants", "reds", "Residue type", "Variants per residue")

    objects = ("Residues", "TMH",  "±5 residues", "Non-TMH residues")
    performance = [g_variants_num/residue_num, tmh_g_variants_num/tmh_residue_num, flank_g_variants_num/flank_residue_num, non_tmh_g_variants_num/non_tmh_residue_num]
    barchart(objects, performance, "TMP_gnomAD_variants", "n", "Residue type", "Variants per residue")

def dayhoff_matrix_normalisation(matrix):
    thornton_matrix = {
        ("A", "A"): 98950,
        ("A", "R"): 21,
        ("A", "N"): 2,
        ("A", "D"): 7,
        ("A", "C"): 13,
        ("A", "Q"): 4,
        ("A", "E"): 6,
        ("A", "G"): 157,
        ("A", "H"): 6,
        ("A", "I"): 43,
        ("A", "L"): 42,
        ("A", "K"): 5,
        ("A", "M"): 10,
        ("A", "F"): 21,
        ("A", "P"): 33,
        ("A", "S"): 194,
        ("A", "T"): 198,
        ("A", "W"): 0,
        ("A", "Y"): 1,
        ("A", "V"): 287,
        # R
        ("R", "A"): 138,
        ("R", "R"): 98590,
        ("R", "N"): 0,
        ("R", "D"): 7,
        ("R", "C"): 13,
        ("R", "Q"): 138,
        ("R", "E"): 20,
        ("R", "G"): 145,
        ("R", "H"): 138,
        ("R", "I"): 26,
        ("R", "L"): 53,
        ("R", "K"): 349,
        ("R", "M"): 125,
        ("R", "F"): 0,
        ("R", "P"): 7,
        ("R", "S"): 33,
        ("R", "T"): 33,
        ("R", "W"): 184,
        ("R", "Y"): 0,
        ("R", "V"): 0,
        # N
        ("N", "A"): 11,
        ("N", "R"): 0,
        ("N", "N"): 99369,
        ("N", "D"): 78,
        ("N", "C"): 6,
        ("N", "Q"): 39,
        ("N", "E"): 0,
        ("N", "G"): 0,
        ("N", "H"): 45,
        ("N", "I"): 22,
        ("N", "L"): 28,
        ("N", "K"): 62,
        ("N", "M"): 17,
        ("N", "F"): 6,
        ("N", "P"): 11,
        ("N", "S"): 179,
        ("N", "T"): 106,
        ("N", "W"): 6,
        ("N", "Y"): 6,
        ("N", "V"): 11,
        # D
        ("D", "A"): 81,
        ("D", "R"): 12,
        ("D", "N"): 162,
        ("D", "D"): 99200,
        ("D", "C"): 0,
        ("D", "Q"): 0,
        ("D", "E"): 139,
        ("D", "G"): 174,
        ("D", "H"): 46,
        ("D", "I"): 12,
        ("D", "L"): 0,
        ("D", "K"): 23,
        ("D", "M"): 12,
        ("D", "F"): 0,
        ("D", "P"): 12,
        ("D", "S"): 0,
        ("D", "T"): 70,
        ("D", "W"): 0,
        ("D", "Y"): 12,
        ("D", "V"): 46,

        ("C", "A"): 61,
        ("C", "R"): 9,
        ("C", "N"): 5,
        ("C", "D"): 0,
        ("C", "C"): 98964,
        ("C", "Q"): 0,
        ("C", "E"): 0,
        ("C", "G"): 61,
        ("C", "H"): 9,
        ("C", "I"): 12,
        ("C", "L"): 52,
        ("C", "K"): 0,
        ("C", "M"): 5,
        ("C", "F"): 160,
        ("C", "P"): 0,
        ("C", "S"): 226,
        ("C", "T"): 61,
        ("C", "W"): 38,
        ("C", "Y"): 108,
        ("C", "V"): 221,

        ("Q", "A"): 29,
        ("Q", "R"): 154,
        ("Q", "N"): 51,
        ("Q", "D"): 0,
        ("Q", "C"): 0,
        ("Q", "Q"): 99158,
        ("Q", "E"): 117,
        ("Q", "G"): 7,
        ("Q", "H"): 190,
        ("Q", "I"): 7,
        ("Q", "L"): 117,
        ("Q", "K"): 44,
        ("Q", "M"): 22,
        ("Q", "F"): 0,
        ("Q", "P"): 37,
        ("Q", "S"): 0,
        ("Q", "T"): 51,
        ("Q", "W"): 0,
        ("Q", "Y"): 0,
        ("Q", "V"): 0,

        ("E", "A"): 64,
        ("E", "R"): 32,
        ("E", "N"): 0,
        ("E", "D"): 128,
        ("E", "C"): 0,
        ("E", "Q"): 171,
        ("E", "E"): 99241,
        ("E", "G"): 225,
        ("E", "H"): 0,
        ("E", "I"): 0,
        ("E", "L"): 0,
        ("E", "K"): 0,
        ("E", "M"): 0,
        ("E", "F"): 0,
        ("E", "P"): 0,
        ("E", "S"): 43,
        ("E", "T"): 21,
        ("E", "W"): 0,
        ("E", "Y"): 0,
        ("E", "V"): 75,

        ("G", "A"): 218,
        ("G", "R"): 30,
        ("G", "N"): 0,
        ("G", "D"): 20,
        ("G", "C"): 18,
        ("G", "Q"): 1,
        ("G", "E"): 29,
        ("G", "G"): 99468,
        ("G", "H"): 1,
        ("G", "I"): 14,
        ("G", "L"): 0,
        ("G", "K"): 0,
        ("G", "M"): 4,
        ("G", "F"): 5,
        ("G", "P"): 10,
        ("G", "S"): 87,
        ("G", "T"): 16,
        ("G", "W"): 7,
        ("G", "Y"): 0,
        ("G", "V"): 72,

        ("H", "A"): 37,
        ("H", "R"): 129,
        ("H", "N"): 49,
        ("H", "D"): 25,
        ("H", "C"): 12,
        ("H", "Q"): 160,
        ("H", "E"): 0,
        ("H", "G"): 6,
        ("H", "H"): 99329,
        ("H", "I"): 18,
        ("H", "L"): 12,
        ("H", "K"): 0,
        ("H", "M"): 6,
        ("H", "F"): 0,
        ("H", "P"): 0,
        ("H", "S"): 0,
        ("H", "T"): 25,
        ("H", "W"): 0,
        ("H", "Y"): 178,
        ("H", "V"): 12,

        ("I", "A"): 38,
        ("I", "R"): 3,
        ("I", "N"): 3,
        ("I", "D"): 1,
        ("I", "C"): 3,
        ("I", "Q"): 1,
        ("I", "E"): 0,
        ("I", "G"): 9,
        ("I", "H"): 3,
        ("I", "I"): 98579,
        ("I", "L"): 237,
        ("I", "K"): 0,
        ("I", "M"): 140,
        ("I", "F"): 57,
        ("I", "P"): 3,
        ("I", "S"): 19,
        ("I", "T"): 130,
        ("I", "W"): 1,
        ("I", "Y"): 3,
        ("I", "V"): 767,

        ("L", "A"): 27,
        ("L", "R"): 5,
        ("L", "N"): 3,
        ("L", "D"): 0,
        ("L", "C"): 7,
        ("L", "Q"): 10,
        ("L", "E"): 0,
        ("L", "G"): 0,
        ("L", "H"): 1,
        ("L", "I"): 172,
        ("L", "L"): 99274,
        ("L", "K"): 1,
        ("L", "M"): 97,
        ("L", "F"): 158,
        ("L", "P"): 23,
        ("L", "S"): 27,
        ("L", "T"): 16,
        ("L", "W"): 13,
        ("L", "Y"): 4,
        ("L", "V"): 161,

        ("K", "A"): 46,
        ("K", "R"): 487,
        ("K", "N"): 101,
        ("K", "D"): 18,
        ("K", "C"): 0,
        ("K", "Q"): 55,
        ("K", "E"): 0,
        ("K", "G"): 0,
        ("K", "H"): 0,
        ("K", "I"): 0,
        ("K", "L"): 9,
        ("K", "K"): 99164,
        ("K", "M"): 37,
        ("K", "F"): 0,
        ("K", "P"): 0,
        ("K", "S"): 9,
        ("K", "T"): 18,
        ("K", "W"): 0,
        ("K", "Y"): 46,
        ("K", "V"): 9,

        ("M", "A"): 31,
        ("M", "R"): 59,
        ("M", "N"): 9,
        ("M", "D"): 3,
        ("M", "C"): 3,
        ("M", "Q"): 9,
        ("M", "E"): 0,
        ("M", "G"): 9,
        ("M", "H"): 3,
        ("M", "I"): 499,
        ("M", "L"): 475,
        ("M", "K"): 12,
        ("M", "M"): 98465,
        ("M", "F"): 25,
        ("M", "P"): 0,
        ("M", "S"): 3,
        ("M", "T"): 99,
        ("M", "W"): 3,
        ("M", "Y"): 16,
        ("M", "V"): 276,

        ("F", "A"): 28,
        ("F", "R"): 0,
        ("F", "N"): 1,
        ("F", "D"): 0,
        ("F", "C"): 45,
        ("F", "Q"): 0,
        ("F", "E"): 0,
        ("F", "G"): 5,
        ("F", "H"): 0,
        ("F", "I"): 88,
        ("F", "L"): 333,
        ("F", "K"): 0,
        ("F", "M"): 11,
        ("F", "F"): 99311,
        ("F", "P"): 0,
        ("F", "S"): 42,
        ("F", "T"): 12,
        ("F", "W"): 3,
        ("F", "Y"): 72,
        ("F", "V"): 49,

        ("P", "A"): 135,
        ("P", "R"): 4,
        ("P", "N"): 8,
        ("P", "D"): 4,
        ("P", "C"): 0,
        ("P", "Q"): 20,
        ("P", "E"): 0,
        ("P", "G"): 28,
        ("P", "H"): 0,
        ("P", "I"): 16,
        ("P", "L"): 147,
        ("P", "K"): 0,
        ("P", "M"): 0,
        ("P", "F"): 0,
        ("P", "P"): 99555,
        ("P", "S"): 36,
        ("P", "T"): 40,
        ("P", "W"): 0,
        ("P", "Y"): 4,
        ("P", "V"): 4,

        ("S", "A"): 360,
        ("S", "R"): 9,
        ("S", "N"): 58,
        ("S", "D"): 0,
        ("S", "C"): 87,
        ("S", "Q"): 13,
        ("S", "E"): 7,
        ("S", "G"): 116,
        ("S", "H"): 0,
        ("S", "I"): 40,
        ("S", "L"): 78,
        ("S", "K"): 2,
        ("S", "M"): 2,
        ("S", "F"): 58,
        ("S", "P"): 16,
        ("S", "S"): 98844,
        ("S", "T"): 244,
        ("S", "W"): 2,
        ("S", "Y"): 40,
        ("S", "V"): 24,

        ("T", "A"): 399,
        ("T", "R"): 10,
        ("T", "N"): 38,
        ("T", "D"): 12,
        ("T", "C"): 26,
        ("T", "Q"): 4,
        ("T", "E"): 4,
        ("T", "G"): 24,
        ("T", "H"): 8,
        ("T", "I"): 296,
        ("T", "L"): 51,
        ("T", "K"): 4,
        ("T", "M"): 63,
        ("T", "F"): 18,
        ("T", "P"): 20,
        ("T", "S"): 265,
        ("T", "T"): 98657,
        ("T", "W"): 2,
        ("T", "Y"): 6,
        ("T", "V"): 95,

        ("W", "A"): 0,
        ("W", "R"): 129,
        ("W", "N"): 5,
        ("W", "D"): 0,
        ("W", "C"): 37,
        ("W", "Q"): 0,
        ("W", "E"): 0,
        ("W", "G"): 23,
        ("W", "H"): 0,
        ("W", "I"): 5,
        ("W", "L"): 92,
        ("W", "K"): 0,
        ("W", "M"): 5,
        ("W", "F"): 9,
        ("W", "P"): 0,
        ("W", "S"): 5,
        ("W", "T"): 5,
        ("W", "W"): 99593,
        ("W", "Y"): 9,
        ("W", "V"): 83,

        ("Y", "A"): 3,
        ("Y", "R"): 0,
        ("Y", "N"): 3,
        ("Y", "D"): 3,
        ("Y", "C"): 73,
        ("Y", "Q"): 0,
        ("Y", "E"): 0,
        ("Y", "G"): 0,
        ("Y", "H"): 92,
        ("Y", "I"): 13,
        ("Y", "L"): 19,
        ("Y", "K"): 16,
        ("Y", "M"): 16,
        ("Y", "F"): 172,
        ("Y", "P"): 3,
        ("Y", "S"): 70,
        ("Y", "T"): 10,
        ("Y", "W"): 6,
        ("Y", "Y"): 99493,
        ("Y", "V"): 6,

        ("V", "A"): 252,
        ("V", "R"): 0,
        ("V", "N"): 2,
        ("V", "D"): 3,
        ("V", "C"): 41,
        ("V", "Q"): 0,
        ("V", "E"): 6,
        ("V", "G"): 46,
        ("V", "H"): 2,
        ("V", "I"): 763,
        ("V", "L"): 220,
        ("V", "K"): 1,
        ("V", "M"): 77,
        ("V", "F"): 32,
        ("V", "P"): 1,
        ("V", "S"): 11,
        ("V", "T"): 41,
        ("V", "W"): 16,
        ("V", "Y"): 2,
        ("V", "V"): 98485,
    }

def run():

    advanced_heatmaps()

    basic_num()

    heatmap_run()
