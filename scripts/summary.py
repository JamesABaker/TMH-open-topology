# Shell Plus Model Imports
from django.contrib.admin.models import LogEntry
from django.contrib.auth.models import Group, Permission, User
from django.contrib.contenttypes.models import ContentType
from django.contrib.sessions.models import Session
from tmh_db.models import Binding_residue, Database_Metadata, Funfam, Funfam_residue, Funfamstatus, Go, Keyword, Pfam, Pfam_residue, Protein, Residue, Structural_residue, Structure, Subcellular_location, Tmh, Tmh_deltag, Tmh_hydrophobicity, Tmh_residue, Tmh_tmsoc, Uniref, Variant
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
import collections

def heatmap_array(var_freq_dict, aa_order):
    var_freq = collections.Counter(var_freq_dict)
    large_array=[]
    for aa_mut in aa_order:
        aa_array=[]
        for aa_wt in aa_order:
            #This query is counter intuitive. The aa_wt is first in the tuple, the aa_mut is second. The aa_mut is first in the loop to make sure it is on the y axis.
            aa_array.append(var_freq[(aa_wt, aa_mut)])
        large_array.append(aa_array)
    #print(np.array(large_array))
    return(np.array(large_array))


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
    heatmap(var_freqs_list, title, aa_list_baezo_order, "d", None)

    title = "TMH±5_disease_normalised_wt_aa_tmh_freq"
    heatmap(normalise_tmh_resid_array(var_freqs_list, aa_list_baezo_order), title, aa_list_baezo_order, "d", None)


    title = "non_TMH±5_disease"
    var_freq=list(Variant.objects.filter(residue__tmh_residue=None).filter(disease_status='d').values_list("aa_wt", "aa_mut"))
    heatmap(heatmap_array(var_freq, aa_list_baezo_order), title, aa_list_baezo_order, "d", None)

    title = "TMH_disease_clinvar"
    clinvar_tmh_var_freq=list(Variant.objects.exclude(residue__tmh_residue__feature_location="TMH").filter(variant_source="ClinVar", disease_status='d').values_list("aa_wt", "aa_mut"))
    clinvar_tmh_disease_var_array = heatmap_array(clinvar_tmh_var_freq, aa_list_baezo_order)
    heatmap(clinvar_tmh_disease_var_array, title, aa_list_baezo_order , "d", None)

    title = "non_TMH±5_disease_clinvar"
    clinvar_non_tmh_var_freq=list(Variant.objects.filter(residue__tmh_residue=None).filter(variant_source="ClinVar", disease_status='d').values_list("aa_wt", "aa_mut"))
    clinvar_non_tmh_disease_var_array = heatmap_array(clinvar_non_tmh_var_freq, aa_list_baezo_order)
    heatmap(clinvar_non_tmh_disease_var_array, title, aa_list_baezo_order, "d", None)

    title = "TMH±5_disease_humsavar"
    var_freq=list(Variant.objects.exclude(residue__tmh_residue=None).filter(variant_source="Humsavar", disease_status='d').values_list("aa_wt", "aa_mut"))
    heatmap(heatmap_array(var_freq, aa_list_baezo_order), title, aa_list_baezo_order, "d", None)

    title = "non_TMH±5_disease_humsavar"
    var_freq=list(Variant.objects.filter(residue__tmh_residue=None).filter(variant_source="Humsavar", disease_status='d').values_list("aa_wt", "aa_mut"))
    heatmap(heatmap_array(var_freq, aa_list_baezo_order), title, aa_list_baezo_order, "d", None)

    title = "TMH_gnomAD"
    gnomad_tmh_var_freq=list(Variant.objects.exclude(aa_mut=F("aa_wt")).filter(residue__tmh_residue__feature_location="TMH").filter(variant_source='gnomAD').values_list("aa_wt", "aa_mut"))
    gnomad_tmh_var_array = heatmap_array(gnomad_tmh_var_freq, aa_list_baezo_order)
    heatmap(gnomad_tmh_var_array, title, aa_list_baezo_order, "n", None)

    title = "non_TMH±5_gnomAD"
    gnomad_non_tmh_var_freq=list(Variant.objects.exclude(aa_mut=F("aa_wt")).filter(residue__tmh_residue=None).filter(variant_source='gnomAD').values_list("aa_wt", "aa_mut"))
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
    "ClinVar": "d"
    }

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

    non_tmh_residue_num = Residue.objects.filter(tmh_residue=None).count()
    print("Non-TMH residues,", non_tmh_residue_num)
    structural_residue_num = Residue.objects.exclude(structural_residue=None).count()
    print("Residues with a map to at least one residue in a structure,", structural_residue_num)
    non_structural_residue_num = Residue.objects.filter(structural_residue=None).count()
    print("Residues with no map to a residue in a structure,",non_structural_residue_num)
    tmh_structural_residue_num = Residue.objects.exclude(tmh_residue=None).exclude(structural_residue=None).count()
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
    tmh_g_variants_num = Variant.objects.exclude(residue__tmh_residue=None).filter(variant_source='gnomAD').count()
    print("TMH gnomAD variants,", tmh_g_variants_num)
    non_tmh_g_variants_num = Variant.objects.filter(residue__tmh_residue=None).filter(variant_source='gnomAD').count()
    print("Non-TMH gnomAD variants,", non_tmh_g_variants_num)

    d_variants_num = Variant.objects.filter(disease_status='d').count()
    print("Disease variants,", d_variants_num)
    tmh_d_variants_num = Variant.objects.exclude(residue__tmh_residue=None).filter(disease_status='d').count()
    print("TMH disease variants,", tmh_d_variants_num)
    d_variants_clinvar_num = Variant.objects.filter(disease_status='d', variant_source="ClinVar").count()
    print("ClinVar disease variants,", d_variants_clinvar_num)
    tmh_d_variants_clinvar_num = Variant.objects.exclude(residue__tmh_residue=None).filter(
        disease_status='d', variant_source="ClinVar").count()
    print("TMH ClinVar disease variants,", tmh_d_variants_clinvar_num)

    tmh_d_variants_humsavar_num = Variant.objects.exclude(residue__tmh_residue=None).filter(
        disease_status='d', variant_source="Humsavar").count()

    d_variants_humsavar_num = Variant.objects.filter(disease_status='d', variant_source="Humsavar").count()
    print("Humsavar disease variants,", d_variants_humsavar_num)
    print("TMH Humsavar disease variants,", tmh_d_variants_humsavar_num)
    non_tmh_d_variants_num = Variant.objects.filter(residue__tmh_residue=None).filter(disease_status='d').count()
    print("Non-TMH disease variants,", non_tmh_d_variants_num)
    non_tmh_d_variants_clinvar_num = Variant.objects.filter(residue__tmh_residue=None).filter(disease_status='d', variant_source="ClinVar").count()
    print("Non-TMH ClinVar disease variants,", non_tmh_d_variants_clinvar_num)
    non_tmh_d_variants_humsavar_num = Variant.objects.filter(residue__tmh_residue=None).filter(disease_status='d', variant_source="Humsavar").count()
    print("Non-TMH Humsavar disease variants,",
          non_tmh_d_variants_humsavar_num)



    disease_status = {"d":"disease", "n": "benign"}
    for state in ["d", "n"]:
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
    print("Disease variants per TMH residue,",
          tmh_d_variants_num / tmh_residue_num)
    print("Disease variants per non-TMH residue,",
          non_tmh_d_variants_num / non_tmh_residue_num)

    fisher_oddsratio, fisher_pvalue = stats.fisher_exact([[tmh_d_variants_num, non_tmh_d_variants_num], [tmh_g_variants_num, non_tmh_g_variants_num]], alternative='two-sided')


    print("Fisher test Disease TMH non-TMH versus gnomAD TMH non-TMH,", fisher_pvalue)

    obs = np.array([[tmh_d_variants_num, non_tmh_d_variants_num], [tmh_g_variants_num, non_tmh_g_variants_num]])
    chi_test = stats.chi2_contingency(obs)

    print("Chi 2 test Disease TMH non-TMH versus gnomAD TMH non-TMH,", chi_test)

    objects = ("Residues", "TMH ±5 residues", "Non-TMH residues")
    performance = [d_variants_num / residue_num, tmh_d_variants_num /
                   tmh_residue_num, non_tmh_d_variants_num / non_tmh_residue_num]
    barchart(objects, performance, "TMP_disease_variants", "d", "Residue type", "Variants per residue")

    objects = ("Residues", "TMH ±5 residues", "Non-TMH residues")
    performance = [d_variants_clinvar_num / residue_num, tmh_d_variants_clinvar_num /
                   tmh_residue_num, non_tmh_d_variants_clinvar_num / non_tmh_residue_num]
    barchart(objects, performance, "ClinVar_disease_variants", "d", "Residue type", "Variants per residue")

    objects = ("Residues", "TMH ±5 residues", "Non-TMH residues")
    performance = [d_variants_humsavar_num / residue_num, tmh_d_variants_humsavar_num /
                   tmh_residue_num, non_tmh_d_variants_humsavar_num / non_tmh_residue_num]
    barchart(objects, performance, "Humsavar_disease_variants", "d", "Residue type", "Variants per residue")

    objects = ("Residues", "TMH ±5 residues", "Non-TMH residues")
    performance = [g_variants_num / residue_num, tmh_g_variants_num /
                   tmh_residue_num, non_tmh_g_variants_num / non_tmh_residue_num]
    barchart(objects, performance, "TMP_gnomAD_variants", "n", "Residue type", "Variants per residue")


def run():

    basic_num()

    heatmap_run()
