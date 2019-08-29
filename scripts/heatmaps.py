from __future__ import division
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from collections import defaultdict
from tmh_db.models import Binding_residue, Database_Metadata, Flank, Funfam, Funfam_residue, Funfamstatus, Go, Keyword, Pfam, Pfam_residue, Protein, Residue, Structural_residue, Structure, Subcellular_location, Tmh, Tmh_deltag, Tmh_hydrophobicity, Tmh_residue, Tmh_tmsoc, Uniref, Variant
from django.core.cache import cache
from django.conf import settings
from django.contrib.auth import get_user_model
from django.db import transaction
from django.db.models import Avg, Case, Count, F, Max, Min, Prefetch, Q, Sum, When, Exists, OuterRef, Subquery
from django.utils import timezone
from django.urls import reverse
from scripts.graphs import *
import numpy as np
from matplotlib.colors import LogNorm


# This has to be a manual matrix.
# It is copied from a PAM1 mutation matrix taken from 1. Jones DT, Taylor WR,
# Thornton JM. A mutation data matrix for transmembrane proteins.
# FEBS Lett. 1994;339(3):269–75.
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

aa_list_baezo_order = ['K', 'R', 'E', 'D', 'Q', 'H', 'N', 'P',
                       'Y', 'W', 'C', 'M', 'T', 'S', 'G', 'V', 'F', 'A', 'I', 'L']

aa_list_alpha = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K',
                 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

# prefetch the residues with variants
Residue.objects.all().prefetch_related("variant")


def advanced_heatmaps():
    '''
    This should be a function that generates some of the more complicated heatmaps by using simplified substitution matrices rather than the clumsy list of lists I was using

    Step 1 - get a list of [("G", "L"), ..., ("L", "D")]
    Step 2 substitution_dictionary(list)
    Step 3 sub_dict_to_heatmap(substitution_dictionary)
    Step 4 heatmap(sub_dict_to_heatmap)
    '''

    # Thornton sub matrix
    title = "1994 Thornton Substitution Matrix"
    heatmap(sub_dict_to_heatmap(thornton_matrix), title,
            aa_list_baezo_order, "Greens", LogNorm())

    # Disease TMHs
    title = "TMH disease variants"
    disease_tmh_variants = list(Variant.objects.filter(residue__tmh_residue__feature_location="TMH").filter(
        disease_status='d').values_list("aa_wt", "aa_mut"))
    print(title, len(disease_tmh_variants))
    disease_subsitutions_in_tmhs = substitution_dictionary(
        disease_tmh_variants)
    heatmap(sub_dict_to_heatmap(disease_subsitutions_in_tmhs),
            title, aa_list_baezo_order, "coolwarm", None)

    # Disease Non-TMHs
    title = "Non-TMH disease variants"
    disease_non_tmh_variants = list(Variant.objects.filter(
        residue__tmh_residue=None, residue__flank_residue=None).filter(disease_status='d').values_list("aa_wt", "aa_mut"))
    print(title, len(disease_non_tmh_variants))
    disease_subsitutions_not_tmhs = substitution_dictionary(
        disease_non_tmh_variants)
    heatmap(sub_dict_to_heatmap(disease_subsitutions_not_tmhs),
            title, aa_list_baezo_order, "coolwarm", None)

    # non TMHs versus TMHs
    title = "TMH disease count divided by non-TMH disease count"
    relative_tmh_to_non_tmh = subs_normalise_by_dic(
        disease_subsitutions_in_tmhs, disease_subsitutions_not_tmhs)
    heatmap(sub_dict_to_heatmap(relative_tmh_to_non_tmh),
            title, aa_list_baezo_order, "coolwarm", None)

    # gnomad TMHs
    title = "TMH gnomAD variants"
    gnomad_variants = list(Variant.objects.exclude(aa_mut=F("aa_wt")).filter(
        residue__tmh_residue__feature_location="TMH").filter(variant_source='gnomAD').values_list("aa_wt", "aa_mut"))
    print(title, len(gnomad_variants))
    gnomad_subsitutions_in_tmhs = substitution_dictionary(gnomad_variants)
    heatmap(sub_dict_to_heatmap(gnomad_subsitutions_in_tmhs),
            title, aa_list_baezo_order, "coolwarm", None)

    # Disease propensity
    title = "TMH disease propensity"
    disease_propensity_in_tmhs = subs_normalise_by_dic(
        disease_subsitutions_in_tmhs, gnomad_subsitutions_in_tmhs)
    disease_propensity_in_tmhs_lists = sub_dict_to_heatmap(
        disease_propensity_in_tmhs)
    # print(disease_propensity_in_tmhs_lists)
    # var_freqs_list, title, aa_list_baezo_order, color, is_it_log
    heatmap(disease_propensity_in_tmhs_lists, title,
            aa_list_baezo_order, "coolwarm", None)

    # Disease variants in TMHs divided by Previous Thornton Dayhoff-like
    title = "Disease count divided by Dayhoff-like TM matrix (Thornton 1994)"
    evolution_dict = subs_normalise_by_dic(
        disease_subsitutions_in_tmhs, thornton_matrix)
    heatmap(sub_dict_to_heatmap(evolution_dict), title,
            aa_list_baezo_order, "coolwarm", None)

    # Single pass disease propensity
    title = "Singlepass TMH disease variants"
    disease_sp_tmh_variants = list(Variant.objects.filter(residue__tmh_residue__feature_location="TMH").filter(
        disease_status='d').filter(residue__protein__total_tmh_number=1).values_list("aa_wt", "aa_mut"))
    print(title, len(disease_sp_tmh_variants))
    disease_sp_tmh_variants_dict = substitution_dictionary(
        disease_sp_tmh_variants)
    heatmap(sub_dict_to_heatmap(disease_sp_tmh_variants_dict),
            title, aa_list_baezo_order, "coolwarm", None)

    title = "Disease propensity in singlepass"
    disease_sp_tmh_variants = list(Variant.objects.filter(residue__tmh_residue__feature_location="TMH").filter(
        disease_status='d').filter(residue__protein__total_tmh_number=1).values_list("aa_wt", "aa_mut"))
    disease_sp_tmh_variants_dict = substitution_dictionary(
        disease_sp_tmh_variants)
    gnomad_sp_tmh_variants = list(Variant.objects.exclude(aa_mut=F("aa_wt")).filter(residue__tmh_residue__feature_location="TMH").filter(
        variant_source='gnomAD').filter(residue__protein__total_tmh_number=1).values_list("aa_wt", "aa_mut"))
    gnomad_sp_tmh_variants_dict = substitution_dictionary(
        gnomad_sp_tmh_variants)
    print(title, "disease:", len(disease_sp_tmh_variants),
          "gnomAD:", len(gnomad_sp_tmh_variants))
    sp_disease_propensity = subs_normalise_by_dic(
        disease_sp_tmh_variants_dict, gnomad_sp_tmh_variants_dict)
    heatmap(sub_dict_to_heatmap(sp_disease_propensity),
            title, aa_list_baezo_order, "Reds", None)

    # Multipass disease propensity

    title = "Multipass TMH disease variants"
    disease_mp_tmh_variants = list(Variant.objects.filter(residue__tmh_residue__feature_location="TMH").filter(
        disease_status='d').exclude(residue__protein__total_tmh_number=1).values_list("aa_wt", "aa_mut"))
    print(title, len(disease_mp_tmh_variants))
    disease_mp_tmh_variants = substitution_dictionary(disease_mp_tmh_variants)
    heatmap(sub_dict_to_heatmap(disease_mp_tmh_variants),
            title, aa_list_baezo_order, "coolwarm", None)

    title = "Disease propensity in multipass"
    print(title)
    disease_mp_tmh_variants = list(Variant.objects.filter(residue__tmh_residue__feature_location="TMH").filter(
        disease_status='d').exclude(residue__protein__total_tmh_number=1).values_list("aa_wt", "aa_mut"))
    disease_mp_tmh_variants_dict = substitution_dictionary(
        disease_mp_tmh_variants)
    gnomad_mp_tmh_variants = list(Variant.objects.exclude(aa_mut=F("aa_wt")).filter(residue__tmh_residue__feature_location="TMH").filter(
        variant_source='gnomAD').exclude(residue__protein__total_tmh_number=1).values_list("aa_wt", "aa_mut"))
    gnomad_mp_tmh_variants_dict = substitution_dictionary(
        gnomad_mp_tmh_variants)
    print(title, "disease:", len(disease_mp_tmh_variants),
          "gnomAD:", len(gnomad_mp_tmh_variants))
    mp_disease_propensity = subs_normalise_by_dic(
        disease_mp_tmh_variants_dict, gnomad_mp_tmh_variants_dict)
    heatmap(sub_dict_to_heatmap(mp_disease_propensity),
            title, aa_list_baezo_order, "Reds", None)

    # Inside flanks
    title = "Disease variants in inside flank"
    disease_inside_flank_variants = list(Variant.objects.filter(
        residue__flank_residue__feature_location="Inside flank").filter(disease_status='d').values_list("aa_wt", "aa_mut"))
    print(title, len(disease_inside_flank_variants))
    disease_subsitutions_inside_flanks = substitution_dictionary(
        disease_inside_flank_variants)
    heatmap(sub_dict_to_heatmap(disease_subsitutions_inside_flanks),
            title, aa_list_baezo_order, "coolwarm", None)

    # Oustide flanks
    title = "Disease variants in outside flank"
    disease_outside_flank_variants = list(Variant.objects.filter(
        residue__flank_residue__feature_location="Outside flank").filter(disease_status='d').values_list("aa_wt", "aa_mut"))
    print(title, len(disease_outside_flank_variants))
    disease_subsitutions_outside_flanks = substitution_dictionary(
        disease_outside_flank_variants)
    heatmap(sub_dict_to_heatmap(disease_subsitutions_outside_flanks),
            title, aa_list_baezo_order, "coolwarm", None)

    # Inside flanks multipass
    title = "Disease variants in inside flank of multipass"
    disease_inside_flank_variants = list(Variant.objects.filter(residue__flank_residue__feature_location="Inside flank").filter(
        disease_status='d').exclude(residue__protein__total_tmh_number=1).values_list("aa_wt", "aa_mut"))
    print(title, len(disease_inside_flank_variants))
    disease_subsitutions_inside_flanks = substitution_dictionary(
        disease_inside_flank_variants)
    heatmap(sub_dict_to_heatmap(disease_subsitutions_inside_flanks),
            title, aa_list_baezo_order, "coolwarm", None)

    title = "Disease propensity in inside flank of multipass"
    gnomad_mp_inside_variants = list(Variant.objects.exclude(aa_mut=F("aa_wt")).filter(residue__flank_residue__feature_location="Inside flank").filter(
        variant_source='gnomAD').exclude(residue__protein__total_tmh_number=1).values_list("aa_wt", "aa_mut"))
    mp_inside_disease_propensity = subs_normalise_by_dic(
        disease_subsitutions_inside_flanks, substitution_dictionary(gnomad_mp_inside_variants))
    heatmap(sub_dict_to_heatmap(mp_inside_disease_propensity),
            title, aa_list_baezo_order, "Reds", None)

    # Oustide flanks multipass
    title = "Disease variants in outside flank of multipass"
    disease_outside_flank_variants = list(Variant.objects.filter(residue__flank_residue__feature_location="Outside flank").filter(
        disease_status='d').exclude(residue__protein__total_tmh_number=1).values_list("aa_wt", "aa_mut"))
    print(title, len(disease_outside_flank_variants))
    disease_subsitutions_outside_flanks = substitution_dictionary(
        disease_outside_flank_variants)
    heatmap(sub_dict_to_heatmap(disease_subsitutions_outside_flanks),
            title, aa_list_baezo_order, "coolwarm", None)

    title = "Disease propensity in outside flank of multipass"
    gnomad_mp_outside_variants = list(Variant.objects.exclude(aa_mut=F("aa_wt")).filter(residue__flank_residue__feature_location="Outside flank").filter(
        variant_source='gnomAD').exclude(residue__protein__total_tmh_number=1).values_list("aa_wt", "aa_mut"))
    mp_outside_disease_propensity = subs_normalise_by_dic(
        disease_subsitutions_outside_flanks, substitution_dictionary(gnomad_mp_outside_variants))
    heatmap(sub_dict_to_heatmap(mp_outside_disease_propensity),
            title, aa_list_baezo_order, "Reds", None)

    # Inside flanks singlepass
    title = "Disease variants in inside flank of singlepass"
    disease_inside_flank_variants = list(Variant.objects.filter(residue__flank_residue__feature_location="Inside flank").filter(
        disease_status='d').filter(residue__protein__total_tmh_number=1).values_list("aa_wt", "aa_mut"))
    print(title, len(disease_inside_flank_variants))
    disease_subsitutions_inside_flanks = substitution_dictionary(
        disease_inside_flank_variants)
    heatmap(sub_dict_to_heatmap(disease_subsitutions_inside_flanks),
            title, aa_list_baezo_order, "coolwarm", None)

    title = "Disease propensity in inside flank of singlepass"
    gnomad_sp_inside_variants = list(Variant.objects.exclude(aa_mut=F("aa_wt")).filter(residue__flank_residue__feature_location="Inside flank").filter(
        variant_source='gnomAD').filter(residue__protein__total_tmh_number=1).values_list("aa_wt", "aa_mut"))
    sp_inside_disease_propensity = subs_normalise_by_dic(
        disease_subsitutions_inside_flanks, substitution_dictionary(gnomad_sp_inside_variants))
    heatmap(sub_dict_to_heatmap(sp_inside_disease_propensity),
            title, aa_list_baezo_order, "Reds", None)

    # Oustide flanks singlepass
    title = "Disease variants in outside flank of singlepass"
    disease_outside_flank_variants = list(Variant.objects.filter(residue__flank_residue__feature_location="Outside flank").filter(
        disease_status='d').filter(residue__protein__total_tmh_number=1).values_list("aa_wt", "aa_mut"))
    print(title, len(disease_outside_flank_variants))
    disease_subsitutions_outside_flanks = substitution_dictionary(
        disease_outside_flank_variants)
    heatmap(sub_dict_to_heatmap(disease_subsitutions_outside_flanks),
            title, aa_list_baezo_order, "coolwarm", None)


    title = "Disease propensity in outside flank of singlepass"
    gnomad_sp_outside_variants = list(Variant.objects.exclude(aa_mut=F("aa_wt")).filter(residue__flank_residue__feature_location="Inside flank").filter(
        variant_source='gnomAD').filter(residue__protein__total_tmh_number=1).values_list("aa_wt", "aa_mut"))
    sp_outside_disease_propensity = subs_normalise_by_dic(
        disease_subsitutions_outside_flanks, substitution_dictionary(gnomad_sp_outside_variants))
    heatmap(sub_dict_to_heatmap(sp_outside_disease_propensity),
            title, aa_list_baezo_order, "Reds", None)

    # Multipass and singlepass by tthe thornton matrix
    title = "Singlepass disease count divided by Dayhoff-like TM matrix (Thornton 1994)"
    evolution_dict = subs_normalise_by_dic(
        disease_sp_tmh_variants_dict, thornton_matrix)
    heatmap(sub_dict_to_heatmap(evolution_dict), title,
            aa_list_baezo_order, "coolwarm", None)

    title = "Multipass disease count divided by Dayhoff-like TM matrix (Thornton 1994)"
    evolution_dict = subs_normalise_by_dic(
        disease_mp_tmh_variants_dict, thornton_matrix)
    heatmap(sub_dict_to_heatmap(evolution_dict), title,
            aa_list_baezo_order, "coolwarm", None)


def amino_acid_count_to_dictionary():
    '''
    This takes a list of amino acids and counts how many of each there are.
    A dictionary of amino acids and their count are returned.
    '''

    aa_count={}
    for reference_amino_acid in aa_list_alpha:
        instances = variant_list.count(reference_amino_acid)
        aa_count[reference_amino_acid]=instances
    return(aa_count)


def frequency_wt_normalisation(sub_dictionary, aa_dictionary):
    '''
    This takes a mutation dictionary and divides each value by the frequency of amino acids in a second dictionary.
    It returns a heatmap array of the values.
    '''
    heatmap_array = []
    for row, mutant_amino_acid in enumerate(aa_list_baezo_order):
        heatmap_array.append([])
        for column, reference_amino_acid in enumerate(aa_list_baezo_order):
                heatmap_array[row].append(sub_dictionary[(reference_amino_acid, mutant_amino_acid)]/aa_dictionary[reference_amino_acid])
    return(heatmap_array)


def substitution_dictionary(variant_list):
    '''
    Returns a substiution dictionary in the format {(X,Y):1, (Y,Z):3}

    variant_list
        a list of tuples from a django query in the format [(X, Y), (X, X), ... ,(Y, Z)]
    '''

    subs = {}
    for reference_amino_acid in aa_list_alpha:
        for mutant_amino_acid in aa_list_alpha:
            # print(variant_list)
            instances = variant_list.count(
                (reference_amino_acid, mutant_amino_acid))
            subs[(reference_amino_acid, mutant_amino_acid)] = instances
    # print(subs)
    return(subs)


def subs_normalise_by_dic(first_dic, second_dic):
    '''
    Returns a new normalised substiution dictionary in the format {(X,Y):5, (Y,Z):10}.
    The first dictionary counts are divided by the second.

    first_dic
        a dictionary of amino acid tuples and counts {(X,Y):10, (Y,Z):50}
    second_dic
        a dictionary of amino acid tuples and counts {(X,Y):2, (Y,Z):5}
    '''

    subs = {}
    for reference_amino_acid in aa_list_alpha:
        for mutant_amino_acid in aa_list_alpha:
            try:
                if first_dic[(reference_amino_acid, mutant_amino_acid)] > 0:
                    normalised_value = first_dic[(
                        reference_amino_acid, mutant_amino_acid)] / second_dic[(reference_amino_acid, mutant_amino_acid)]
                else:
                    # In the event of both datasets containing 0, the normalised output should also be 0, not 1.
                    normalised_value = 0
            # This stops the script failing if there is a 0 in the second dictionary
            except ZeroDivisionError:
                normalised_value = 0

            subs[(reference_amino_acid, mutant_amino_acid)] = normalised_value
    # print(subs)
    return(subs)


def sub_dict_to_heatmap(sub_dictionary):
    '''
    Builds a heatmap array (list of lists) ordered by baeza delgado 2013 paper:
    'K', 'R', 'E', 'D', 'Q', 'H', 'N', 'P','Y', 'W', 'C', 'M', 'T', 'S', 'G', 'V', 'F', 'A', 'I', 'L'
    '''
    # I like this order of amino acids. It is a sort of functional family tree of the amino acids that generally decreases in polarity

    heatmap_array = []
    for row, mutant_amino_acid in enumerate(aa_list_baezo_order):
        heatmap_array.append([])
        for column, reference_amino_acid in enumerate(aa_list_baezo_order):
            heatmap_array[row].append(
                sub_dictionary[(reference_amino_acid, mutant_amino_acid)])
    return(heatmap_array)


def run():
    advanced_heatmaps()
