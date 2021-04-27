from scripts.populate_general_functions import *
import numpy as np
import matplotlib.pyplot as plt
from datetime import date
from datetime import datetime

# Shell Plus Django Imports

date = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

plt.rcdefaults()


def barchart(var_aa_dictionary, title):
    print(f"{title}")
    aa_x = aa_baezo_order()
    y_pos = np.arange(len(aa_x))
    performance = []
    for i in aa_baezo_order():
        performance.append(var_aa_dictionary[i])
        print(f"{i}, {var_aa_dictionary[i]}")
    plt.bar(y_pos, performance, align="center", alpha=0.5)
    plt.xticks(y_pos, aa_x)
    plt.ylabel(f"Count")
    plt.xlabel(f"Amino acid type")
    plt.title(f"{title}")
    filename = f"images/{title}_barchart_aa_count_{date}.png"
    plt.savefig(filename, bbox_inches="tight")
    plt.clf()


def barchart_dictionary(query_set_of_variants):
    aa_dict = {}
    list_of_variants = list(query_set_of_variants)
    #   print(list_of_variants)
    for aa in aa_baezo_order():
        aa_count = list_of_variants.count(aa)

        aa_dict[aa] = aa_count
    return aa_dict


# Inside flanks
inside_flank_residues = (
    Residue.objects.filter(
        flank_residue__feature_location="Inside flank",
        flank_residue__flank__tmh__meta_tmh=True,
        protein__total_tmh_number=1,
    )
    .distinct("pk")
    .values_list("amino_acid_type", flat=True)
)
barchart(barchart_dictionary(inside_flank_residues), "Singlepass inside flank residues")

# Outside flanks
outside_flank_residues = (
    Residue.objects.filter(
        flank_residue__feature_location="Outside flank",
        flank_residue__flank__tmh__meta_tmh=True,
        protein__total_tmh_number=1,
    )
    .distinct("pk")
    .values_list("amino_acid_type", flat=True)
)
barchart(barchart_dictionary(outside_flank_residues), "Singlepass outside flank residues")

# Inside flanks
inside_flank_residues = (
    Residue.objects.filter(
        flank_residue__feature_location="Inside flank",
        flank_residue__flank__tmh__meta_tmh=True,
        protein__total_tmh_number__gt=1,
    )
    .distinct("pk")
    .values_list("amino_acid_type", flat=True)
)
barchart(barchart_dictionary(inside_flank_residues), "Multipass inside flank residues")

# Outside flanks
outside_flank_residues = (
    Residue.objects.filter(
        flank_residue__feature_location="Outside flank",
        flank_residue__flank__tmh__meta_tmh=True,
        protein__total_tmh_number__gt=1,
    )
    .distinct("pk")
    .values_list("amino_acid_type", flat=True)
)
barchart(barchart_dictionary(outside_flank_residues), "Multipass outside flank residues")


# Single pass

singlepass_residues = (
    Residue.objects.filter(
        protein__total_tmh_number=1, tmh_residue__tmh_id__meta_tmh=True
    )
    .distinct("pk")
    .values_list("amino_acid_type", flat=True)
)
barchart(barchart_dictionary(singlepass_residues), "Singlepass TMH residues")


# multipass

multipass_residues = (
    Residue.objects.filter(
        protein__total_tmh_number__gt=1, tmh_residue__tmh_id__meta_tmh=True
    )
    .distinct("pk")
    .values_list("amino_acid_type", flat=True)
)
barchart(barchart_dictionary(multipass_residues), "Multipass TMH residues")


# Helix
helix_residues = (
    Residue.objects.filter(non_tmh_helix_residue__nont_tmh_helix_id__helix_start__gte=0)
    .exclude(tmh_residue__tmh_id__meta_tmh=True)
    .distinct("pk")
    .values_list("amino_acid_type", flat=True)
)
barchart(barchart_dictionary(helix_residues), "Non-TMH helix residues")


# Pore
pore_residues = (
    Residue.objects.filter(
        tmh_residue__feature_location="TMH",
        structural_residue__pore_residue=True,
        tmh_residue__tmh_id__meta_tmh=True,
    )
    .distinct("pk")
    .values_list("amino_acid_type", flat=True)
)
barchart(barchart_dictionary(pore_residues), "Pore residues")

# Simple
simple_residues = (
    Residue.objects.filter(
        tmh_residue__tmh_id__meta_tmh=True, tmh_residue__tmh_id__tmh_tmsoc__test_score__lte="-3.45"
    )
    .distinct("pk")
    .values_list("amino_acid_type", flat=True)
)
barchart(barchart_dictionary(simple_residues), "Simple TMH residues")

# Complex
complex_residues = (
    Residue.objects.filter(
        tmh_residue__tmh_id__meta_tmh=True,
        tmh_residue__tmh_id__tmh_tmsoc__test_score__gte="0.44"
    )
    .distinct("pk")
    .values_list("amino_acid_type", flat=True)
)
barchart(barchart_dictionary(complex_residues), "Complex TMH residues")

# spont
spont_residues = (
    Residue.objects.filter(
        tmh_residue__tmh_id__meta_tmh=True,
        tmh_residue__tmh_id__tmh_deltag__test_score__lte="-0.018"
    )
    .distinct("pk")
    .values_list("amino_acid_type", flat=True)
)
barchart(barchart_dictionary(spont_residues), "Spontaneous TMH residues")

# Non spont
non_spont_residues = (
    Residue.objects.filter(
        tmh_residue__tmh_id__meta_tmh=True,
        tmh_residue__tmh_id__tmh_deltag__test_score__gte="2.601"
    )
    .distinct("pk")
    .values_list("amino_acid_type", flat=True)
)
barchart(barchart_dictionary(non_spont_residues), "Non-spontaneous TMH residues")


def run():
    print("complete")
