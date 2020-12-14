from __future__ import division

import matplotlib
import numpy as np
import pytz
from django.db import models

from scripts.populate_general_functions import *
from tmh_db.models import Database_Metadata
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

# env LDFLAGS="-I/usr/local/opt/openssl/include -L/usr/local/opt/openssl/lib" pip install_vars psycopg2
matplotlib.use("TkAgg")
import matplotlib.pyplot as plt
import random

count_threshold = 50
evidence_source = "UniProt"

kyte = {
    "A": 1.8,
    "R": -4.5,
    "N": -3.5,
    "D": -3.5,
    "C": 2.5,
    "Q": -3.5,
    "E": -3.5,
    "G": -0.4,
    "H": -3.2,
    "I": 4.5,
    "L": 3.8,
    "K": -3.9,
    "M": 1.9,
    "F": 2.8,
    "P": -1.6,
    "S": -0.8,
    "T": -0.7,
    "W": -0.9,
    "Y": -1.3,
    "V": 4.2,
}


aa_type = {
    "A": "Hydrophobic",
    "R": "Charged",
    "N": "Polar",
    "D": "Charged",
    "C": "Special",
    "Q": "Polar",
    "E": "Charged",
    "G": "Special",
    "H": "Polar",
    "I": "Hydrophobic",
    "L": "Hydrophobic",
    "K": "Charged",
    "M": "Hydrophobic",
    "F": "Hydrophobic",
    "P": "Special",
    "S": "Polar",
    "T": "Polar",
    "W": "Hydrophobic",
    "Y": "Hydrophobic",
    "V": "Hydrophobic",
}

aa_charge = {"R": "Positive", "D": "Negative", "E": "Negative", "K": "Positive"}


def delta_hydro(aa_wildtype, aa_varianttype):
    return kyte[aa_wildtype] + kyte[aa_varianttype]


def mutation_type(aa_wildtype, aa_varianttype):
    type = (aa_type[aa_wildtype], aa_type[aa_varianttype])
    return type


def clean_positions(position_list):
    clean_list = []
    for entry in position_list:
        # print(list(entry))
        if entry[0] is None:
            pass
        else:
            clean_list.append(list(entry))

    return clean_list


def plot_positions_disease_propensity(
    all_positions_disease, all_positions_gnomad, colour
):
    count_coordinates = []
    position_coordinates = []

    d_all_positions_z = []
    for i in all_positions_disease:
        d_all_positions_z.append(i[0])

    b_all_positions_z = []
    for n in all_positions_gnomad:
        b_all_positions_z.append(n[0])

    for number in range(min(d_all_positions_z), max(d_all_positions_z)):
        d_count = 0
        b_count = 0
        for z in d_all_positions_z:
            if z == number:
                d_count = d_count + 1

        for z in b_all_positions_z:
            if z == number:
                b_count = b_count + 1

        if b_count >= count_threshold and d_count >= count_threshold:
            print(number, d_count, b_count)
            count_coordinates.append(d_count / b_count)
            position_coordinates.append(number)
        else:
            pass

    # normalised_coordinates=normalise_count(count_coordinates)
    plt.scatter(position_coordinates, count_coordinates, c=colour)

    # plt.show()


def plot_positions_hydro_change(all_positions_disease, all_positions_gnomad, colour):

    d_all_positions_z = []
    for i in all_positions_disease:
        d_all_positions_z.append(i[0])

    b_all_positions_z = []
    for n in all_positions_gnomad:
        b_all_positions_z.append(n[0])

    positions = []
    x_axis = []
    for number in range(min(d_all_positions_z), max(d_all_positions_z)):
        entry = [number, []]
        positions.append(entry)
        x_axis.append(number)

    for variant_info in all_positions_disease:
        print(variant_info)
        coordinate = variant_info[0]
        start_aa = variant_info[1]
        variant_aa = variant_info[2]
        delta_hydro = kyte[start_aa] - kyte[variant_aa]
        for position_number, axis_coordinate in enumerate(positions):
            if axis_coordinate[0] == coordinate:
                positions[position_number][1].append(delta_hydro)

    print(positions)
    y_axis = []
    for position in positions:
        avg_delata_hydro = np.mean(position[1])
        y_axis.append(avg_delata_hydro)

    # normalised_coordinates=normalise_count(count_coordinates)
    plt.scatter(x_axis, y_axis, c=colour)

    # plt.show()


def get_variant_z_disease_tmh():
    # print(Variant.objects.filter(residue__protein__uniprot_id=aline).values_list("residue__sequence_position"))
    return list(
        Variant.objects.filter(
            disease_status="d", residue__tmh_residue__evidence=evidence_source
        )
        .distinct("pk")
        .values_list(
            "residue__tmh_residue__amino_acid_location_in_to_out", "aa_wt", "aa_mut"
        )
    )


def get_variant_z_disease_flank():
    # print(Variant.objects.filter(residue__protein__uniprot_id=aline).values_list("residue__sequence_position"))
    return list(
        Variant.objects.filter(
            disease_status="d", residue__flank_residue__evidence=evidence_source
        )
        .distinct("pk")
        .values_list(
            "residue__flank_residue__amino_acid_location_in_to_out", "aa_wt", "aa_mut"
        )
    )


def get_variant_z_benign_tmh():
    # print(Variant.objects.filter(residue__protein__uniprot_id=aline).values_list("residue__sequence_position"))
    return list(
        Variant.objects.filter(
            disease_status="n", residue__tmh_residue__evidence=evidence_source
        )
        .distinct("pk")
        .values_list(
            "residue__tmh_residue__amino_acid_location_in_to_out", "aa_wt", "aa_mut"
        )
    )


def get_variant_z_benign_flank():
    # print(Variant.objects.filter(residue__protein__uniprot_id=aline).values_list("residue__sequence_position"))
    return list(
        Variant.objects.filter(
            disease_status="n", residue__flank_residue__evidence=evidence_source
        )
        .distinct("pk")
        .values_list(
            "residue__flank_residue__amino_acid_location_in_to_out", "aa_wt", "aa_mut"
        )
    )


def add_to_disease_prop_plot(list_disease, list_benign, colour):
    disease_clean_positions_list = clean_positions(list_disease)
    benign_clean_positions_list = clean_positions(list_benign)
    plt.title("Disease propensity variants across the TMH")
    plt.xlabel("Distance in residues from TMH center")
    plt.ylabel("Frequency normalised to sum of type")
    plot_positions_disease_propensity(
        disease_clean_positions_list, benign_clean_positions_list, colour
    )


def add_to_delta_hydro_plot(list_disease, list_benign, colour):
    disease_clean_positions_list = clean_positions(list_disease)
    benign_clean_positions_list = clean_positions(list_benign)
    plt.title("Hydrophobicity change across the TMH")
    plt.xlabel("Distance in residues from TMH center")
    plt.ylabel("Frequency normalised to sum of type")
    plot_positions_hydro_change(
        disease_clean_positions_list, benign_clean_positions_list, colour
    )

    # print(clean_positions_list)


def run():
    d_tmh_positions = get_variant_z_disease_tmh()
    b_tmh_positions = get_variant_z_benign_tmh()

    d_flank_positions = get_variant_z_disease_flank()
    b_flank_positions = get_variant_z_benign_flank()

    add_to_disease_prop_plot(d_tmh_positions, b_tmh_positions, "firebrick")
    add_to_disease_prop_plot(d_flank_positions, b_flank_positions, "lightcoral")
    plt.savefig("disease_propensity_z.png")

    add_to_delta_hydro_plot(d_tmh_positions, b_tmh_positions, "firebrick")
    add_to_delta_hydro_plot(d_flank_positions, b_flank_positions, "lightcoral")
    plt.savefig("delta_hydro_z.png")

    # plt.title("Disease propensity variants across the TMH")
    # plt.xlabel("Distance in residues from TMH center")
    # plt.ylabel("Frequency normalised to sum of type")

    # plt.savefig('disease_propensity_z.png')

    # plot_positions_hydro_change(disease_clean_positions_list,benign_clean_positions_list, colour)
