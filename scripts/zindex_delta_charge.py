from __future__ import division

import matplotlib
import numpy as np
import pytz
from django.db import models
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

from scripts.populate_general_functions import *

# env LDFLAGS="-I/usr/local/opt/openssl/include -L/usr/local/opt/openssl/lib" pip install_vars psycopg2
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import random

evidence_source="UniProt"
# The original kyte & doolittle scale

#kyte = {'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5, 'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5, 'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6, 'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2}

#charge
charge = {'A': 0, 'R': 1, 'N': 0, 'D': -1, 'C': 0, 'Q': 0, 'E': -1, 'G': 0, 'H': 0, 'I': 0, 'L': 0, 'K': 1, 'M': 0, 'F': 0, 'P': 0, 'S': 0, 'T': 0, 'W': 0, 'Y': 0, 'V': 0}
def plot_violins(x, y, this_color):
    print("Plotting bubbles")
    biggest_bubble=0



    for x_position, points in enumerate(y):
        for charge_difference in set(points):
            if charge_difference == 0:
                pass
            else:
                y_value=points.count(charge_difference)

                area = (y_value)**2

                plt.scatter(x[x_position], charge_difference, s=area,c=this_color, alpha=0.5)
    # plt.show()



def get_variant_z_disease_tmh():
#print(Variant.objects.filter(residue__protein__uniprot_id=aline).values_list("residue__sequence_position"))
    vars=list(Variant.objects.filter(disease_status='d', variant_source="ClinVar", residue__tmh_residue__tmh_id__meta_tmh=True).distinct('pk').values_list("residue__tmh_residue__amino_acid_location_in_to_out", "aa_wt", "aa_mut"))
    print(len(vars), "disease tmh vars")
    return(vars)

def get_variant_z_disease_flank():
#print(Variant.objects.filter(residue__protein__uniprot_id=aline).values_list("residue__sequence_position"))

    vars=list(Variant.objects.filter(disease_status='d', variant_source="ClinVar", residue__flank_residue__flank__tmh__meta_tmh=True).distinct('pk').values_list("residue__flank_residue__amino_acid_location_in_to_out", "aa_wt", "aa_mut"))
    print(len(vars), "disease flank vars")
    return(vars)


def get_variant_z_benign_tmh():
#print(Variant.objects.filter(residue__protein__uniprot_id=aline).values_list("residue__sequence_position"))
    vars=list(Variant.objects.filter( variant_source="gnomAD3", residue__tmh_residue__tmh_id__meta_tmh=True).exclude(aa_mut=F("aa_wt")).distinct('pk').values_list("residue__tmh_residue__amino_acid_location_in_to_out", "aa_wt", "aa_mut"))
    print(len(vars), "benign tmh vars")
    return(vars)
def get_variant_z_benign_flank():
#print(Variant.objects.filter(residue__protein__uniprot_id=aline).values_list("residue__sequence_position"))

    vars=list(Variant.objects.filter(variant_source="gnomAD3", residue__flank_residue__flank__tmh__meta_tmh=True).exclude(aa_mut=F("aa_wt")).distinct('pk').values_list("residue__flank_residue__amino_acid_location_in_to_out", "aa_wt", "aa_mut"))
    print(len(vars), "benign flank vars")

    return(vars)

def delta_hydro(aa_wildtype, aa_varianttype):
    delta=abs(charge[aa_wildtype]-charge[aa_varianttype])
    if charge[aa_wildtype] > charge[aa_varianttype]:
        delta=0-delta

    return(delta)

def hydro_list(list):
    new_list=[]
    for i in list:
        postion=i[0]
        aa_wildtype=i[1]
        aa_varianttype=i[2]
        delta_hydrophobicity = delta_hydro(aa_wildtype, aa_varianttype)
        new_list.append([postion, delta_hydrophobicity])
    return(new_list)

def sorting_list(a_list):
    x=[]
    for items in a_list:
        x.append(int(items[0]))

    x=sorted(list(set(x)))
    filled_in_x=[]
    for each_x in range(min(x),max(x)):
        filled_in_x.append(each_x)
    x=filled_in_x
    y=[]

    for each_x in x:
        xy=[]
        isx=False
        for items in a_list:

            if items[0]==each_x:

                #print(items[1])
                xy.append(float(items[1]))
                isx=True
        #print(each_x, xy)
        if isx ==False:
            #print(each_x)
            y.append([0])
        else:
            y.append(xy)
    return(x,y)

def plot_sort_x_y(a_list, a_color):
    #build_x_positions:
    x_and_y=sorting_list(a_list)
    x=x_and_y[0]
    y=x_and_y[1]
    plot_violins(x,y, a_color)

def disease_prop_normalise(disease_list, benign_list):
    dx_and_y=sorting_list(disease_list)
    dx=dx_and_y[0]
    dy=dx_and_y[1]

    bx_and_y=sorting_list(benign_list)
    bx=bx_and_y[0]
    by=bx_and_y[1]
    #print(dy)
    ### Charge = 0
    benign_value=0
    for x_position_in_list, dy_values in enumerate(dy):
        x_position = dx[x_position_in_list]
        for number, an_x in enumerate(bx):
            #print(x_position, an_x)
            if an_x ==x_position:
                benign_value=by[number].count(0)

        disease_value=dy_values.count(0)

        if disease_value>0 and benign_value>0:
            normalise_value=disease_value/benign_value
            #print(x_position, charge_change, area)
            print(x_position,normalise_value)
            plt.scatter(x_position,normalise_value , c="gray", alpha=0.5)

    benign_value=0

    ### Charge =1
    for x_position_in_list, dy_values in enumerate(dy):
        x_position = dx[x_position_in_list]

        for number, an_x in enumerate(bx):
            #print(x_position, an_x)
            if an_x ==x_position:
                if number <= len(by[number]):
                    # How many are at the position minus how many at that position are 0
                    benign_value=by[number].count(1)

        # How many are at the position minus how many at that position are 0
        disease_value=dy_values.count(1)

        if disease_value>0 and benign_value>0:
            normalise_value=disease_value/benign_value
            #print(x_position, charge_change, area)
            print(x_position,normalise_value)
            plt.scatter(x_position,normalise_value , c="plum", alpha=0.5)

        ### Charge =2
        #for x_position_in_list, dy_values in enumerate(dy):
        #    x_position = dx[x_position_in_list]
#
        #    for number, an_x in enumerate(bx):
        #        print(x_position, an_x)
        #        if an_x ==x_position:
        #            if number <= len(by[number]):
        #                # How many are at the position minus how many at that #position are 0
        #                benign_value=by[number].count(2)
#
        #    # How many are at the position minus how many at that position are 0
        #    disease_value=dy_values.count(2)
#
        #    if disease_value>0 and benign_value>0:
        #        normalise_value=disease_value/benign_value
        #        #print(x_position, charge_change, area)
        #        plt.scatter(x_position,normalise_value , c="indigo", alpha=0.5)


def clean_positions(position_list):
    clean_list=[]
    for entry in position_list:
        #print(list(entry))
        if entry[0] == None:
            pass
        else:
            clean_list.append(list(entry))

    return(clean_list)

def run():


    list = get_variant_z_disease_tmh()
    clean_list= clean_positions(list)
    flank_list = get_variant_z_disease_flank()
    flank_clean_list= clean_positions(flank_list)
    clean_list=clean_list+flank_clean_list
    disease_hydrophobicity_list = hydro_list(clean_list)
    plot_sort_x_y(disease_hydrophobicity_list, "firebrick")

    name="delta_charge_z_disease"
    plt.title("Charge change across the TMH")
    plt.xlabel("Distance in residues from TMH center")
    plt.ylabel("ΔQ")
    plt.xlim(-20, 20)
    plt.ylim(-3, 3)
    file_name = str(name + ".pdf")
    plt.savefig(file_name, dpi=700)
    file_name = str(name + ".png")
    plt.savefig(file_name, dpi=700)
    plt.close()

    list = get_variant_z_benign_tmh()
    clean_list= clean_positions(list)
    flank_list = get_variant_z_benign_flank()
    flank_clean_list= clean_positions(flank_list)
    clean_list=clean_list+flank_clean_list
    benign_hydrophobicity_list = hydro_list(clean_list)
    plot_sort_x_y(benign_hydrophobicity_list, "steelblue")

    # Save bitmap and vector images of figure, then flush figure.
    name="delta_charge_z_benign"
    plt.title("Charge change across the TMH")
    plt.xlabel("Distance in residues from TMH center")
    plt.ylabel("ΔQ")
    plt.xlim(-20, 20)
    plt.ylim(-3, 3)
    file_name = str(name + ".pdf")
    plt.savefig(file_name, dpi=700)
    file_name = str(name + ".png")
    plt.savefig(file_name, dpi=700)
    plt.close()

    disease_prop_normalise(disease_hydrophobicity_list, benign_hydrophobicity_list)
    name="delta_charge_z_propensity"
    plt.title("Disease propensity of charge change across the TMH")
    plt.xlabel("Distance in residues from TMH center")
    plt.ylabel("z disease propensity")
    plt.xlim(-20, 20)
    #plt.ylim(-3, 3)
    file_name = str(name + ".pdf")
    plt.savefig(file_name, dpi=700)
    file_name = str(name + ".png")
    plt.savefig(file_name, dpi=700)
    plt.close()
