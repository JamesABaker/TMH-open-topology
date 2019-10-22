from __future__ import division
# env LDFLAGS="-I/usr/local/opt/openssl/include -L/usr/local/opt/openssl/lib" pip install_vars psycopg2
from django.db import models
from tmh_db.models import Database_Metadata, Subcellular_location, Uniref, Go, Structure, Structural_residue, Funfam_residue, Funfamstatus, Protein, Residue, Tmh, Tmh_deltag, Tmh_hydrophobicity, Tmh_residue, Tmh_tmsoc, Variant, Keyword, Binding_residue
from django.db.models import Avg, Case, Count, F, Max, Min, Prefetch, Q, Sum, When, Exists, OuterRef, Subquery
import pytz
import numpy as np
from scripts.populate_general_functions import *
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import random

evidence_source="UniProt"
# The original kyte & doolittle scale

kyte = {'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5, 'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5, 'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6, 'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2}

#charge
#kyte = {'A': 0, 'R': 1, 'N': 0, 'D': -1, 'C': 0, 'Q': 0, 'E': -1, 'G': 0, 'H': 0, 'I': 0, 'L': 0, 'K': 1, 'M': 0, 'F': 0, 'P': 0, 'S': 0, 'T': 0, 'W': 0, 'Y': 0, 'V': 0}

def plot_violins(x, y, this_color):
    print("Plotting violins")
    fig, axes = plt.subplots(nrows=1, ncols=1, figsize=(20, 4))

    parts = axes.violinplot(y, x, showmeans=True,  showmedians=True, showextrema=False, points=1000)
    y_size=[]
    for n in y:
        y_size.append(len(n))
    # axes[0].set_title('violin plot')
    # Colour palette is colour blind friendly according to Wong B 2011 Points
    # of view: Color blindness Nature Methods 8:441.
    for x_position, pc in enumerate(parts['bodies']):
        pc.set_facecolor(this_color)
        size=len(y[x_position])/max(y_size)
        print(size)
        pc.set_alpha(size)
        #pc.set_edgecolor('black')
        pc.set_linewidth(0.2)

    parts['cmedians'].set_color('black')
    parts['cmeans'].set_color('#e69f00')


    # plt.show()



def get_variant_z_disease_tmh():
#print(Variant.objects.filter(residue__protein__uniprot_id=aline).values_list("residue__sequence_position"))
    return(list(Variant.objects.filter(disease_status='d', residue__tmh_residue__evidence=evidence_source).distinct('pk').values_list("residue__tmh_residue__amino_acid_location_in_to_out", "aa_wt", "aa_mut")))

def get_variant_z_disease_flank():
#print(Variant.objects.filter(residue__protein__uniprot_id=aline).values_list("residue__sequence_position"))
    return(list(Variant.objects.filter(disease_status='d', residue__flank_residue__evidence=evidence_source).distinct('pk').values_list("residue__flank_residue__amino_acid_location_in_to_out", "aa_wt", "aa_mut")))

def get_variant_z_benign_tmh():
#print(Variant.objects.filter(residue__protein__uniprot_id=aline).values_list("residue__sequence_position"))
    return(list(Variant.objects.filter(disease_status='n', residue__tmh_residue__evidence=evidence_source).exclude(aa_mut=F("aa_wt")).distinct('pk').values_list("residue__tmh_residue__amino_acid_location_in_to_out", "aa_wt", "aa_mut")))

def get_variant_z_benign_flank():
#print(Variant.objects.filter(residue__protein__uniprot_id=aline).values_list("residue__sequence_position"))
    return(list(Variant.objects.filter(disease_status='n', residue__flank_residue__evidence=evidence_source).exclude(aa_mut=F("aa_wt")).distinct('pk').values_list("residue__flank_residue__amino_acid_location_in_to_out", "aa_wt", "aa_mut")))

def delta_hydro(aa_wildtype, aa_varianttype):
    delta=abs(kyte[aa_wildtype]-kyte[aa_varianttype])
    #if aa_wildtype > aa_varianttype:
    #    delta=0-delta

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

def sort_x_y(a_list, a_color):
    #build_x_positions:
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
        print(each_x, xy)
        if isx ==False:
            print(each_x)
            y.append([0])
        else:
            y.append(xy)




    plot_violins(x,y, a_color)



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
    hydrophobicity_list = hydro_list(clean_list)
    sort_x_y(hydrophobicity_list, "firebrick")

    name="delta_hydro_z_disease"
    plt.title("Hydrophobicity change across the TMH")
    plt.xlabel("Distance in residues from TMH center")
    plt.ylabel("ΔH")
    plt.xlim(-20, 20)
    plt.ylim(0, 10)
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
    hydrophobicity_list = hydro_list(clean_list)
    sort_x_y(hydrophobicity_list, "steelblue")




    # Save bitmap and vector images of figure, then flush figure.
    name="delta_hydro_z_benign"
    plt.title("Hydrophobicity change across the TMH")
    plt.xlabel("Distance in residues from TMH center")
    plt.ylabel("ΔH")
    plt.xlim(-20, 20)
    plt.ylim(0, 10)
    file_name = str(name + ".pdf")
    plt.savefig(file_name, dpi=700)
    file_name = str(name + ".png")
    plt.savefig(file_name, dpi=700)
    plt.close()
