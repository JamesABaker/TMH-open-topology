# Shell Plus Model Imports
from django.contrib.admin.models import LogEntry
from django.contrib.auth.models import Group, Permission, User
from django.contrib.contenttypes.models import ContentType
from django.contrib.sessions.models import Session
from tmh_db.models import Database_Metadata, Disease, Flank, Flank_residue, Funfam, FunfamResidue, Go, Keyword, Non_tmh_helix, Non_tmh_helix_residue, Protein, Residue, Signal_peptide, Signal_residue, Structural_residue, Structure, SubcellularLocation, Tmh, Tmh_deltag, Tmh_hydrophobicity, Tmh_residue, Tmh_tmsoc, Uniref, Variant
# Shell Plus Django Imports
from django.core.cache import cache
from django.conf import settings
from django.contrib.auth import get_user_model
from django.db import transaction
from django.db.models import Avg, Case, Count, F, Max, Min, Prefetch, Q, Sum, When, Exists, OuterRef, Subquery
from django.utils import timezone
from django.urls import reverse

def csv_parser(csv_item):
    '''
    Takes a line and removes commas from quoted lines which should be cells.
    Returns a line
    '''
    line=csv_item.split('"')
    #print(line)
    for n, i in enumerate(line):
        if '"' in i:
            line[n]=[s.strip(',') for s in line[n]]
    line="".join(line)
    currentline = line.split(",")
    return(currentline)

def physico_change(aa_wt, aa_mut):
    '''
    Takes the two amino acids and does a crude job at working out the change.
    ''' 
    '''
    aa_chem={
    'A': "Hydrophobic",
    'R': "Positive",
    'N': "Polar",
    'D': "Negative",
    'C': "Cysteine",  
    'Q': "Polar",
    'E': "Negative",
    'G': "Glycine",
    'H': "Aromatic",
    'I': "Hydrophobic",
    'L': "Hydrophobic",
    'K': "Charged",
    'M': "Hydrophobic",
    'F': "Aromatic",
    'P': "Proline",
    'S': "Polar",
    'T': "Threonine",
    'W': "Aromatic",
    'Y': "Aromatic",
    'V': "Hydrophobic"}
    if aa_chem[aa_wt] == aa_chem[aa_mut]:
        return("small biochemical change")
    elif aa_chem[aa_wt] != aa_chem[aa_mut]:
        return("large biochemical change")
    '''
    
    blosum62 = {
    ('W', 'F'): 1, ('L', 'R'): -2, ('S', 'P'): -1, ('V', 'T'): 0,
    ('Q', 'Q'): 5, ('N', 'A'): -2, ('Z', 'Y'): -2, ('W', 'R'): -3,
    ('Q', 'A'): -1, ('S', 'D'): 0, ('H', 'H'): 8, ('S', 'H'): -1,
    ('H', 'D'): -1, ('L', 'N'): -3, ('W', 'A'): -3, ('Y', 'M'): -1,
    ('G', 'R'): -2, ('Y', 'I'): -1, ('Y', 'E'): -2, ('B', 'Y'): -3,
    ('Y', 'A'): -2, ('V', 'D'): -3, ('B', 'S'): 0, ('Y', 'Y'): 7,
    ('G', 'N'): 0, ('E', 'C'): -4, ('Y', 'Q'): -1, ('Z', 'Z'): 4,
    ('V', 'A'): 0, ('C', 'C'): 9, ('M', 'R'): -1, ('V', 'E'): -2,
    ('T', 'N'): 0, ('P', 'P'): 7, ('V', 'I'): 3, ('V', 'S'): -2,
    ('Z', 'P'): -1, ('V', 'M'): 1, ('T', 'F'): -2, ('V', 'Q'): -2,
    ('K', 'K'): 5, ('P', 'D'): -1, ('I', 'H'): -3, ('I', 'D'): -3,
    ('T', 'R'): -1, ('P', 'L'): -3, ('K', 'G'): -2, ('M', 'N'): -2,
    ('P', 'H'): -2, ('F', 'Q'): -3, ('Z', 'G'): -2, ('X', 'L'): -1,
    ('T', 'M'): -1, ('Z', 'C'): -3, ('X', 'H'): -1, ('D', 'R'): -2,
    ('B', 'W'): -4, ('X', 'D'): -1, ('Z', 'K'): 1, ('F', 'A'): -2,
    ('Z', 'W'): -3, ('F', 'E'): -3, ('D', 'N'): 1, ('B', 'K'): 0,
    ('X', 'X'): -1, ('F', 'I'): 0, ('B', 'G'): -1, ('X', 'T'): 0,
    ('F', 'M'): 0, ('B', 'C'): -3, ('Z', 'I'): -3, ('Z', 'V'): -2,
    ('S', 'S'): 4, ('L', 'Q'): -2, ('W', 'E'): -3, ('Q', 'R'): 1,
    ('N', 'N'): 6, ('W', 'M'): -1, ('Q', 'C'): -3, ('W', 'I'): -3,
    ('S', 'C'): -1, ('L', 'A'): -1, ('S', 'G'): 0, ('L', 'E'): -3,
    ('W', 'Q'): -2, ('H', 'G'): -2, ('S', 'K'): 0, ('Q', 'N'): 0,
    ('N', 'R'): 0, ('H', 'C'): -3, ('Y', 'N'): -2, ('G', 'Q'): -2,
    ('Y', 'F'): 3, ('C', 'A'): 0, ('V', 'L'): 1, ('G', 'E'): -2,
    ('G', 'A'): 0, ('K', 'R'): 2, ('E', 'D'): 2, ('Y', 'R'): -2,
    ('M', 'Q'): 0, ('T', 'I'): -1, ('C', 'D'): -3, ('V', 'F'): -1,
    ('T', 'A'): 0, ('T', 'P'): -1, ('B', 'P'): -2, ('T', 'E'): -1,
    ('V', 'N'): -3, ('P', 'G'): -2, ('M', 'A'): -1, ('K', 'H'): -1,
    ('V', 'R'): -3, ('P', 'C'): -3, ('M', 'E'): -2, ('K', 'L'): -2,
    ('V', 'V'): 4, ('M', 'I'): 1, ('T', 'Q'): -1, ('I', 'G'): -4,
    ('P', 'K'): -1, ('M', 'M'): 5, ('K', 'D'): -1, ('I', 'C'): -1,
    ('Z', 'D'): 1, ('F', 'R'): -3, ('X', 'K'): -1, ('Q', 'D'): 0,
    ('X', 'G'): -1, ('Z', 'L'): -3, ('X', 'C'): -2, ('Z', 'H'): 0,
    ('B', 'L'): -4, ('B', 'H'): 0, ('F', 'F'): 6, ('X', 'W'): -2,
    ('B', 'D'): 4, ('D', 'A'): -2, ('S', 'L'): -2, ('X', 'S'): 0,
    ('F', 'N'): -3, ('S', 'R'): -1, ('W', 'D'): -4, ('V', 'Y'): -1,
    ('W', 'L'): -2, ('H', 'R'): 0, ('W', 'H'): -2, ('H', 'N'): 1,
    ('W', 'T'): -2, ('T', 'T'): 5, ('S', 'F'): -2, ('W', 'P'): -4,
    ('L', 'D'): -4, ('B', 'I'): -3, ('L', 'H'): -3, ('S', 'N'): 1,
    ('B', 'T'): -1, ('L', 'L'): 4, ('Y', 'K'): -2, ('E', 'Q'): 2,
    ('Y', 'G'): -3, ('Z', 'S'): 0, ('Y', 'C'): -2, ('G', 'D'): -1,
    ('B', 'V'): -3, ('E', 'A'): -1, ('Y', 'W'): 2, ('E', 'E'): 5,
    ('Y', 'S'): -2, ('C', 'N'): -3, ('V', 'C'): -1, ('T', 'H'): -2,
    ('P', 'R'): -2, ('V', 'G'): -3, ('T', 'L'): -1, ('V', 'K'): -2,
    ('K', 'Q'): 1, ('R', 'A'): -1, ('I', 'R'): -3, ('T', 'D'): -1,
    ('P', 'F'): -4, ('I', 'N'): -3, ('K', 'I'): -3, ('M', 'D'): -3,
    ('V', 'W'): -3, ('W', 'W'): 11, ('M', 'H'): -2, ('P', 'N'): -2,
    ('K', 'A'): -1, ('M', 'L'): 2, ('K', 'E'): 1, ('Z', 'E'): 4,
    ('X', 'N'): -1, ('Z', 'A'): -1, ('Z', 'M'): -1, ('X', 'F'): -1,
    ('K', 'C'): -3, ('B', 'Q'): 0, ('X', 'B'): -1, ('B', 'M'): -3,
    ('F', 'C'): -2, ('Z', 'Q'): 3, ('X', 'Z'): -1, ('F', 'G'): -3,
    ('B', 'E'): 1, ('X', 'V'): -1, ('F', 'K'): -3, ('B', 'A'): -2,
    ('X', 'R'): -1, ('D', 'D'): 6, ('W', 'G'): -2, ('Z', 'F'): -3,
    ('S', 'Q'): 0, ('W', 'C'): -2, ('W', 'K'): -3, ('H', 'Q'): 0,
    ('L', 'C'): -1, ('W', 'N'): -4, ('S', 'A'): 1, ('L', 'G'): -4,
    ('W', 'S'): -3, ('S', 'E'): 0, ('H', 'E'): 0, ('S', 'I'): -2,
    ('H', 'A'): -2, ('S', 'M'): -1, ('Y', 'L'): -1, ('Y', 'H'): 2,
    ('Y', 'D'): -3, ('E', 'R'): 0, ('X', 'P'): -2, ('G', 'G'): 6,
    ('G', 'C'): -3, ('E', 'N'): 0, ('Y', 'T'): -2, ('Y', 'P'): -3,
    ('T', 'K'): -1, ('A', 'A'): 4, ('P', 'Q'): -1, ('T', 'C'): -1,
    ('V', 'H'): -3, ('T', 'G'): -2, ('I', 'Q'): -3, ('Z', 'T'): -1,
    ('C', 'R'): -3, ('V', 'P'): -2, ('P', 'E'): -1, ('M', 'C'): -1,
    ('K', 'N'): 0, ('I', 'I'): 4, ('P', 'A'): -1, ('M', 'G'): -3,
    ('T', 'S'): 1, ('I', 'E'): -3, ('P', 'M'): -2, ('M', 'K'): -1,
    ('I', 'A'): -1, ('P', 'I'): -3, ('R', 'R'): 5, ('X', 'M'): -1,
    ('L', 'I'): 2, ('X', 'I'): -1, ('Z', 'B'): 1, ('X', 'E'): -1,
    ('Z', 'N'): 0, ('X', 'A'): 0, ('B', 'R'): -1, ('B', 'N'): 3,
    ('F', 'D'): -3, ('X', 'Y'): -1, ('Z', 'R'): 0, ('F', 'H'): -1,
    ('B', 'F'): -3, ('F', 'L'): 0, ('X', 'Q'): -1, ('B', 'B'): 4
    }
    try:
        bloscore=blosum62[(aa_wt, aa_mut)]
    except:
        bloscore=blosum62[(aa_mut, aa_wt)]
    return(bloscore)




def binding_check(variant_database_object):
    binding_residues=Residue.objects.filter(variant=variant_database_object, binding_residue=True)
    if binding_residues.count()>0:
        return("Binding residue")
    else:
        return(None)

def tmh_check(variant_database_object):
    tmh=Tmh.objects.filter(meta_tmh=True, tmh_residue__residue__variant=variant_database_object)
    if tmh.count()>0:
        return("Tmh")
    else:
        return(None) 

def file_output(variant_object, neuro_status):
    file_dictionary={
    "Neurological":"neuro_disease_variants.csv",
    "Non-neurological":"non_neuro_disease_variants.csv"
    }

    physico_status=physico_change(variant_object.aa_wt, variant_object.aa_mut)

    tmh_status=tmh_check(variant_object)

    bind_status=binding_check(variant_object)

    line_for_file=f"{variant_object.variant_source_id},{variant_object.residue.protein.uniprot_id},{variant_object.residue.sequence_position},{variant_object.variant_source},{variant_object.aa_wt}, {variant_object.aa_mut},{physico_status}, {tmh_status}, {bind_status}\n"

    # print(line_for_file)
    try:
        with open(file_dictionary[neuro_status], "a") as file_object:
            # Append 'hello' at the end of file
            file_object.write(line_for_file)
    except FileNotFound:
        with open("unknown_neuro_status_disease_variant.csv", "a") as file_object:
            file_object.write(line_for_file)

def uniprot_both_list_check(uniprot_id, file_neuro_non_neuro):
    '''
    This aims to check the file to exclude ids that appear in neuro AND non neuro
    '''
    neuro=False
    non_neuro=False
    with open(file_neuro_non_neuro, "r") as filestream:
        for line in filestream:
            line=csv_parser(line) 
            if line[4]== uniprot_id:
                if line[1] == "Neurological":
                    neuro=True
                if line[1] == "Non-neurological":
                    non_neuro=True
    if neuro == True and non_neuro == True:
        return(True)
    else:
        return(False)           



disease_list="scripts/disease_split.csv"
def run():
    variant_objects=Variant.objects.filter(variant_source="ClinVar", disease_status="d").distinct('pk')
    for a_variant in variant_objects:
        #print(a_variant.residue.protein.uniprot_id)
        with open(disease_list, "r") as filestream:
            for line in filestream:
                line=csv_parser(line)
                if a_variant.residue.protein.uniprot_id == line[4] and uniprot_both_list_check(line[4], disease_list)==False:
                    #print(line[0],line[1])
                    file_output(a_variant, line[1]) 
                elif a_variant.residue.protein.uniprot_id == line[4] and uniprot_both_list_check(line[4], disease_list)==True:
                    print(f'{line[4]} ...\n... found in neuro and non neuro')
            
