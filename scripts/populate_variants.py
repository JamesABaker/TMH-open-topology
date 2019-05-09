from __future__ import division
import requests
import urllib
from requests import get
import shutil
import numpy as np
import os
import collections
import time
import gzip
import subprocess
import json
from subprocess import check_output
import re
import sys
import defusedxml.ElementTree as ET
import Bio
from Bio import SeqIO
from Bio import SwissProt
from Bio.SeqUtils.ProtParam import ProteinAnalysis
# env LDFLAGS="-I/usr/local/opt/openssl/include -L/usr/local/opt/openssl/lib" pip install psycopg2
from django.conf import settings
from django.db import models
from tmh_db.models import Database_Metadata, Uniref, Go, Structure, Structural_residue, Funfam_residue, Funfamstatus, Protein, Residue, Tmh, Tmh_deltag, Tmh_hydrophobicity, Tmh_residue, Tmh_tmsoc, Variant, Keyword, Binding_residue
from datetime import datetime, timedelta
from django.utils import timezone
from datetime import date
import pytz
from lxml import etree


print("Usage:\npython manage.py runscript populate --traceback")

# How many days should be allowed to not enforce updates
time_threshold = 7
today = date.today()
todaysdate = today.strftime("%d_%m_%Y")


def download(url, file_name):
    '''
    Downloads the content of a url to a local file.
    '''
    # open in binary mode
    with open(file_name, "wb") as file:
        # get request
        response = get(url)
        # write to file
        file.write(response.content)


def clean_query(query):
    '''
    This aims to generate a clean ascii query of a viable UniProt ID from a
     dirty input like a user input.
    '''

    illegal_characters = ["!", "\n", " ", "@"]
    for char in illegal_characters:
        query = query.replace(char, "")
    a_clean_query = query
    # print("Clean query result:", a_clean_query)
    return(a_clean_query)


def input_query_process(input_query):
    input_queries = []
    for query_number, a_query in enumerate(input_query):
        a_query = clean_query(a_query)
        print("Checking cache/downloading", a_query, ",",
              query_number + 1, "of", len(input_query), "records...")

        input_queries.append(a_query)

    input_query_set = set(input_queries)
    return([input_queries, input_query_set])

def disease_class(disease_type):
    '''
    Sorts:
        ?Affects
        ?association
        #Benign
        #Benign/Likely_benign
        ?Conflicting_interpretations_of_pathogenicity
        drug_response
        #Likely_benign
        *Likely_pathogenic
        ?no_interpretation_for_the_single_variant
        ?not_provided
        ?other
        *Pathogenic
        *Pathogenic/Likely_pathogenic
        ?protective
        ?risk_factor
        ?Uncertain_significance
    '''
    # Sometimes spaces are used instead of "_" s.
    disease_type = str(disease_type.replace(" ", "_"))
    disease = ["Disease", "Likely_pathogenic",
               "Pathogenic", "Pathogenic/Likely_pathogenic"]
    benign = ["Unclassified", "Polymorphism", "Affects", "association", "Benign", "Benign/Likely_benign", "Likely_benign", "Conflicting_interpretations_of_pathogenicity",
              "drug_response", "no_interpretation_for_the_single_variant", "not_provided",  "other", "protective", "risk_factor", "Uncertain_significance"]
    if str(disease_type) in disease:
        pathogenicity = "d"
    elif str(disease_type) in benign:
        pathogenicity = "n"
    else:
        pathogenicity = "u"
        print("Unknown pathogenicity:", str(disease_type))
    return(pathogenicity)


def clinvar_variant_check(clinvar_variants, clinvar_summary):
    '''
    Checks if a tmh has any variants in the variant file and spews out a list of
    variants and their position in the tmh. Unlike the generic gnomAD function,
    this crossreferences clinvar VarMap tsv against the clinvar summary file.
    '''

    variant_source = "ClinVar"
    variants_in_tmh = []
    var_database_entry = clinvar_variants

    try:
        # Yes, I know all caps is bad, but this is just way easier than reformatting every time SnipClip headers change.
        CHROMOSOME = str(var_database_entry[0])
        COORDS = str(var_database_entry[1])
        USER_BASE = str(var_database_entry[2])
        USER_VARIANT = str(var_database_entry[3])
        ENSEMBL_BASE = str(var_database_entry[4])
        VEP_CODING_BASE = str(var_database_entry[5])
        GENE = str(var_database_entry[6])
        GENE_ACC = str(var_database_entry[7])
        REFSEQ_GENE_ACC = str(var_database_entry[8])
        TRANSCRIPT = str(var_database_entry[9])
        REFSEQ_TRANSCRIPT = str(var_database_entry[10])
        STRAND_DIR = str(var_database_entry[11])
        CODON_CHANGE = str(var_database_entry[12])
        VEP_AA = str(var_database_entry[13])
        UNIPROT_AA = str(var_database_entry[14])
        AA_CHANGE = str(var_database_entry[15])
        POLYPHEN_SCORE = str(var_database_entry[16])
        SIFTS_SCORE = str(var_database_entry[17])
        UNIPROT_ACCESSION = str(var_database_entry[18])
        PROTEIN_NAME = str(var_database_entry[19])
        SEQ_NO = str(var_database_entry[20])
        CHANGE_TYPE = str(var_database_entry[21])
        ALL_TRANSCRIPTS = str(var_database_entry[22])
        NOTE = str(var_database_entry[23])
        GNOMAD_AF = str(var_database_entry[24])
        NEGATIVE = str(var_database_entry[25])
        USER_ID = str(var_database_entry[26])
        SYNONYMOUS = str(var_database_entry[27])
        HAVE_PDB = str(var_database_entry[28])
        PDB_UNIPROT_MATCH = str(var_database_entry[29])
        CLOSEST_PDB_CODE = str(var_database_entry[30])
        PDB_CHAIN = str(var_database_entry[31])
        PDB_PROTEIN_NAME = str(var_database_entry[32])
        PDB_EXPT_TYPE = str(var_database_entry[33])
        PDB_RESOLUTION = str(var_database_entry[34])
        PDB_RFACT = str(var_database_entry[35])
        PDB_UNIPROT_ACC = str(var_database_entry[36])
        PDB_IDENTITY = str(var_database_entry[37])
        PDB_SW_SCORE = str(var_database_entry[38])
        PDB_E_VALUE = str(var_database_entry[39])
        RES_NAME = str(var_database_entry[40])
        RES_NUM = str(var_database_entry[41])
        SST = str(var_database_entry[42])
        CAT_RES = str(var_database_entry[43])
        DISULPHIDE = str(var_database_entry[44])
        NTO_DNA = str(var_database_entry[45])
        NTO_LIGAND = str(var_database_entry[46])
        NTO_METAL = str(var_database_entry[47])
        NTO_PROTEIN = str(var_database_entry[48])
        NPDB_RES = str(var_database_entry[49])
        LIGANDS = str(var_database_entry[50])
        METALS = str(var_database_entry[51])
        PFAM_DOMAIN = str(var_database_entry[52])
        PFAM_NAME = str(var_database_entry[53])
        CATH_DOMAIN = str(var_database_entry[54])
        CATH_NAME = str(var_database_entry[55])
        RES_CONSERVATION = str(var_database_entry[56])
        NCONS_SEQS = str(var_database_entry[57])
        DISEASES = str(var_database_entry[58])
        DISEASE_VARIANTS = str(var_database_entry[59])
        NVARIANTS = str(var_database_entry[60])
        NAT_VARIANTS = str(var_database_entry[61])

        var_record_location = SEQ_NO
        variant_source_id = USER_ID
        uniprot_record = UNIPROT_ACCESSION
        variant_type = "Unknown"
        variant_review = "Unknown"
        aa_wt = UNIPROT_AA
        # VarMap shows the change as X/N
        if len(AA_CHANGE) == 3 and "/" in AA_CHANGE:
            aa_mut = AA_CHANGE.split("/")[1]
        else:
            aa_mut = AA_CHANGE
        disease_status = ""
        disease_comments = ""

        # Is the variant disease causing?

        for i in clinvar_summary:
                #print("Is", int(i[-1]), "equal to", int(USER_ID), "?" )
            if int(i[-1]) == int(variant_source_id):  #  (variant id is last column in summary)
                    # print("clinvar summary and snipclip finally found a hit for variant ",int(var_record_id))

                disease_status = disease_class(i[6])
                disease_comments = i[24]

        var_to_database(uniprot_record, var_record_location, aa_wt,
                        aa_mut, disease_status, disease_comments, variant_source, variant_source_id)

    except(IndexError):
        #print("Not enough datapoints in line.")
        pass


def gnomad_variant_check(gnomad_variants):
    '''
    Checks if a tmh has any variants in the variant file and spews out a list of
    variants and their position in the tmh.
    '''

    # This could be worth investigating if isoforms are an issue
    # ISOFORMS!!!!
    # This bit is fiddly since there are isoforms. First, we need to establish if the record is the right line to save some time.
    #id_match = False
    # if query_id == UNIPROT_ACCESSION:

    #    var_record_id = UNIPROT_ACCESSION
    #    var_record_location = SEQ_NO
    #    id_match = True

    # elif query_id in ALL_TRANSCRIPTS:
    #    id_match = True

    #    # ' / ' deliniates isoforms. ',' deliniates items in isoforms.
    #    # Example:
    #    #   ENST00000379389,-,P05161,21,S/N,*,ISG15,0,0.73,Missense variant / ENST00000458555,-,P05161,-,-,*,ISG15,0,0.73,Upstream gene variant

    #    list_of_transcripts = str(ALL_TRANSCRIPTS).split(" / ")
    #    for isoform in list_of_transcripts:
    #        this_isoform = isoform.split(",")
    #        if query_id == this_isoform[2]:
    #            var_record_id = this_isoform[2]
    #            var_record_location = this_isoform[3]
    #        else:
    #            pass
    # else:
    #    pass

    try:
        variant_source = "gnomAD"
        var_database_entry = gnomad_variants

        # Yes, I know all caps is bad, but this is just way easier than reformatting every time SnipClip headers change.
        CHROMOSOME = str(var_database_entry[0])
        COORDS = str(var_database_entry[1])
        USER_BASE = str(var_database_entry[2])
        USER_VARIANT = str(var_database_entry[3])
        ENSEMBL_BASE = str(var_database_entry[4])
        VEP_CODING_BASE = str(var_database_entry[5])
        GENE = str(var_database_entry[6])
        GENE_ACC = str(var_database_entry[7])
        REFSEQ_GENE_ACC = str(var_database_entry[8])
        TRANSCRIPT = str(var_database_entry[9])
        REFSEQ_TRANSCRIPT = str(var_database_entry[10])
        STRAND_DIR = str(var_database_entry[11])
        CODON_CHANGE = str(var_database_entry[12])
        VEP_AA = str(var_database_entry[13])
        UNIPROT_AA = str(var_database_entry[14])
        AA_CHANGE = str(var_database_entry[15])
        POLYPHEN_SCORE = str(var_database_entry[16])
        SIFTS_SCORE = str(var_database_entry[17])
        UNIPROT_ACCESSION = str(var_database_entry[18])
        PROTEIN_NAME = str(var_database_entry[19])
        SEQ_NO = str(var_database_entry[20])
        CHANGE_TYPE = str(var_database_entry[21])
        ALL_TRANSCRIPTS = str(var_database_entry[22])
        NOTE = str(var_database_entry[23])
        GNOMAD_AF = str(var_database_entry[24])
        NEGATIVE = str(var_database_entry[25])
        USER_ID = str(var_database_entry[26])
        SYNONYMOUS = str(var_database_entry[27])
        HAVE_PDB = str(var_database_entry[28])
        PDB_UNIPROT_MATCH = str(var_database_entry[29])
        CLOSEST_PDB_CODE = str(var_database_entry[30])
        PDB_CHAIN = str(var_database_entry[31])
        PDB_PROTEIN_NAME = str(var_database_entry[32])
        PDB_EXPT_TYPE = str(var_database_entry[33])
        PDB_RESOLUTION = str(var_database_entry[34])
        PDB_RFACT = str(var_database_entry[35])
        PDB_UNIPROT_ACC = str(var_database_entry[36])
        PDB_IDENTITY = str(var_database_entry[37])
        PDB_SW_SCORE = str(var_database_entry[38])
        PDB_E_VALUE = str(var_database_entry[39])
        RES_NAME = str(var_database_entry[40])
        RES_NUM = str(var_database_entry[41])
        SST = str(var_database_entry[42])
        CAT_RES = str(var_database_entry[43])
        DISULPHIDE = str(var_database_entry[44])
        NTO_DNA = str(var_database_entry[45])
        NTO_LIGAND = str(var_database_entry[46])
        NTO_METAL = str(var_database_entry[47])
        NTO_PROTEIN = str(var_database_entry[48])
        NPDB_RES = str(var_database_entry[49])
        LIGANDS = str(var_database_entry[50])
        METALS = str(var_database_entry[51])
        PFAM_DOMAIN = str(var_database_entry[52])
        PFAM_NAME = str(var_database_entry[53])
        CATH_DOMAIN = str(var_database_entry[54])
        CATH_NAME = str(var_database_entry[55])
        RES_CONSERVATION = str(var_database_entry[56])
        NCONS_SEQS = str(var_database_entry[57])
        DISEASES = str(var_database_entry[58])
        DISEASE_VARIANTS = str(var_database_entry[59])
        NVARIANTS = str(var_database_entry[60])
        NAT_VARIANTS = str(var_database_entry[61])

        var_record_location = SEQ_NO
        var_record_id = USER_ID
        uniprot_record = UNIPROT_ACCESSION
        variant_type = "Unknown"
        variant_review = "Unknown"
        aa_wt = UNIPROT_AA
        if len(AA_CHANGE) == 3 and "/" in AA_CHANGE:
            aa_mut = AA_CHANGE.split("/")[1]
        else:
            aa_mut = AA_CHANGE
        disease_status = "gnomAD"
        disease_comments = ""
        user_id = USER_ID

        # Is the variant disease causing?

        # for i in clinvar_summary:
        #        #print("Is", int(i[-1]), "equal to", int(USER_ID), "?" )
        #    if int(i[-1]) == int(var_record_id):  #  (variant id is last column in summary)
        #            # print("clinvar summary and snipclip finally found a hit for variant ",int(var_record_id))

        #        disease_status = i[6]
        #        disease_comments = i[24]

        var_to_database(uniprot_record, var_record_location, aa_wt,
                        aa_mut, disease_status, disease_comments, variant_source, user_id)

    except(IndexError):
        #print("Not enough datapoints in line.")
        pass


def humsavar_variant_check(humsavar_variant):
    print(humsavar_variant)
    humsavar_gene = humsavar_variant[0]
    uniprot_record = humsavar_variant[1]
    humsavar_variant_id = humsavar_variant[2]
    humsavar_variant_change = humsavar_variant[3]
    humsavar_variant_disease_type = humsavar_variant[4]
    humsavar_variant_gene_position = humsavar_variant[5]
    humsavar_variant_comment = humsavar_variant[6]

    variant_source = "Humsavar"
    filename = str("scripts/external_datasets/uniprot_bin/" +
                   uniprot_record + ".txt")
    input_format = "swiss"
    subcellular_location = "TOPO_DOM"

    # print("Checking ", query_id, "in humsavar.txt.")
    for record in SeqIO.parse(filename, input_format):
        for i, feature in enumerate(record.features):
            if feature.type == 'VARIANT':
                if str(humsavar_variant_id) == str(feature.id):
                    #variant_types=[str('Disease'), str('Polymorphism'), str('Unclassified')]
                    variant_review = "SwissProt"

                    for char_num, char in enumerate(str(feature.qualifiers)):
                        if char == "-":
                            # This is some hideous code that will break at the slightest change to how variants are sorted.
                            # FT   VARIANT     838    838       R -> H (in CORD6; dbSNP:rs61750173). is a ypical line that
                            # Bio parses to {'description': 'R -> G (in dbSNP:rs742493).
                            # {ECO:0000269|PubMed:14769797, ECO:0000269|PubMed:15489334}.'}. Here I take advantage of the
                            # preceding "'" and proceding " " to identify point changes in variants.
                            # Before we figure if it's TRANSMEM or not, here, we catch the variant for point mutations.
                            # At some point this needs to be rewritted to handle other types of variant.

                            if "->" in str(feature.qualifiers) and str(feature.qualifiers)[char_num + 1] == ">" and str(feature.qualifiers)[char_num - 3] == "'" and str(feature.qualifiers)[char_num + 4] == " ":
                                # print(feature.id)
                                # print(feature.qualifiers)
                                aa_wt = str(feature.qualifiers)[char_num - 2]
                                aa_mut = str(feature.qualifiers)[char_num + 3]

                                # We want as much information to be passed onto the next table.
                                disease_status = str(disease_class(
                                    humsavar_variant_disease_type))
                                disease_comments = str(
                                    humsavar_variant_disease_type + ";" + humsavar_variant_comment)
                                var_record_location = feature.location.start + 1  # This might need +1?
                                variant_source_id = humsavar_variant_id
                                var_to_database(uniprot_record, var_record_location, aa_wt,
                                                aa_mut, disease_status, disease_comments, variant_source, variant_source_id)


def var_to_database(uniprot_record, var_record_location, aa_wt, aa_mut, disease_status, disease_comments, variant_source, variant_source_id):
    '''
    Adds the variant from various external databases and flat files to the database in a standardised way.
    '''
    if var_record_location == "-":
        print("Unkown sequence location. Possibly intron: ", uniprot_record, var_record_location,
              aa_wt, "->", aa_mut, disease_status, disease_comments, variant_source, variant_source_id)
    elif aa_wt == "-":
        print("Wildtype amino acid not defined. Assuming this is not an SNP: ", uniprot_record,
              var_record_location, aa_wt, "->", aa_mut, disease_status, disease_comments, variant_source, variant_source_id)
    elif len(aa_wt) > 1 or len(aa_mut) > 1:
        print("More than a single residue changed. Assuming this is not an SNP: ", uniprot_record,
              var_record_location, aa_wt, "->", aa_mut, disease_status, disease_comments, variant_source, variant_source_id)
    else:

        protein = Protein.objects.get(uniprot_id=uniprot_record)

        if int(var_record_location) <= len(str(protein.full_sequence)):

            residue_variant = Residue.objects.get(
                protein=protein, sequence_position=var_record_location)
            if str(residue_variant.amino_acid_type) == str(aa_wt):
                print("Adding ", uniprot_record, var_record_location, aa_wt, "->", aa_mut,
                      disease_status, disease_comments, variant_source, "to database variant table.")
                record_for_database, created = Variant.objects.update_or_create(
                    residue=residue_variant,
                    aa_wt=aa_wt,
                    aa_mut=aa_mut,

                    disease_status=disease_status,
                    disease_comments=disease_comments,
                    variant_source=variant_source,
                    variant_source_id=variant_source_id,
                    defaults={
                    }
                )
            else:
                print("Mismatch between wild-type amino acids. UniProt:", str(residue_variant.amino_acid_type), str(variant_source), ":", str(
                    aa_wt), "for record", uniprot_record, var_record_location, aa_wt, "->", aa_mut, disease_status, disease_comments, variant_source)
        else:
            print("Variant position exceeds the length of the protein. Protein length:", len(str(protein.full_sequence)), "Variant position:",
                  var_record_location, "for record", uniprot_record, var_record_location, aa_wt, "->", aa_mut, disease_status, disease_comments, variant_source)



def run():
    '''
    This is what django runs. This is effectively the canonical script,
    even though django forces it to be in a function.
    This will go through several databases and extract the TMH boundaries from proteins,
    and then identify which variants are in those TMHs.
    $ python3 manage.py runscript populate --traceback
    '''

    ### Canonical script starts here ###

    # In full scale mode it will take a long time which may not be suitable for development.
    #input_query = get_uniprot()
    # Here we will just use a watered down list of tricky proteins. Uncomment this line for testing the whole list.
    input_query = ["P01850", "P22760", "Q5K4L6","Q7Z5H4", "P32897", "Q9NR77", "P31644", "Q9NS61"]


    # Also, parse the variant files which can be massive.
    # humsavar table
    print(input_query)
    print("Starting TMH database population script...")
    os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'tmh_database.settings')

    ### Download uniprot files ###
    inputs = input_query_process(input_query)
    input_queries = inputs[0]
    input_query_set = inputs[1]
    ### Variant tables ###

    print("Reading the variant tables...")

    ## Humsavar ##

    humsavar_file = "scripts/external_datasets/humsavar.txt"
    st = os.stat(humsavar_file)
    humsavar_file_age = (time.time() - st.st_mtime)
    print("Downloading humsavar.txt from UniProt.")
    humsavar_url = 'https://www.uniprot.org/docs/humsavar.txt'

    if humsavar_file_age > time_threshold:
        download(humsavar_url, humsavar_file)

    humsavar_list = []
    with open(humsavar_file) as f:
        lines = f.read().splitlines()
        for line_number, i in enumerate(lines):

            i = i.replace('  ', ' ')
            humsavar_variant = i.split()
            if line_number > 50 and line_number < len(lines) - 5:
                if humsavar_variant[1] in input_query_set:

                    # fixes issue 40. Some IDs are longer than the column width and use the space we are using to split.
                    if len(humsavar_variant) == 6:
                        if humsavar_variant[5][-1] == "-" and len(humsavar_variant[5]) > 1:
                            humsavar_variant[5] = humsavar_variant[5][:-1]
                            humsavar_variant.append(humsavar_variant[5][-1:])
                    if len(humsavar_variant) > 8:
                        humsavar_variant[7: -
                                         1] = [''.join(humsavar_variant[7:-1])]

                    humsavar_list.append(humsavar_variant)

    for humsavar_variant in humsavar_list:
        humsavar_variant_check(humsavar_variant)
    # print(humsavar_list)

    ## ClinVar ##

    # Load the  varsite tsv file from snip clip.
    clinvar_results = []
    clinvar_results_set = set()

    print("Loading ClinVar tsv from VarMap to memory.")
    with open("scripts/external_datasets/clinvar_snipclipa13_02_2019.tsv") as inputfile:
        for line_number, var_database_entry in enumerate(inputfile):
            if line_number == 0:
                pass
            else:
                # print(var_database_entry)
                var_database_entry = var_database_entry.strip().split('\t')
                UNIPROT_ACCESSION = str(var_database_entry[18])
                USER_ID = str(var_database_entry[26])

                if str(UNIPROT_ACCESSION) in input_query_set:
                    print("Storing variant", USER_ID, "for",
                          UNIPROT_ACCESSION, "to memory.")
                    clinvar_results.append(var_database_entry)
                    clinvar_results_set.add(clean_query(str(USER_ID)))

    print(clinvar_results_set)
    print(len(clinvar_results),
          "ClinVar variants found in database that will be checked.")

    # Load the clinvar summary file
    clinvar_summary_lines = []
    print("Loading the variant summaries from ClinVar. This holds information on disease states in clinvar.")
    with open("scripts/external_datasets/variant_summary.txt") as inputfile:
        for line_number, summary_variant in enumerate(inputfile):
            if line_number > 0:
                summary_variant = summary_variant.strip().split('\t')
                # print(str(summary_variant[-1]))
                if clean_query(str(summary_variant[-1])) in clinvar_results_set:
                    clinvar_summary_lines.append(summary_variant)
            else:
                pass

    print(len(clinvar_summary_lines), "summaries fetched of",
          len(clinvar_results), "ClinVar variants.")

    # We now have a list of the clinvar variants and a list of the clinvar summary.
    # Other tsvs form VarMap hopefully won't need this.
    for clinvar_variant in clinvar_results:
        clinvar_variant_check(clinvar_variant, clinvar_summary_lines)

    ## gnomAD ##

    print("Loading gnomAD tsv file to memory. This may take a while...")
    gnomad_results = []
    with open("scripts/external_datasets/gnomAD_varsite.tsv") as inputfile:
        for line_number, var_database_entry in enumerate(inputfile):
            if line_number == 0:
                pass
            else:
                # print(var_database_entry)
                var_database_entry = var_database_entry.strip().split('\t')
                UNIPROT_ACCESSION = str(var_database_entry[18])

                if UNIPROT_ACCESSION in input_query_set:
                    gnomad_results.append(var_database_entry)

    print(len(gnomad_results),
          "variants relating to query list found in gnomAD. Adding SNPs to database...")

    for gnomad_variant in gnomad_results:
        gnomad_variant_check(gnomad_variant)
