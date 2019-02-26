from __future__ import division
import requests
import urllib.request
import numpy as np
import os
import subprocess
import re
import sys
import xml.etree.ElementTree as ET
from Bio import SeqIO
#from django.db import models
#from django.conf import settings
#from django.utils import timezone
from django.db import models
from tmh_db.models import Protein


def clinvar_check(tmh_info):
    '''
    Checks if a tmh has any variants in the variant file and spews out a list of
    variants and their position in the tmh.
    '''
    var_source = "ClinVar"
    variants_in_tmh = []
    for var_database_entry in var_results:
        # print(var_database_entry)
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

        except(IndexError):
            #print("Not enough datapoints in line.")
            pass
            # This list should get bigger as scores etc are added.

        # This is a really messy way to get the structure of the tmh data.
        # But heck, I'd rather be human readable than have less code!
        # tmh_info=[query_id, tmh_start, tmh_stop, tmh_topology, , membrane_location, n_ter_seq, tmh_sequence, c_ter_seq, tmh_number, len(list_of_tmhs)]

        query_id = str(tmh_info[0])
        tmh_start = int(tmh_info[1])
        tmh_stop = int(tmh_info[2])
        tmh_topology = str(tmh_info[3]).strip()
        evidence_type = str(tmh_info[4])
        n_location = str(tmh_info[5])
        n_ter_seq = str(tmh_info[6])
        tmh_seq = str(tmh_info[7])
        c_ter_seq = str(tmh_info[8])
        tmh_number = int(tmh_info[9])
        tmd_total = int(tmh_info[10])

        var_record_location = str(SEQ_NO)
        var_record_id = USER_ID
        AA_CHANGE

        '''
        This could be worth investigating if isoforms are an issue
        # ISOFORMS!!!!
        # This bit is fiddly since there are isoforms. First, we need to establish if the record is the right line to save some time.
        id_match = False
        if query_id == UNIPROT_ACCESSION:

            var_record_id = UNIPROT_ACCESSION
            var_record_location = SEQ_NO
            id_match = True

        elif query_id in ALL_TRANSCRIPTS:
            id_match = True

            # ' / ' deliniates isoforms. ',' deliniates items in isoforms.
            # Example:
            #   ENST00000379389,-,P05161,21,S/N,*,ISG15,0,0.73,Missense variant / ENST00000458555,-,P05161,-,-,*,ISG15,0,0.73,Upstream gene variant

            list_of_transcripts = str(ALL_TRANSCRIPTS).split(" / ")
            for isoform in list_of_transcripts:
                this_isoform = isoform.split(",")
                if query_id == this_isoform[2]:
                    var_record_id = this_isoform[2]
                    var_record_location = this_isoform[3]
                else:
                    pass
        else:
            pass
        '''
        if query_id == UNIPROT_ACCESSION:
            variant_type = "Unknown"
            variant_review = "Unknown"

            # Is the variant disease causing?

            for i in clinvar_lines:
                #print("Is", int(i[-1]), "equal to", int(USER_ID), "?" )
                if int(i[-1]) == int(var_record_id):  #  (variant id is last column in summary)
                    # print("clinvar summary and snipclip finally found a hit for variant ",int(var_record_id))

                    variant_type = i[6]
                    variant_review = i[24]
                    #print(variant_type, variant_review)

                    # AlleleID
                    # Type
                    # Name
                    # GeneID
                    # GeneSymbol
                    # HGNC_ID
                    # ClinicalSignificance
                    # ClinSigSimple
                    # LastEvaluated
                    # RS# (dbSNP)
                    #nsv/esv (dbVar)
                    # RCVaccession
                    # PhenotypeIDS
                    # PhenotypeList
                    # Origin
                    # OriginSimple
                    # Assembly
                    # ChromosomeAccession
                    # Chromosome
                    # Start
                    # Stop
                    # ReferenceAllele
                    # AlternateAllele
                    # Cytogenetic
                    # ReviewStatus
                    # NumberSubmitters
                    # Guidelines
                    # TestedInGTR
                    # OtherIDs
                    # SubmitterCategories
                    # VariationID

                    # n_ter_seq_all =  str(n_ter_seq + n_ter_seq_all)
                    # tmh_seq_all = str(tmh_seq + tmh_seq_all)
                    # c_ter_seq_all = str(c_ter_seq + c_ter_seq_all)

                    # Lets sort the flanks into inside outside
                    if tmh_topology == str("Inside"):
                        in_seq = str(n_ter_seq)
                        out_seq = str(c_ter_seq)
                    elif tmh_topology == str("Outside"):
                        in_seq = str(c_ter_seq)
                        out_seq = str(n_ter_seq)
                    else:
                        # print("Topology missing.")
                        pass

                    if "-" not in str(AA_CHANGE):
                        aa_change = AA_CHANGE.split("/")
                        starting_residue = str(aa_change[0])
                        mutation_residue = str(aa_change[1])
                    else:
                        starting_residue = "-"
                        mutation_residue = "-"
                    # We want as much information to be passed onto the next table.

                    variant_class = str(disease_class(variant_type))
                    # Z coordinate is the absolute distance from the central TMH residue to the variant position. I wonder if direction information is useful at this point.
                    if var_record_location == "-":  # Lets do some exceotions without "try"
                        pass
                    else:
                        z_coord_n_c = int(
                            int(tmh_start + (abs(tmh_stop - tmh_start) / 2)) - int(var_record_location))
                        if tmh_topology == str("Inside"):
                            z_coord_in_out = z_coord_n_c
                        elif tmh_topology == str("Outside"):
                            z_coord_in_out = 0 - z_coord_n_c

                        #print (int(var_record_location), tmh_start, int(var_record_location), tmh_stop)
                        if int(var_record_location) > tmh_start and int(var_record_location) <= tmh_stop:
                            variant = [query_id, var_record_location, starting_residue, mutation_residue,
                                       "TM", z_coord_n_c, z_coord_in_out, tmh_number, tmd_total, var_source, variant_class, variant_type, variant_review]
                            variants_in_tmh.append(variant)
                            print_list(variant)
                        elif int(var_record_location) >= tmh_start - len(n_ter_seq) and int(var_record_location) <= tmh_start:
                            if tmh_topology == str("Inside"):
                                variant = [query_id, var_record_location, starting_residue, mutation_residue,
                                           "N-terminal inside flank", z_coord_n_c, z_coord_in_out, tmh_number, tmd_total, var_source, variant_class, variant_type, variant_review]
                                variants_in_tmh.append(variant)
                                print_list(variant)
                            elif tmh_topology == str("Outside"):
                                variant = [query_id, var_record_location, starting_residue, mutation_residue,
                                           "N-terminal outside flank", z_coord_n_c, z_coord_in_out, tmh_number, tmd_total, var_source, variant_class, variant_type, variant_review]
                                variants_in_tmh.append(variant)
                                print_list(variant)
                        elif int(var_record_location) >= tmh_stop and int(var_record_location) <= tmh_stop + len(c_ter_seq):
                            if tmh_topology == str("Inside"):
                                variant = [query_id, var_record_location, starting_residue, mutation_residue,
                                           "C-terminal outside flank", z_coord_n_c, z_coord_in_out, tmh_number, tmd_total, var_source, variant_class, variant_type, variant_review]
                                variants_in_tmh.append(variant)
                                print_list(variant)
                            elif tmh_topology == str("Outside"):
                                variant = [query_id, var_record_location, starting_residue, mutation_residue,
                                           "C-terminal inside flank", z_coord_n_c, z_coord_in_out, tmh_number, tmd_total, var_source, variant_class, variant_type, variant_review]
                                variants_in_tmh.append(variant)
                                print_list(variant)

    return(variants_in_tmh)


def humsavar_var_check(tmh_info):
    var_source = "Humsavar"
    variants_in_tmh = []
    query_id = str(tmh_info[0])
    tmh_start = int(tmh_info[1])
    tmh_stop = int(tmh_info[2])
    tmh_topology = str(tmh_info[3]).strip()
    evidence_type = str(tmh_info[4])
    n_location = str(tmh_info[5])
    n_ter_seq = str(tmh_info[6])
    tmh_seq = str(tmh_info[7])
    c_ter_seq = str(tmh_info[8])
    tmh_number = int(tmh_info[9])
    tmd_total = int(tmh_info[10])

    filename = str("uniprot_bin/" + query_id + ".txt")
    input_format = "swiss"
    feature_type = "TRANSMEM"
    subcellular_location = "TOPO_DOM"

    # print("Checking ", query_id, "in humsavar.txt.")
    for record in SeqIO.parse(filename, input_format):
        tmhs = 0
        for i, feature in enumerate(record.features):
            if feature.type == 'VARIANT':
                for entry in humsavar_list:
                    # print(entry)
                    if str(entry[2]) == str(feature.id):
                        #variant_types=[str('Disease'), str('Polymorphism'), str('Unclassified')]
                        variant_type = str(entry[4])
                        variant_review = "SwissProt"

                        if tmh_start <= feature.location.end and tmh_stop >= feature.location.start:
                            # print("Transmembrane variant!")
                            # print(original_residue,"to", variant_residue)
                            # print(record.id)
                            tm_variant = True

                            #variant_record(original_residue, variant_residue)

                            var_record_location = str(feature.location.end)

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
                                        original_residue = str(feature.qualifiers)[char_num - 2]
                                        variant_residue = str(feature.qualifiers)[char_num + 3]

                            # We want as much information to be passed onto the next table.
                            variant_class = str(disease_class(variant_type))
                            # Z coordinate is the absolute distance from the central TMH residue to the variant position. I wonder if direction information is useful at this point.
                            z_coord_n_c = int(
                                int(tmh_start + (abs(tmh_stop - tmh_start) / 2)) - int(var_record_location))
                            if tmh_topology == str("Inside"):
                                z_coord_in_out = z_coord_n_c
                            elif tmh_topology == str("Outside"):
                                z_coord_in_out = 0 - z_coord_n_c

                            if int(var_record_location) > tmh_start and int(var_record_location) <= tmh_stop:
                                variant = [query_id, var_record_location, original_residue, variant_residue, "TM", z_coord_n_c,
                                           z_coord_in_out, tmh_number, tmd_total, var_source, variant_class, variant_type, variant_review]
                                variants_in_tmh.append(variant)
                                print_list(variant)
                            elif int(var_record_location) >= tmh_start - len(n_ter_seq) and int(var_record_location) <= tmh_start:
                                if tmh_topology == str("Inside"):
                                    variant = [query_id, var_record_location, original_residue, variant_residue, "N-terminal inside flank",
                                               z_coord_n_c, z_coord_in_out, tmh_number, tmd_total, var_source, variant_class, variant_type, variant_review]
                                    variants_in_tmh.append(variant)
                                    print_list(variant)
                                elif tmh_topology == str("Outside"):
                                    variant = [query_id, var_record_location, original_residue, variant_residue, "N-terminal outside flank",
                                               z_coord_n_c, z_coord_in_out, tmh_number, tmd_total, var_source, variant_class, variant_type, variant_review]
                                    variants_in_tmh.append(variant)
                                    print_list(variant)
                            elif int(var_record_location) >= tmh_stop and int(var_record_location) <= tmh_stop + len(c_ter_seq):
                                if tmh_topology == str("Inside"):
                                    variant = [query_id, var_record_location, original_residue, variant_residue, "C-terminal outside flank",
                                               z_coord_n_c, z_coord_in_out, tmh_number, tmd_total, var_source, variant_class, variant_type, variant_review]
                                    variants_in_tmh.append(variant)
                                    print_list(variant)
                                elif tmh_topology == str("Outside"):
                                    variant = [query_id, var_record_location, original_residue, variant_residue, "C-terminal inside flank",
                                               z_coord_n_c, z_coord_in_out, tmh_number, tmd_total, var_source, variant_class, variant_type, variant_review]
                                    variants_in_tmh.append(variant)
                                    print_list(variant)

    return(variants_in_tmh)


### Canonical script starts here ###
def run():
    var_file = str(sys.argv[1]) # Snipclip file
    # query_list = str(sys.argv[2]) # A lit of uniprot ids... This I guess should a DB query for all the protein IDs
    query_list = list(Protein.objects.order_by('uniprot_id'))
    print(query_list)

    clinvar_file = 'variant_summary.txt'

    print("Loading optimised flat file structure to memory...")
    # Grab the query input list file into a python list
    input_file = query_list
    #input_query = input_file.readlines()


    # humsavar table
    humsavar_list = []
    with open('humsavar.txt') as f:
        lines = f.read().splitlines()
        for i in lines:
            i = i.replace('  ', ' ')
            humsavar_list.append(i.split())
    # print(humsavar_list)

    # Load the  varsite tsv file from snip clip.
    var_results = []
    with open(var_file) as inputfile:
        for line_number, line in enumerate(inputfile):
            if line_number == 0:
                pass
            else:
                for a_query in input_query:
                    if clean_query(a_query) in str(line):
                        var_results.append(line.strip().split('\t'))
    # print(var_results)


    # Load the clinvar summary file
    clinvar_lines = []
    with open(clinvar_file) as inputfile:
        for line_number, line in enumerate(inputfile):
            if line_number == 0:
                pass
            else:
                for variant in var_results:
                    # print(str(variant[26])) # This should be simply the file ID
                    if clean_query(variant[26]) in str(line):
                        clinvar_lines.append(line.strip().split('\t'))
