from __future__ import division
import random
import os
import time

from datetime import date
from datetime import datetime
from datetime import timedelta

import pytz
from Bio import SeqIO
from django.conf import settings
from django.db import models
from django.utils import timezone

from scripts.populate_general_functions import *
# env LDFLAGS="-I/usr/local/opt/openssl/include -L/usr/local/opt/openssl/lib" pip install psycopg2
bulk_variants_to_add=[]
bulk_diseases=[]

def varmap_columns_and_keys(column_headers):
    column_headers = column_headers.split()
    varmap_col_dictionary = {}
    for column_number, column_title in enumerate(column_headers):
        varmap_col_dictionary[column_title] = column_number
    # print(varmap_col_dictionary)
    return varmap_col_dictionary

def var_to_database(uniprot_record, var_record_location, aa_wt, aa_mut, disease_status, disease_comments, variant_source, variant_source_id, diseases=[], germline=None):
    '''
    Adds the variant from various external databases and flat files to the database in a standardised way.
    '''
    global bulk_diseases
    global bulk_variants_to_add
    if var_record_location == "-":
        print("Unkown sequence location. Possibly intron: ", uniprot_record, var_record_location,
              aa_wt, "->", aa_mut, disease_status, disease_comments, variant_source, variant_source_id)
    elif aa_wt == "-":
        print("Wildtype amino acid not defined. Assuming this is not an SNP: ", uniprot_record,
              var_record_location, aa_wt, "->", aa_mut, disease_status, disease_comments, variant_source, variant_source_id)
    elif len(aa_wt) > 1 or len(aa_mut) > 1:
        print("More than a single residue changed. Assuming this is not an SNP: ", uniprot_record,
              var_record_location, aa_wt, "->", aa_mut, disease_status, disease_comments, variant_source, variant_source_id)
    elif "*" in str(aa_mut):
        print("Stop codon introduced. This will change more than one residue:", uniprot_record,
              var_record_location, aa_wt, "->", aa_mut, disease_status, disease_comments, variant_source, variant_source_id)
    else:

        print(f'updating {uniprot_record}')
        try:
            protein = Protein.objects.get(uniprot_id=uniprot_record)

            if int(var_record_location) <= len(str(protein.full_sequence)):

                residue_variant = Residue.objects.get(
                    protein=protein, sequence_position=var_record_location)
                if str(residue_variant.amino_acid_type) == str(aa_wt):
                    print("Adding ", uniprot_record, var_record_location, aa_wt, "->", aa_mut,
                          disease_status, disease_comments, variant_source, "to database variant table.")
                    #record_for_database, created = Variant.objects.update_or_create(
                    #    residue=residue_variant,
                    #    aa_wt=aa_wt,
                    #    aa_mut=aa_mut,

                    #    disease_status=disease_status,
                    #    disease_comments=disease_comments,
                    #    variant_source=variant_source,
                    #    variant_source_id=variant_source_id,
                    #    germline=germline,
                    #    defaults={
                    #    }
                    #)
                    bulk_variants_to_add.append(
                        Variant(
                            residue=residue_variant,
                            aa_wt=aa_wt,
                            aa_mut=aa_mut,

                            disease_status=disease_status,
                            disease_comments=disease_comments,
                            variant_source=variant_source,
                            variant_source_id=variant_source_id,
                            germline=germline,
                        )
                        )

                    '''
                    disease_variants=Variant.objects.filter(variant_source_id=variant_source_id, variant_source=variant_source).distinct('pk')
                
=======
                    disease_variants=Variant.objects.filter(variant_source_id=variant_source_id, variant_source=variant_source).distinct('pk')
>>>>>>> 7d0581ca1391f95dc71047acd29287fc53ec84f4
                    for var in disease_variants:
                        for dis in diseases:
                            record_for_database, created = Disease.objects.update_or_create(
                                disease_name=dis)
                            disease_type=Disease.objects.get(disease_name=dis) 
                            disease_type.implicated_variants.add(var)
<<<<<<< HEAD
                    '''
                else:
                    print("Mismatch between wild-type amino acids. UniProt:", str(residue_variant.amino_acid_type), str(variant_source), ":", str(
                        aa_wt), "for record", uniprot_record, var_record_location, aa_wt, "->", aa_mut, disease_status, disease_comments, variant_source)
            else:
                print("Variant position exceeds the length of the protein. Protein length:", len(str(protein.full_sequence)), "Variant position:",
                      var_record_location, "for record", uniprot_record, var_record_location, aa_wt, "->", aa_mut, disease_status, disease_comments, variant_source)
        except(Protein.DoesNotExist):
            print(f'{uniprot_record} was not found in database. Skipping entry.')
        if len(bulk_variants_to_add)>1000:
            Variant.objects.bulk_create(bulk_variants_to_add) 
            bulk_variants_to_add=[]                     
        '''
        if len(bulk_diseases)>1000:
            Disease.objects.bulk_create(bulk_diseases)
        '''

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
    disease = ["Disease", "Likely_pathogenic","Pathogenic", "Pathogenic/Likely_pathogenic"]
    benign = ["Benign", "Benign/Likely_benign", "Likely_benign", ]
    unknown = ["no_assertion_criteria_provided", "Unclassified", "Polymorphism", "Affects", "association", "Conflicting_interpretations_of_pathogenicity", "drug_response", "no_interpretation_for_the_single_variant", "not_provided",  "other", "protective", "risk_factor", "Uncertain_significance"]

    if str(disease_type) in disease:
        pathogenicity = "d"
    elif str(disease_type) in benign:
        pathogenicity = "n"
    elif str(disease_type) in unknown:
        pathogenicity = "u"
    else:
        print(f'Unknown pathogenic status: {disease_type}')
        pathogenicity = "e"
    return(pathogenicity)

def humsavar(uniprot_input_set):

    ## Humsavar ##
    time_threshold = 7

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
                if humsavar_variant[1] in uniprot_input_set:

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

    return()


def humsavar_variant_check(humsavar_variant):
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


def varmap_process(varmap_list, source, varmap_index):
    for varmap_item in varmap_list:
        varmap_item=varmap_item.strip().split('\t') 
        variant_source = source

        var_record_location = varmap_item[varmap_index["SEQ_NO"]]
        var_record_id = varmap_item[varmap_index["USER_ID"]]
        uniprot_record = varmap_item[varmap_index["UNIPROT_ACCESSION"]]
        if "-" in uniprot_record:
            uniprot_record=uniprot_record.split("-")
            uniprot_record=uniprot_record[0]
        variant_review = "Unknown"
        aa_wt = varmap_item[varmap_index["UNIPROT_AA"]]
        if len(varmap_item[varmap_index["AA_CHANGE"]]) == 3 and "/" in varmap_item[varmap_index["AA_CHANGE"]]:
            aa_mut = varmap_item[varmap_index["AA_CHANGE"]].split("/")[1]
        else:
            aa_mut = varmap_item[varmap_index["AA_CHANGE"]]
        disease_status = ""
        disease_comments = ""
        user_id = varmap_item[varmap_index["USER_ID"]]
        germline_status=None 
        # Updated ClinVar format
        parsed_id=varmap_clinvar_id_parse(user_id)
        disease_status = disease_class(''.join(parsed_id["disease_status"]))
        disease_list=parsed_id["diseases"]


        #print(uniprot_record, var_record_location, aa_wt, aa_mut, disease_status, variant_source, user_id)

        var_to_database(uniprot_record, var_record_location, aa_wt, aa_mut,
                        disease_status, disease_comments, variant_source, user_id, germline=germline_status, diseases=disease_list)

def varmap_clinvar_id_parse(varmap_id):
    varmap_id_dict={}
    varmap_id=varmap_id.split('yYy')
    varmap_clinvar_order={
        0:'qc',
        1:'clinvar_id', 
        2:'clinvar_allele_id',
        3:'diseases', 
        4:'disease_db_ids', 
        5:'CLNHGVS', 
        6: 'CLNREVSTAT', 
        7:'disease_status', 
        8:'CLNSIGCONF',
        9:'CLNVC',
        10:'CLNVCSO',
        11:'CLNVI',
        12:'DBVARID',
        13:'GENEINFO',
        14:'MC',
        15:'ORIGIN',
        16:'RS',
        17:'SSR',       
        19:'CLNDNINCL*',
        20:'CLNDISDBINCL*',
        21:'CLNSIGINCL*',
        22:'CLNVC'}

    for n, i in enumerate(varmap_id):
        i=i.split('zZz')
        if n in varmap_clinvar_order:
            varmap_id_dict[varmap_clinvar_order[n]]="".join(i)
    return(varmap_id_dict)

def varmap_header_dict(varmap_file):
    with open(varmap_file, encoding="ISO-8859-1") as inputfile:
        for line_number, var_database_entry in enumerate(inputfile):

            # This is the header line
            if line_number == 0:
                varmap_headers = varmap_columns_and_keys(var_database_entry)
    return(varmap_headers)


def run():
    input_query = input_query_get()
    # Also, parse the variant files which can be massive.
    # humsavar table
    print(input_query)
    print("Starting TMH database population script...")
    os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'tmh_database.settings')

    ### Download uniprot files ###
    inputs = input_query_process(input_query)
    #input_queries = inputs[0]
    input_query_set = inputs[1]

    ### VarMap ###

    varmap_files = {
        "clinvar": "scripts/external_datasets/clinvar_restitched.tsv",
        "clinvar_test": "scripts/external_datasets/test_clinvar.tsv",
    }

    clinvar_results_list=[]
    with open(varmap_files["clinvar"], encoding="ISO-8859-1") as inputfile:
        for line_number, var_database_entry in enumerate(inputfile):
            if line_number > 0:
                clinvar_results_list.append(var_database_entry)
    random.shuffle(clinvar_results_list)
    print(varmap_header_dict(varmap_files["clinvar"]))
    print(clinvar_results_list[0])
    varmap_process(clinvar_results_list, "ClinVar", varmap_header_dict(varmap_files["clinvar"]))
    Variant.objects.bulk_create(bulk_variants_to_add)

    ### Humsavar ###
    humsavar(input_query_set)

