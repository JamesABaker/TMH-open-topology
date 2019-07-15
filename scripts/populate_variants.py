from __future__ import division
import os
# env LDFLAGS="-I/usr/local/opt/openssl/include -L/usr/local/opt/openssl/lib" pip install psycopg2
from django.conf import settings
from django.db import models
from tmh_db.models import Database_Metadata, Uniref, Go, Structure, Structural_residue, Funfam_residue, Funfamstatus, Protein, Residue, Tmh, Tmh_deltag, Tmh_hydrophobicity, Tmh_residue, Tmh_tmsoc, Variant, Keyword, Binding_residue
from datetime import datetime, timedelta
from django.utils import timezone
from datetime import date
import pytz
from scripts.populate_general_functions import *
import time
import Bio
from Bio import SeqIO


def varmap_columns_and_keys(column_headers):
    column_headers = column_headers.split('\t')
    varmap_col_dictionary = {}
    for column_number, column_title in enumerate(column_headers):
        varmap_col_dictionary[column_title] = column_number
    # print(varmap_col_dictionary)
    return varmap_col_dictionary


def stripped_variant_list(varmap_file, input_query_set):

    varmap_results = []
    varmap_results_set = set()

    with open(varmap_file) as inputfile:
        for line_number, var_database_entry in enumerate(inputfile):

            # This is the header line
            if line_number == 0:
                varmap_headers = varmap_columns_and_keys(var_database_entry)
                uniprot_accession_index = varmap_headers["UNIPROT_ACCESSION"]
                varmap_user_id_index = varmap_headers["USER_ID"]
            else:
                var_database_entry = var_database_entry.strip().split('\t')
                uniprot_accesstion = str(var_database_entry[uniprot_accession_index])
                varmap_user_id = str(var_database_entry[varmap_user_id_index])

                if str(uniprot_accesstion) in input_query_set:
                    #print("Storing variant", varmap_user_id, "for", uniprot_accesstion, "to memory.")
                    varmap_results.append(var_database_entry)
                    varmap_results_set.add(clean_query(str(varmap_user_id)))
    return(varmap_results, varmap_results_set)


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
    elif "*" in str(aa_mut):
        print("Stop codon introduced. This will change more than one residue:", uniprot_record,
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


def varmap_process(varmap_list, source, varmap_index, *clinvar_summary):
    for varmap_item in varmap_list:
        variant_source = source
        variants_in_tmh = []
        var_record_location = varmap_item[varmap_index["SEQ_NO"]]
        var_record_id = varmap_item[varmap_index["USER_ID"]]
        uniprot_record = varmap_item[varmap_index["UNIPROT_ACCESSION"]]
        variant_type = "Unknown"
        variant_review = "Unknown"
        aa_wt = varmap_item[varmap_index["UNIPROT_AA"]]
        if len(varmap_item[varmap_index["AA_CHANGE"]]) == 3 and "/" in varmap_item[varmap_index["AA_CHANGE"]]:
            aa_mut = varmap_item[varmap_index["AA_CHANGE"]].split("/")[1]
        else:
            aa_mut = varmap_item[varmap_index["AA_CHANGE"]]
        disease_status = ""
        disease_comments = ""
        user_id = varmap_item[varmap_index["USER_ID"]]

        # Is the variant disease causing?
        # This only applies to clinvar
        for list in clinvar_summary:
            for summary_line in list: #No idea why this [0] is needed. A list in a list should be what it is, not a list in a list in a list.
                # print(summary_line)
                    #print("Is", int(i[-1]), "equal to", int(USER_ID), "?" )
                try:
                    if int(summary_line[-1]) == int(user_id):  #  (variant id is last column in summary)
                        # print("clinvar summary and snipclip finally found a hit for variant ",int(var_record_id))

                        disease_status = disease_class(summary_line[6])
                        disease_comments = summary_line[24]
                except(ValueError):
                    print(user_id, "should be a number. Summary line:", summary_line)

        var_to_database(uniprot_record, var_record_location, aa_wt, aa_mut,
                        disease_status, disease_comments, variant_source, user_id)


def varmap_header_dict(varmap_file):
    with open(varmap_file) as inputfile:
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
    input_queries = inputs[0]
    input_query_set = inputs[1]

    ### VarMap ###
    # VarMap files can be big. Preprocessing them saves a lot of time

    varmap_files = {
        "clinvar": "scripts/external_datasets/clinvar_vartmh28_06_2019.tsv",
        "gnomad": "scripts/external_datasets/gnomAD_varsite.tsv"
    }

    #print(stripped_variant_list(varmap_files["clinvar"], input_query_set))
    #print(stripped_variant_list(varmap_files["gnomad"], input_query_set))

    clinvar_for_processing = stripped_variant_list(
        varmap_files["clinvar"], input_query_set)
    gnomad_for_processing = stripped_variant_list(
        varmap_files["gnomad"], input_query_set)

    clinvar_results_list = clinvar_for_processing[0]
    clinvar_results_set = clinvar_for_processing[1]

    clinvar_summary_lines = []
    print("Loading the variant summaries from ClinVar. This holds information on disease states in clinvar.")
    with open("scripts/external_datasets/variant_summary_25_06_2019.txt") as inputfile:
        for line_number, summary_variant in enumerate(inputfile):
            if line_number > 0:
                summary_variant = summary_variant.strip().split('\t')
                # print(str(summary_variant[-1]))
                # This relies on the last column being the id column
                if clean_query(str(summary_variant[-1])) in clinvar_results_set:
                    clinvar_summary_lines.append(summary_variant)
            else:
                pass


    gnomad_results_list = gnomad_for_processing[0]
    gnomad_results_set = gnomad_for_processing[1]

    varmap_process(clinvar_results_list, "ClinVar", varmap_header_dict(varmap_files["clinvar"]), clinvar_summary_lines)

    varmap_process(gnomad_results_list, "gnomAD", varmap_header_dict(varmap_files["gnomad"]))

    ### Humsavar ###
    humsavar(input_query_set)
