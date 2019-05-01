from __future__ import division
import requests
from requests import get
import shutil
import numpy as np
import os
import time
import subprocess
import json
from subprocess import check_output
import re
import sys
import defusedxml.ElementTree as ET
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
# env LDFLAGS="-I/usr/local/opt/openssl/include -L/usr/local/opt/openssl/lib" pip install psycopg2
from django.conf import settings
from django.db import models
from tmh_db.models import Database_Metadata, Funfam_residue, Funfamstatus, Protein, Residue, Tmh, Tmh_deltag, Tmh_hydrophobicity, Tmh_residue, Tmh_tmsoc, Variant, Keyword
from datetime import datetime, timedelta
from django.utils import timezone
from datetime import date
import pytz

# Complex query example. How many variants are in the TMH?
# >>> Variant.objects.exclude(residue__tmh_residue=None).filter(disease_status='d').count()
# 1153
# >>> Variant.objects.exclude(residue__tmh_residue=None).filter(disease_status='d', variant_source="ClinVar").count()
# 394
# >>> Variant.objects.exclude(residue__tmh_residue=None).filter(disease_status='d', variant_source="Humsavar").count()
# 759


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


def uniprot_bin(query_id):
    try:
        filename = str(
            "scripts/external_datasets/uniprot_bin/" + query_id + ".txt")
        file = open(filename, "r")
        file_test = file.readlines
    # If the file is not found, an attempt is made to grab the file from the internet.
    except(FileNotFoundError):
        print("File not found.", query_id, ".")
        uniprot_url = str(f'https://www.uniprot.org/uniprot/{query_id}.txt')
        uniprot_bin = str(
            "scripts/external_datasets/uniprot_bin/" + query_id + ".txt")
        download(uniprot_url, uniprot_bin)


def uniprot_table(query_id):
    filename = str("scripts/external_datasets/uniprot_bin/" +
                   query_id + ".txt")
    input_format = "swiss"
    feature_type = "TRANSMEM"
    tm_protein = False
    print("Checking UniProt for TM annotation in", query_id, ".")
    for record in SeqIO.parse(filename, input_format):


        list_of_tmhs = []
        tmh_count = 0
        # features locations is a bit annoying as the start location needs +1 to match the sequence IO, but end is the correct sequence value.
        for i, f in enumerate(record.features):
            if f.type == feature_type:
                if "UnknownPosition" in str(f.location.start) or "UnknownPosition" in str(f.location.end):
                    print(record.id, "Unknown position for TMH in record")
                    # This doesn't mean that there are coordinates, just that there is a TM somewhere in the protein.
                    tm_protein = True
                else:
                    list_of_tmhs.append(int(f.location.start) + 1)
                    list_of_tmhs.append(int(f.location.end))
                    tmh_count = tmh_count + 1
                    tm_protein = True
        sequence = record.seq

    record_for_database, created = Protein.objects.update_or_create(
        uniprot_id=query_id)

    print("TM annotation found in", query_id, ".")
    target_protein = Protein.objects.get(uniprot_id=query_id)

    record_for_database, created = Protein.objects.update_or_create(
        uniprot_id=query_id,
        defaults={
            "total_tmh_number": tmh_count,
            "full_sequence": str(sequence),
            "updated_date": timezone.now()
        }
    )

    residue_table(query_id, sequence)

    for keyword in record.annotations["keywords"]:
        keyword_to_database(keyword, query_id)
    #print(record.cross_references)
    #for map in record.annotations["features"]:
    #    go_to_database(go_id, query_id)

def residue_table(query_id, sequence):
    protein = Protein.objects.get(uniprot_id=query_id)
    # Now we build the residue table.
    # Are there ever skips of unknow length? This could affect TMH number.

    # This method will not update. A separate out of date script should be used to check if this needs to be removed and updated.
    existing = set(Residue.objects.values_list(
        "protein__uniprot_id", "sequence_position"))
    residues_to_create = []

    for residue_number, a_residue in enumerate(sequence):
        if not (query_id, residue_number + 1) in existing:
            residues_to_create.append(
                Residue(
                    protein=protein,
                    sequence_position=residue_number + 1,
                    amino_acid_type=a_residue,
                )
            )
    Residue.objects.bulk_create(residues_to_create)


def uniprot_tm_check(query_id):
    '''
    This fetches the uniprot id from either a local bin or the internet and
    checks the annotation for TRANSMEM regions.
    '''
    print("Checking", query_id, "in UniProt.")
    evidence_type = str("UniProt")
    tmh_list = []

    # The UniProt bin contains lots of uniprot files.
    # This bin should either routinely be flushed, or have some sort of timestamp.

    # These are the parameters used by the biopython Seq.IO module

    filename = str("scripts/external_datasets/uniprot_bin/" +
                   query_id + ".txt")
    input_format = "swiss"
    feature_type = "TRANSMEM"
    subcellular_location = "TOPO_DOM"
    # subcellular_location = "TOPO_DOM"
    # avoid_features = ["TRANSMEM", "INTRAMEM"]

    # We iterate through each record, parsed by biopython.
    # First we need a list of all the TMH stop start positions
    for record in SeqIO.parse(filename, input_format):
        list_of_tmhs = []
        total_tmh_number = 0
        # features locations is a bit annoying as the start location needs +1 to match the sequence IO, but end is the correct sequence value.
        for i, f in enumerate(record.features):
            if f.type == feature_type:
                if "UnknownPosition" in str(f.location.start) or "UnknownPosition" in str(f.location.end):
                    print(record.id, "Unknown position for TMH in record")
                    # this stops unknown tmhs masking polytopic as single pass
                    total_tmh_number = total_tmh_number + 1
                else:
                    list_of_tmhs.append(int(f.location.start) + 1)
                    list_of_tmhs.append(int(f.location.end))
                    total_tmh_number = total_tmh_number + 1

        if list_of_tmhs:  # Checks if it is a TM protein.
            # The protein for the database
            # This feels like it could go somewhere else...
            print("TMHs found in", query_id)

        else:
            print("No TMHs found in UniProt text file for", query_id)

        # We can also check if any isoforms are in or near the TM region
        for i, f in enumerate(record.features):
            if f.type == "VAR_SEQ":
                for n, x in enumerate(record.features):
                    # Remember, feature type is set to transmembrane regions
                    if x.type == feature_type:
                        if int(x.location.start) + 1 - 5 <= int(f.location.end) and int(f.location.end) <= int(x.location.end) + 5:
                            # print("An isoform will interfere with record", record.id)
                            pass
                        elif int(x.location.start) + 1 - 5 < int(f.location.start) + 1 and int(f.location.start) + 1 <= int(x.location.end) + 5:
                            # print("An isoform will interfere with record", record.id)
                            pass

    # Now we can go through the record and write each TMH to the database (separate function)
    for record in SeqIO.parse(filename, input_format):
        full_sequence = str(record.seq)
        new_record = True
        tmd_count = 0
        for i, f in enumerate(record.features):

            if f.type == feature_type:

                if "UnknownPosition" in str(f.location.start) or "UnknownPosition" in str(f.location.end):
                    # Let's count the tmhs even though they do not have a position and will not be in the database.
                    tmd_count = tmd_count + 1
                else:

                    tmd_count = tmd_count + 1
                    tmh_number = tmd_count
                    n_terminal_start = "none"
                    record_present = True
                    full_sequence = str(record.seq)
                    # This is the human readable value. It should not be used for slices
                    tmh_start = int(f.location.start) + 1
                    tmh_stop = int(f.location.end)
                    # Slices should not use +1 on start.
                    tmh_sequence = str(
                        record.seq[(f.location.start):(f.location.end)])

                    ### N terminal clash ###
                    n_clash = False
                    # print(list_of_tmhs)
                    for position in list_of_tmhs:
                        # checks if another tmh/non-flank feature is near
                        # print(int(f.location.start - 5) , position , int(f.location.start))
                        if int(f.location.start - 10) <= position < int(f.location.start):
                            n_ter_seq = str(
                                record.seq[int(position + abs(position - int(f.location.start)) / 2):(f.location.start)])
                            n_clash = True
                            print("N-clash detected")

                    if int(f.location.start) - 5 <= 0 and n_clash is False:
                        n_ter_seq = str(record.seq[0:(f.location.start)])
                    elif int(f.location.start) - 5 > 0 and n_clash is False:
                        n_ter_seq = str(
                            record.seq[(f.location.start - 5):(f.location.start)])

                    ### C terminal clash ###
                    c_clash = False
                    for position in list_of_tmhs:
                        # checks if another tmh/non-flank feature is near
                        if int(f.location.end) < position <= int(f.location.end) + 10:
                            c_ter_seq = str(record.seq[int(f.location.end):int(
                                (abs(int(f.location.end) - position) / 2) + int(f.location.end))])
                            c_clash = True
                            print("C-clash detected")

                    if int(f.location.end) + 5 <= len(record.seq) and c_clash is False:
                        c_ter_seq = str(
                            record.seq[(f.location.end):(f.location.end + 5)])
                    elif int(f.location.end) + 5 > len(record.seq) and c_clash is False:
                        c_ter_seq = str(
                            record.seq[(f.location.end):(len(record.seq))])

                    # A list of common locations. These need sorting into inside/outside locations
                    locations = ["Chloroplast intermembrane", "Cytoplasmic", "Extracellular", "Intravacuolar", "Intravirion", "Lumenal", "Lumenal, thylakoid", "Lumenal, vesicle", "Mitochondrial intermembrane",
                                 "Mitochondrial matrix", "Periplasmic", "Peroxisomal", "Peroxisomal matrix", "Nuclear", "Perinuclear space", "Stromal", "Vacuolar", "Vesicular", "Virion surface"]
                    tmh_topology = "Unknown"
                    membrane_location = "Unknown"
                    if n_terminal_start == "none" and tmh_start > 1:
                        previous_feautre_location = tmh_start - 1

                        for index, a_features in enumerate(record.features):
                            tmh_topology = None
                            membrane_location = ''
                            if 'UnknownPosition' in str(a_features.location.start) or 'UnknownPosition' in str(a_features.location.end):
                                pass
                            else:

                                # Using topo_dom will only work if there are no short loops. Short loops could be assumed, or labelled as no topology.
                                if a_features.type == subcellular_location and a_features.location.start < previous_feautre_location and a_features.location.end > previous_feautre_location:
                                    inside_locations = [
                                        "Cytoplasmic", "Mitochondrial matrix"]
                                    outside_locations = [
                                        "Extracellular", "Lumenal", "Periplasmic", "Mitochondrial intermembrane"]
                                    print("echo", n_terminal_start, tmh_start)
                                    # Here we should set inside and outside, and then use that for even/odd with checks.

                                    for location in inside_locations:
                                        if location in str(a_features.qualifiers):
                                            tmh_topology = "Inside"
                                            membrane_location = location
                                    for location in outside_locations:
                                        if location in str(a_features.qualifiers):
                                            tmh_topology = "Outside"
                                            membrane_location = location

                    tmh_list.append([query_id, tmh_number, total_tmh_number, tmh_start, tmh_stop, tmh_topology,
                                     evidence_type, membrane_location, n_ter_seq, tmh_sequence, c_ter_seq, evidence_type, full_sequence])

        corrected_tmh_list = topology_tidy(tmh_list)
        tmh_to_database(corrected_tmh_list)

        return(tmh_list)


def topdb_check(query_id, topdb):
    '''
    Checks the TOPDB xml file for transmem regions mapped to the UniProt search ID.
    '''
    evidence_type = str("TOPDB")

#    for item in root.findall("item"):
#        ElementTree.dump(item)

    for node in topdb.findall('.//TOPDB'):
        # Clears the sequence in case of a blank or dodgy record.
        sequence = None
        membrane_location = None

        records = node.getchildren()
        for features in records:
            if str(features.tag) == str("Sequence"):
                for seq_info in features:
                    if str(seq_info.tag) == str("Seq"):
                        sequence = str(seq_info.text).replace("\n", "")
                        sequence = sequence.replace(" ", "")
                        # print(sequence)
        for features in records:
            if str(features.tag) == str("Membrane"):

                membrane_location = str(features.text)
                # print("Membrane known: ", membrane_location)
                # print(sequence)

        for features in records:
            if str(features.tag) == str("CrossRef"):
                id_types = features.getchildren()
                for id_type in id_types:
                    if id_type.tag == str("UniProt"):
                        acs = id_type.getchildren()
                        for ids in acs:
                            if str(ids.text) == query_id:
                                tmh_list = []
                                for feature in records:
                                    if str(feature.tag) == str("Topology"):
                                        topology = feature.getchildren()
                                        for item in topology:
                                            if str(item.tag) == str("Regions"):
                                                tmhs = item.getchildren()
                                                for feature_number, feature in enumerate(tmhs):
                                                    tmh_details = feature.attrib

                                                    # get total number of tmhs
                                                    total_tmh_number = 0
                                                    for a_feature_number, a_feature in enumerate(tmhs):
                                                        if str(tmh_details["Loc"]) == str("Membrane"):
                                                            total_tmh_number = total_tmh_number + 1
                                                    tmh_number = 0
                                                    if str(tmh_details["Loc"]) == str("Membrane"):
                                                        tmh_number = tmh_number + 1
                                                        tmh_start = int(
                                                            tmh_details["Begin"])
                                                        tmh_stop = int(
                                                            tmh_details["End"])
                                                        ### NO FLANK CLASH CHECKS! NUMBERS WILL BE WRONG!!! ###
                                                        n_ter_seq = sequence[tmh_start
                                                                             - 5:tmh_start]
                                                        tmh_sequence = sequence[tmh_start:tmh_stop]
                                                        c_ter_seq = sequence[tmh_stop:tmh_stop + 5]

                                                        # preparing any non established variables for standard tmh recording.
                                                        full_sequence = sequence

                                                        tmh_list.append([query_id, tmh_number, total_tmh_number, tmh_start + 1, tmh_stop, tmh_topology,
                                                                         evidence_type, membrane_location, n_ter_seq, tmh_sequence, c_ter_seq, evidence_type, full_sequence])
                                                    # Although it is about as elegant as a sledgehammer,
                                                    # this catches the previous non tmh environment.
                                                    tmh_topology = tmh_details["Loc"]
                                tmh_to_database(tmh_list)
                                return(tmh_list)
        sequence = None
        membrane_location = None


def go_to_database(go_id, uniprot_id):
    go_for_database, created = Go.objects.get_or_create(go_id=go_id)
    target_protein = Protein.objects.get(uniprot_id=uniprot_id)
    go_for_database.proteins.add(target_protein)


def keyword_to_database(keyword, uniprot_id):
    print("Mapping keyword to", uniprot_id, ":", keyword)
    keyword_for_database, created = Keyword.objects.get_or_create(keyword=keyword)
    target_protein = Protein.objects.get(uniprot_id=uniprot_id)
    keyword_for_database.proteins.add(target_protein)


def topology_tidy(tmh_list):

    # Is there any information on what odd/even tmh numbers should topologically be in terms of subcellular location and inside/outside?
    topo_odd="Unknown"
    topo_even="Unknown"
    location_odd="Unknown"
    location_even="Unknown"

    for tmh_number, tmh_properties in enumerate(tmh_list):
        topology = tmh_properties[5]
        location = tmh_properties[7]
        if tmh_number % 2 == 0:
            # Even
            if topology is not None:
                topo_even = topology
                location_even = location

        else:
            # Odd
            if topology is not None:
                topo_odd = topology
                location_odd = location

    for tmh_number, tmh_properties in enumerate(tmh_list):
        if tmh_number % 2 == 0:
            # Even
            tmh_list[tmh_number][5] = topo_even
            tmh_list[tmh_number][7] = location_even
        else:
            # Odd
            tmh_list[tmh_number][5] = topo_odd
            tmh_list[tmh_number][7] = location_odd
    return(tmh_list)


def tmh_to_database(tmh_list):
    '''
    This takes the input standardised by the other functions of a TMH and adds them to the database.
    The function also has some integrity checks.
    '''
    # Now we have a complete list of the TMHs.
    for tmh_number_iteration, a_tmh in enumerate(tmh_list):
        query_id = a_tmh[0]
        tmh_number = a_tmh[1]
        tmh_total_number = a_tmh[2]
        tmh_start = a_tmh[3]
        tmh_stop = a_tmh[4]
        tmh_topology = a_tmh[5]
        evidence_type = a_tmh[6]
        # At this point we have all the membrane locations, and some may be dead. Integrity needs to be run to ensure this makes sense, at least in an IO sense.
        membrane_location = a_tmh[7]
        n_ter_seq = a_tmh[8].replace("\n", "")
        tmh_sequence = a_tmh[9].replace("\n", "")
        c_ter_seq = a_tmh[10].replace("\n", "")
        evidence = a_tmh[11]
        full_sequence = a_tmh[12].replace("\n", "")

        if tmh_topology is None:
            tmh_topology = "Unknown"

        if tmh_topology in "Inside":
            n_terminal_inside = "Inside"
        elif tmh_topology in "Outside":
            n_terminal_inside = "Outside"
        else:
            n_terminal_inside = "Unknown"
        tmh_protein = Protein.objects.get(uniprot_id=query_id)
        tmh_unique_id = str(query_id + "." + str(tmh_number) + "." + evidence)

        print(tmh_unique_id)
        print(a_tmh)

        # The TMH for the database
        record_for_database, created = Tmh.objects.update_or_create(
            protein=tmh_protein,
            tmh_id=tmh_unique_id,
            defaults={
                "tmh_sequence": tmh_sequence,
                "tmh_start": tmh_start,
                "tmh_stop": tmh_stop,
                "tmh_evidence": evidence,
                "membrane_type": membrane_location,
                "tmh_number": tmh_number,
                "tmh_total_number": tmh_total_number,
                "n_terminal_inside": n_terminal_inside
            }
        )

        # Now we run the calculations on anything that works at the TM level.
        tmsoc(tmh_unique_id, full_sequence, tmh_sequence, tmh_start, tmh_stop)
        deltag(tmh_unique_id, tmh_sequence)
        hydrophobicity(tmh_unique_id, full_sequence,
                       tmh_sequence, tmh_start, tmh_stop)

        # Now we will add residues to the TM residues table
        transmembrane_helix = Tmh.objects.get(tmh_id=tmh_unique_id)
        # Now we build the residue table.
        # Are there ever skips of unknow length? This could affect TMH number.

        # This method will not update. A separate out of date script should be used to check if this needs to be removed and updated.

        sequences_to_add = str(n_ter_seq + tmh_sequence + c_ter_seq)

        for tmh_residue_number, a_residue in enumerate(sequences_to_add):

            # Where is the residue in reference to the TMH
            sequence_position = int(
                tmh_start - len(n_ter_seq)) + tmh_residue_number

            if sequence_position < tmh_start:  # This doesn't make sense
                #"N flank"
                if n_terminal_inside == "Inside":
                    feature_location = "Inside flank"
                elif n_terminal_inside == "Outside":
                    feature_location = "Outside flank"
                else:
                    feature_location = "Unknown"
            elif sequence_position > tmh_stop:
                #"C flank"
                if n_terminal_inside == "Inside":
                    feature_location = "Outside flank"
                elif n_terminal_inside == "Outside":
                    feature_location = "Inside flank"
                else:
                    feature_location = "Unknown"
            elif sequence_position >= tmh_start and sequence_position <= tmh_stop:
                feature_location = "TMH"

            # What is the exact residue position.

            # This avoids odd numbers rounding to 0 twice at -1 and +1.

            if len(sequences_to_add) % 2 == 0:
                amino_acid_location_n_to_c = tmh_residue_number - \
                    len(sequences_to_add) / 2 + \
                    len(n_ter_seq)  # These correction numbers are wrong!
            else:
                amino_acid_location_n_to_c = tmh_residue_number - \
                    (len(sequences_to_add) + 1 / 2) + len(n_ter_seq)

            if "Inside" in n_terminal_inside:
                amino_acid_location_in_to_out = amino_acid_location_n_to_c
            elif "Outside" in n_terminal_inside:
                amino_acid_location_in_to_out = 0 - amino_acid_location_n_to_c
            elif "Unknown" in n_terminal_inside:
                amino_acid_location_in_to_out = None

            # specific_residue = Residue.objects.get(
            #    unique_together=[tmh_protein, sequence_position])
            print("Adding TM residue at position", sequence_position,
                  "in ",  query_id, "from", tmh_unique_id)

            specific_residue = Residue.objects.get(
                protein=tmh_protein, sequence_position=int(sequence_position))

            record_for_database, created = Tmh_residue.objects.update_or_create(
                residue=specific_residue,
                tmh_id=transmembrane_helix,
                defaults={
                    "amino_acid_type": a_residue,
                    "evidence": evidence,
                    "amino_acid_location_n_to_c": amino_acid_location_n_to_c,
                    "amino_acid_location_in_to_out": amino_acid_location_in_to_out,
                    # This is either TMH or flank. TMH, inside flank, outside flank.
                    "feature_location":  feature_location
                }
            )


def tmsoc(tmh_unique_id, full_sequence, tmh_sequence, tmh_start, tmh_stop):

        # Writing input files for TMSOC
    with open("scripts/external_scripts/tmsoc/inputseq.fasta", 'w') as temp_tmh_fasta:
        temp_tmh_fasta.write(str(">" + tmh_unique_id + "\n"))
        sequence_counter = 0

        # Break needed every 60 characters.
        for aa in full_sequence:
            if sequence_counter < 59:
                temp_tmh_fasta.write(aa)
                sequence_counter = sequence_counter + 1
            elif sequence_counter >= 59:
                temp_tmh_fasta.write(str(aa + "\n"))
                sequence_counter = 0

    # Coordinates from the original TMSOC project look like
    # 37,63 74,96
    with open("scripts/external_scripts/tmsoc/TMsegments.txt", 'w') as temp_tmh_seq:
        temp_tmh_seq.write(str(tmh_start) + "," + str(tmh_stop) + " ")

    # Running TMSOC

    tmsoc_result = check_output(["/usr/bin/perl", "scripts/external_scripts/tmsoc/TMSOC.pl", "scripts/external_scripts/tmsoc/inputseq.fasta",
                                 "scripts/external_scripts/tmsoc/TMsegments.txt"])  # stdout=subprocess.PIPE)
    tmsoc_result = tmsoc_result.decode("utf-8")
    tmsoc_result = tmsoc_result.split("\n")
    for line in tmsoc_result:
        if str(tmh_sequence + ";") in line:
            result = line.split(";")
            sequence = result[0]
            start_stop = result[1]
            score_one = result[2]
            score_two = result[3]
            z_score = result[4]
            simple_twighlight_complex = result[5]

            print("Writing TMSOC results to database: ",
                  z_score, simple_twighlight_complex)

            # Add to database
            # Fetch the tmh foreign key
            tmh_protein = Tmh.objects.get(tmh_id=tmh_unique_id)
            # What are we recording
            # Record to the database
            record_for_database, created = Tmh_tmsoc.objects.update_or_create(
                tmh=tmh_protein,
                defaults={
                    "test_type": "TMSOC",
                    "test_result": simple_twighlight_complex,
                    "test_score": float(z_score)
                }
            )


def deltag(tmh_unique_id, tmh_sequence):
    with open("scripts/external_scripts/dgpred/inputseq.fasta", 'w') as temp_tmh_fasta:
        temp_tmh_fasta.write(str(">" + tmh_unique_id + "\n" + tmh_sequence))

    deltag_result = check_output(["/usr/bin/perl", "scripts/external_scripts/dgpred/myscanDG.pl",
                                  "scripts/external_scripts/dgpred/inputseq.fasta"])  # stdout=subprocess.PIPE)
    print(deltag_result)
    deltag_result = deltag_result.decode("utf-8")
    deltag_result = deltag_result.split("\n")
    deltag_result = deltag_result[8].split(" ")
    deltag_result = deltag_result[1]
    print("deltag=", deltag_result)

    tmh_protein = Tmh.objects.get(tmh_id=tmh_unique_id)
    # What are we recording
    # Record to the database
    record_for_database, created = Tmh_deltag.objects.update_or_create(
        tmh=tmh_protein,
        defaults={
            "test_type": "Delta G",
            "test_score": float(deltag_result)
        }
    )


def window_slice(list_for_slicing, window_length, start_slice, end_slice, full_sequence_length):
    '''
    This checks that the slice needed does not exceed the final residue.
    '''
    # print(list_for_slicing, window_length,
    #      start_slice, end_slice, full_sequence_length)
    if int(window_length / 2) + start_slice >= full_sequence_length:
        windowed_values = list_for_slicing[start_slice +
                                           int(window_length / 2) - 1:]

    elif int(window_length / 2) + end_slice >= full_sequence_length:
        windowed_values = list_for_slicing[int(
            len(list_for_slicing) - window_length / 2) - 1:]

    elif int(window_length / 2) + end_slice < full_sequence_length:
        windowed_values = list_for_slicing[int(
            start_slice + window_length / 2) - 1:int(end_slice + window_length / 2) - 1]
    # print(windowed_values) # These values should be printed to a graph and stored in a file.
    return(windowed_values)


def hydrophobicity(tmh_unique_id, full_sequence, tmh_sequence, tmh_start, tmh_stop):
    window_length = 5
    edge = 1

    tmh_sequence_analysis = ProteinAnalysis(str(tmh_sequence))
    print(len(tmh_sequence))
    full_sequence_analysis = ProteinAnalysis(str(full_sequence))

    aromaticity = tmh_sequence_analysis.aromaticity()
    print("Aromaticity:", aromaticity)
    flexibility = np.mean(tmh_sequence_analysis.flexibility())
    print("Flexibility:", flexibility)

    ww = {'A': 0.33, 'R': 1.00, 'N': 0.43, 'D': 2.41, 'C': 0.22, 'Q': 0.19, 'E': 1.61, 'G': 1.14, 'H': -0.06, 'I': -0.81,
          'L': -0.69, 'K': 1.81, 'M': -0.44, 'F': -0.58, 'P': -0.31, 'S': 0.33, 'T': 0.11, 'W': -0.24, 'Y': 0.23, 'V': -0.53}
    ww_window = full_sequence_analysis.protein_scale(ww, window_length, edge)
    ww_window = window_slice(ww_window, window_length,
                             tmh_start, tmh_stop, len(full_sequence))
    ww_avg = np.mean(ww_window)
    print("White Wimley:", ww_avg)

    kyte = {'A': 1.8, 'R': -4.5, 'N': -3.5, 'D': -3.5, 'C': 2.5, 'Q': -3.5, 'E': -3.5, 'G': -0.4, 'H': -3.2, 'I': 4.5,
            'L': 3.8, 'K': -3.9, 'M': 1.9, 'F': 2.8, 'P': -1.6, 'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'V': 4.2}
    kyte_window = full_sequence_analysis.protein_scale(
        kyte, window_length, edge)
    kyte_window = window_slice(
        kyte_window, window_length, tmh_start, tmh_stop, len(full_sequence))
    kyte_avg = np.mean(kyte_window)
    print("Kyte:", kyte_avg)

    # These numbers need chaning
    eisenberg = {'A': 0.620, 'R': -2.530, 'N': -0.780, 'D': -0.900, 'C': 0.290, 'Q': -0.850, 'E': -0.740, 'G': 0.480, 'H': -0.400,
                 'I': 1.380, 'L': 1.060, 'K': -1.500, 'M': 0.640, 'F': 1.190, 'P': 0.120, 'S': -0.180, 'T': -0.050, 'W': 0.810, 'Y': 0.260, 'V': 1.080}
    eisenberg_window = full_sequence_analysis.protein_scale(
        eisenberg, window_length, edge)
    eisenberg_window = window_slice(
        eisenberg_window, window_length, tmh_start, tmh_stop, len(full_sequence))
    eisenberg_avg = np.mean(eisenberg_window)
    print("Eisenberg:", eisenberg_avg)
    # if len(kyte_window)== len(tmh_sequence):
    #    print("YAY!")
    # else:
    #    print("DOH!")
    tmh_protein = Tmh.objects.get(tmh_id=tmh_unique_id)
    # What are we recording
    # Record to the database
    record_for_database, created = Tmh_hydrophobicity.objects.update_or_create(
        tmh=tmh_protein,
        defaults={
            "aromaticity": aromaticity,
            "flexibility": flexibility,
            "kyte_avg": kyte_avg,
            "ww_avg": ww_avg,
            "eisenberg_avg": eisenberg_avg,
            "kyte_window": str(kyte_window),
            "ww_window": str(ww_window),
            "eisenberg_window": str(eisenberg_window),
        }
    )


def get_uniprot():
    '''
    Downloads UniProt IDs from Human transmembrane proteins from UniProt.org.
    '''
    # Grab the input list
    print("Fetching UniProt TM protein IDs")
    uniprot_list_url = "https://www.uniprot.org/uniprot/?query=reviewed%3Ayes+AND+organism%3A%22Homo+sapiens+%28Human%29+%5B9606%5D%22+AND+annotation%3A%28type%3Atransmem%29&sort=score&columns=id,&format=tab"
    # "https://www.uniprot.org/uniprot/?query=reviewed%3Ayes+AND+annotation%3A(type%3Atransmem)&sort=score&columns=id,&format=tab"
    # uniprot_list = 'https://www.uniprot.org/uniprot/?query=reviewed%3Ayes+AND+organism%3A"Homo+sapiens+(Human)+[9606]"+AND+annotation%3A(type%3Atransmem)&sort=score&columns=id,&format=tab'
    uniprot_list_file = "scripts/external_datasets/uniprot_bin/uniprot_list" + \
        todaysdate + ".txt"

    download(uniprot_list_url, uniprot_list_file)

    # This saves the request to a file for reasons beyond me.
    # So we now need to open the file to recover the items as a list.
    with open(uniprot_list_file) as f:
        # Somehow this has already made a list.
        lines = f

        input_query = list(lines)
        # Entry is the first line, which will break later code as it is not a valid uniprot id!
        input_query = input_query[1:]

    return(input_query)


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


def input_query_process(input_query):
    input_queries = []
    for query_number, a_query in enumerate(input_query):
        a_query = clean_query(a_query)
        print("Checking cache/downloading", a_query, ",",
              query_number + 1, "of", len(input_query), "records...")

        input_queries.append(a_query)

    input_query_set = set(input_queries)
    return([input_queries, input_query_set])


def funfam_submit(a_query):
    print(a_query)
    protein = Protein.objects.get(uniprot_id=a_query)
    record_for_database, created = Funfamstatus.objects.update_or_create(
        protein=protein,
        defaults={
        }
    )

    # check our database for submission key to funfam
    funfam_key = Funfamstatus.objects.get(protein=protein).submission_key
    funfam_completed_date = Funfamstatus.objects.get(
        protein=protein).completed_date
    print("Funfam key for query", a_query, ":", funfam_key)


    # Convert the UniProt binned file to a fasta.
    fasta_file = f"./scripts/external_datasets/fasta_bin/{a_query}.fasta"
    # print(fasta_file)
    with open(str(f"./scripts/external_datasets/uniprot_bin/{a_query}.txt"), "rU") as input_handle:
        with open(str(fasta_file), "w") as output_handle:
            sequences = SeqIO.parse(input_handle, "swiss")
            count = SeqIO.write(sequences, output_handle, "fasta")

        # submit the query to CATH funfam
        base_url = 'http://www.cathdb.info/search/by_funfhmmer'

        with open(fasta_file) as x:
            fasta_contents = x.read()
            data = {'fasta': fasta_contents, "queue": "hmmscan_funvar"}
            headers = {'accept': 'application/json'}
            if len(str(fasta_contents)) > 2000:  # Quick FIX!!!
                pass
            else:

                r = requests.post(base_url, data=data, headers=headers)
                funfam_submission_code = r.json()
                print(funfam_submission_code)
                funfam_key = funfam_submission_code['task_id']

                print("submitted task: " + funfam_key)

            record_for_database, created = Funfamstatus.objects.update_or_create(
                protein=protein,
                defaults={
                    "submission_key": funfam_key,
                }
            )
        return(funfam_key)


def funfam_result(a_query, funfam_submission_code):
    base_url = f'http://www.cathdb.info/search/by_funfhmmer/check/{funfam_submission_code}'

    headers = {'accept': 'application/json'}

    r = requests.get(base_url, headers=headers)
    # Result is something like this: {'success': 0, 'data': {'date_completed': '', 'status': 'queued', 'worker_hostname': '', 'id': '715b00cba220424897cb09df7e81129f', 'date_started': ''}, 'message': 'queued'}
    funfam_status = r.json()
    print(funfam_status)

    # Results can take a while to complete. Best to just add those that have finished. A week in and everything should have settled down.
    if funfam_status['success'] == 1:
        headers = {'accept': 'application/json'}
        results_url = f'http://www.cathdb.info/search/by_funfhmmer/results/{funfam_submission_code}'
        r = requests.get(results_url, headers=headers)
        if str(r) == "<Response [204]>":
            print("No funfam hits for ", a_query)
        elif str(r) == "<Response [200]>":
            funfam_api_result = r.json()
            protein = Protein.objects.get(uniprot_id=a_query)
            funfam_to_update = Funfamstatus.objects.get(protein=protein)
            # record_for_database, created = Funfamstatus.objects.update(
            #    protein=protein,
            #    completed_date=timezone.now(),

            #    #"funfam":

            # )
            funfam_to_update = Funfamstatus.objects.get(protein=protein)
            funfam_to_update.completed_date = timezone.now()  # change field
            funfam_to_update.save()  # this will update only

            # This next bit isn't perfect. We assign each residue in the region the funfam score and id. There may be a way to elegantly put in another table, but given the queries we are asking, this will suffice.
            for key, value in funfam_api_result.items():
                print("Key:", key, "Value:", value)
                print("\n")


def sifts_mapping(a_query):

    # Download via the API
    print("Fetching sifts information")
    sifts_file = str(
        "scripts/external_datasets/sifts_mapping/" + a_query + ".json")
    sifts_url = f"https://www.ebi.ac.uk/pdbe/api/mappings/all_isoforms/{a_query}"
    download(sifts_url, sifts_file)
    # PARSE and add to database
    data = "data"


def tmh_input(input_query):
    print("Extracting TMH bounadries from...")
    for query_number, a_query in enumerate(input_query):
        a_query = clean_query(a_query)
        print("\nExtracting TMH boundaries for", a_query, ",",
              query_number + 1, "of", len(input_query), "records.")
        # print(clean_query(a_query))
        ### OPM needs adding to here also. ###
        # mptopo_tm_check(a_query)
        uniprot_tm_check(a_query)
        topdb_check(a_query, topdb)


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
    input_query = get_uniprot()
    # Here we will just use a watered down list of tricky proteins. Uncomment this line for testing the whole list.
    #input_query = ["P22760", "Q5K4L6", "Q7Z5H4", "P32897", "Q9NR77", "P31644", "Q9NS61"]

    # Parse the xml static files since this is the slowest part.
    # Ignore this for now -  we need to sort out uniprot before anything else!
    topdb = ET.parse('scripts/external_datasets/topdb_all.xml')
    # mptopo = ET.parse('mptopoTblXml.xml')

    # Also, parse the variant files which can be massive.
    # humsavar table
    print(input_query)
    print("Starting TMH database population script...")
    os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'tmh_database.settings')

    ### Download uniprot files ###
    inputs = input_query_process(input_query)
    input_queries = inputs[0]
    input_query_set = inputs[1]

    ### If UniProt says it is a TMH, add it to the protein table ###

    for query_number, a_query in enumerate(input_query):
        print("Checking UniProt bin for", a_query)
        a_query = clean_query(a_query)
        uniprot_bin(a_query)
        print("Adding UniProt record", a_query, " to table,",
              query_number + 1, "of", len(input_query), "records...")
        uniprot_table(a_query)


    ### TMH Tables ###
    tmh_input(input_query)


    # Now get all TM information from these and build a residue table flat file.
    # Check residues in TMHs are consistent. Ditch anything that does not match uniprot and print to log.
    # now generate flat files for VarSite.

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


    ### Redundancy tables ###

    # The funfams need to be submitted, then checked for status and results.
    # This submits all the ids to the funfams and gets job ids.
    print("Finding closest funfams for records...")
    uniprotid_funfam_dict = {}
    for a_query in input_query:
        a_query = clean_query(a_query)
        this_funfam = funfam_submit(a_query)
        uniprotid_funfam_dict.update({a_query: this_funfam})

    # This uses the job id to wait until the job is complete and fetch the result.
    for a_query in input_query:
        a_query = clean_query(a_query)
        funfam = funfam_result(a_query, uniprotid_funfam_dict[a_query])

    # Populate structures
    for a_query in input_query:
        a_query = clean_query(a_query)
        sifts_mapping(a_query)


    print("This is the end of the script. It seems like there were no script breaking errors along the way.")
