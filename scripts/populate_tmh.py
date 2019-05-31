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
from tmh_db.models import Database_Metadata, Subcellular_location, Uniref, Go, Structure, Structural_residue, Funfam_residue, Funfamstatus, Protein, Residue, Tmh, Tmh_deltag, Tmh_hydrophobicity, Tmh_residue, Tmh_tmsoc, Variant, Keyword, Binding_residue
from datetime import datetime, timedelta
from django.utils import timezone
from datetime import date
import pytz
from scripts.populate_general_functions import *

print("Usage:\npython manage.py runscript populate --traceback")

# How many days should be allowed to not enforce updates
time_threshold = 7
today = date.today()
todaysdate = today.strftime("%d_%m_%Y")


def uniprot_bin(query_id):
    try:
        filename = str(f"scripts/external_datasets/uniprot_bin/{query_id}.txt")
        file = open(filename, "r")
        file_test = file.readlines
    # If the file is not found, an attempt is made to grab the file from the internet.
    except(FileNotFoundError):
        print("File not found:", filename)
        uniprot_url = str(f'https://www.uniprot.org/uniprot/{query_id}.txt')
        uniprot_bin = str(
            f"scripts/external_datasets/uniprot_bin/{query_id}.txt")
        download(uniprot_url, uniprot_bin)


def uniprot_table(query_id):
    filename = str(f"scripts/external_datasets/uniprot_bin/{query_id}.txt")
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

    binding_residues_to_table(filename)
    subcellular_location(filename)

    for keyword in record.annotations["keywords"]:
        keyword_to_database(keyword, query_id)

    record = SwissProt.read(open(filename))
    cross_reference = (record.cross_references)
    for map in cross_reference:
        if map[0] == 'GO':
            go_to_database(map[1], query_id)


def binding_residues_to_table(filename):
    for record in SeqIO.parse(filename, "swiss"):
        protein = Protein.objects.get(uniprot_id=record.id)
        for i, f in enumerate(record.features):
            if f.type == "BINDING":
                for position in range(f.location.start, f.location.end):

                    specific_residue = Residue.objects.get(
                        protein=protein, sequence_position=int(position))
                    record_for_database, created = Binding_residue.objects.update_or_create(
                        residue=specific_residue,
                        comment=f.qualifiers,)

def subcellular_location(filename):
    for record in SeqIO.parse(filename, "swiss"):
        protein = Protein.objects.get(uniprot_id=record.id)
        for i, f in enumerate(record.features):
            if f.type == "TOPO_DOM":
                subcellular_location_for_database, created = Subcellular_location.objects.get_or_create(location=f.qualifiers["description"])
                subcellular_location_for_database.proteins.add(protein)


def go_to_database(go_id, uniprot_id):
    print("Mapping GO", go_id, "to", uniprot_id)
    go_for_database, created = Go.objects.get_or_create(go_id=go_id)
    target_protein = Protein.objects.get(uniprot_id=uniprot_id)
    go_for_database.proteins.add(target_protein)

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

    filename = str(f"scripts/external_datasets/uniprot_bin/{query_id}.txt")
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
                                uniprot_ref_sequence = Residue.objects.filter(
                                    protein__uniprot_id=query_id).count()
                                if len(sequence) == uniprot_ref_sequence:

                                    add_topdb = True
                                elif len(sequence) != uniprot_ref_sequence:
                                    print("TOPD Uniprot length mismatch in",
                                          query_id, ". Aborting TOPDB TMH record")
                                    add_topdb = False
                                # print(sequence)
                                if add_topdb is True:
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




def keyword_to_database(keyword, uniprot_id):
    print("Mapping keyword to", uniprot_id, ":", keyword)
    keyword_for_database, created = Keyword.objects.get_or_create(
        keyword=keyword)
    target_protein = Protein.objects.get(uniprot_id=uniprot_id)
    keyword_for_database.proteins.add(target_protein)


def topology_tidy(tmh_list):

    # Is there any information on what odd/even tmh numbers should topologically be in terms of subcellular location and inside/outside?
    topo_odd = "Unknown"
    topo_even = "Unknown"
    location_odd = "Unknown"
    location_even = "Unknown"

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
                    (len(sequences_to_add) - 1 / 2) + len(n_ter_seq)

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


def tmh_input(input_query):
    print("Extracting TMH bounadries from...")
    # Parse the xml static files since this is the slowest part.
    # Ignore this for now -  we need to sort out uniprot before anything else!
    topdb_url="http://topdb.enzim.hu/?m=download&file=topdb_all.xml"
    topdb_file = 'scripts/external_datasets/topdb_all.xml'
    download(topdb_url, topdb_file)
    topdb = ET.parse(topdb_file)
    # mptopo = ET.parse('mptopoTblXml.xml')
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
