from __future__ import division
import requests
import urllib
from requests import get
import shutil
import numpy as np
import os
import time
import gzip
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
from django.db import models
from tmh_db.models import Database_Metadata, Subcellular_location, Flank, Flank_residue, Uniref, Go, Structure, Structural_residue, Funfam_residue, Funfamstatus, Protein, Residue, Tmh, Tmh_deltag, Tmh_hydrophobicity, Tmh_residue, Tmh_tmsoc, Variant, Keyword, Binding_residue
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
flank_length=5


def uniprot_bin(query_id):
    try:
        filename = str(f"scripts/external_datasets/uniprot_bin/{query_id}.txt")
        file = open(filename, "r")
        file_test = file.readlines
    # If the file is not found, an attempt is made to grab the file from the internet.
    except(FileNotFoundError):
        print("File not found:", filename)
        uniprot_url = str(f'https://www.uniprot.org/uniprot/{query_id}.txt')
        uniprot_bin = str(f"scripts/external_datasets/uniprot_bin/{query_id}.txt")
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

    tmh_input([query_id])

    binding_residues_to_table(filename)
    subcellular_location(filename)

    for keyword in record.annotations["keywords"]:
        keyword_to_database(keyword, query_id)

    record = SwissProt.read(open(filename))
    cross_reference = (record.cross_references)
    for go_mapping in cross_reference:
        if go_mapping[0] == 'GO':
            go_to_database(go_mapping[1], query_id)


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
                subcellular_location_for_database, created = Subcellular_location.objects.get_or_create(
                    location=f.qualifiers["description"].split(".")[0])
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
    avoid_features = ["TRANSMEM", "INTRAMEM"]

    # We iterate through each record, parsed by biopython.
    # First we need a list of all the TMH stop start positions
    # Now we can go through the record and write each TMH to the database (separate function)
    for record in SeqIO.parse(filename, input_format):
        total_tmh_number = total_tmh_uniprot(record)
        list_of_tmhs = []
        full_sequence = str(record.seq)
        new_record = True
        tmd_count = 0
        for i, f in enumerate(record.features):

            if f.type == feature_type:

                if "UnknownPosition" in str(f.location.start) or "UnknownPosition" in str(f.location.end):
                    # Let's count the tmhs even though they do not have a position and will not be in the database.
                    tmd_count = tmd_count + 1
                else:
                    # The variables need to be flushed
                    # query_id = None
                    tmh_start = None
                    tmh_stop = None
                    tmh_topology = None
                    #evidence_type = None
                    membrane_location = None
                    n_ter_seq = None
                    tmh_sequence = None
                    c_ter_seq = None
                    # tmh_number=None
                    tmd_count = tmd_count + 1
                    tmh_number = tmd_count
                    n_terminal_start = "none"

                    if "Helical" in str(f.qualifiers["description"]):
                        tm_type = "Helix"
                    elif "Strand" in str(f.qualifiers["description"]):
                        tm_type = "Beta"
                    else:
                        tm_type = "Unknown"

                    full_sequence = str(record.seq)
                    # This is the human readable value. It should not be used for slices
                    tmh_start = int(f.location.start) + 1
                    tmh_stop = int(f.location.end)
                    # Slices should not use +1 on start.
                    tmh_sequence = str(
                        record.seq[(f.location.start):(f.location.end)])

                    ### N terminal clash ###
                    # print(list_of_tmhs)

                    #Checks that the TMHs don't get mashed at the begining or end.
                    if int(f.location.start) - flank_length <= 0:
                        n_ter_seq = str(record.seq[0:(f.location.start)])
                    elif int(f.location.start) - flank_length > 0:
                        n_ter_seq = str(record.seq[(f.location.start - flank_length):(f.location.start)])

                    if int(f.location.end) + flank_length <= len(record.seq) :
                        c_ter_seq = str(record.seq[(f.location.end):(f.location.end + flank_length)])
                    elif int(f.location.end) + flank_length > len(record.seq):
                        c_ter_seq = str(record.seq[(f.location.end):(len(record.seq))])

                    # A list of common locations. These need sorting into inside/outside locations
                    io_dictionary_odd_even = uniprot_topo_check(record)

                    tmh_topology = io_dictionary_odd_even[odd_or_even(tmh_number)]

                    membrane_location = uniprot_membrane_location(record)
                    tmh_info=[query_id, tmh_number, total_tmh_number, tmh_start, tmh_stop, tmh_topology, evidence_type, membrane_location, n_ter_seq, tmh_sequence, c_ter_seq, evidence_type, full_sequence, tm_type]
                    #print(tmh_info)
                    tmh_list.append(tmh_info)
        tmh_list=integrity_check(tmh_list)
        tmh_list=clash_correction(tmh_list)
        tmh_to_database(tmh_list)
        return(tmh_list)


def total_tmh_uniprot(record):
    '''
    How many TRANSMEM regions in the TMP. This can include Beta.
    '''
    total_tmh = 0
    for i, f in enumerate(record.features):
        if f.type == "TRANSMEM":
            total_tmh = total_tmh + 1
    return(total_tmh)


def odd_or_even(number):
    '''
    Is a number odd or even? Takes a number and returns "even" or "odd".
    '''
    if number % 2 == 0:
        return("even")
    elif number % 2 != 0:
        return("odd")


def uni_subcellular_location(feature_description):
    '''
    Returns inside or outside from a record feature description in UniProt. This will need to be tweaked depending on what subcelullar organelles are present in Uniprot.
    '''

    locations = ["Chloroplast intermembrane", "Cytoplasmic", "Extracellular", "Intravacuolar", "Intravirion", "Lumenal", "Lumenal, thylakoid", "Lumenal, vesicle", "Mitochondrial intermembrane",
                 "Mitochondrial matrix", "Periplasmic", "Peroxisomal", "Peroxisomal matrix", "Nuclear", "Perinuclear space", "Stromal", "Vacuolar", "Vesicular", "Virion surface"]

    inside_locations = set(["Cytoplasmic", "Mitochondrial matrix", "Nuclear"])
    outside_locations = set(
        ["Extracellular", "Lumenal", "Periplasmic", "Mitochondrial intermembrane", "Perinuclear space", "Vesicular"])
    for i in inside_locations:
        if i in feature_description:
            return("Inside")
    for i in outside_locations:
        if i in feature_description:
            return("Outside")


def io_flip(io_value):
    '''
    Inverts inside to outside and vice versa. Also takes "Unknown" as a variable and returns unknown.
    '''

    if io_value == "Inside":
        return("Outside")
    elif io_value == "Outside":
        return("Inside")
    elif io_value == "Unknown":
        return("Unknown")


def odd_even_io(ordered_topo_list):
    '''
    Returns if odd or even has n-terminal-inside from an ordered topology list in the format:
    [['Outside', 0], ['TM', 51], ['Inside', 76], ['TM', 88], ['Outside', 114], ['TM', 124]]
    '''
    topo_only_list = []
    for n, i in enumerate(ordered_topo_list):
        if 'TM' in i:
            pass
        elif 'TM' in i and n == 0:
            topo_only_list.append(io_flip(ordered_topo_list[n + 1][0]))
        else:
            topo_only_list.append(i[0])
    #print(topo_only_list)

    if len(list(topo_only_list)) > 1:
        if 'Outside' in str(topo_only_list[0]):
            io_dict = {"even": "Inside",
                       "odd": "Outside"}
        elif 'Inside' in str(topo_only_list[0]):
            io_dict = {"even": "Outside",
                       "odd": "Inside"}
        elif topo_only_list[0] is None:
            if 'Inside' in str(topo_only_list[1]):
                io_dict = {"even": "Inside",
                           "odd": "Outside"}
            elif 'Outside' in str(topo_only_list[1]):
                io_dict = {"even": "Outside",
                           "odd": "Inside"}

    else:
        io_dict = {"even": "Unknown",
                   "odd": "Unknown"}
    return(io_dict)


def uniprot_membrane_location(record):
    '''
    Gets the TOPO_DOMs from a uniprot txt file string and returns a human
    readable list of locations
    '''

    topology_list = []
    for i, f in enumerate(record.features):
        if f.type == "TOPO_DOM":
            topology_list.append(f.qualifiers["description"].split(".")[0])
    locations = list(dict.fromkeys(topology_list))

    return(clean_query(str(locations)))


def location_to_number(exact_position):
    if "UnknownPosition" in str(exact_position):
        return(0)
    else:
        return(int(exact_position))


def uniprot_topo_check(record):

    topology_list = []
    # This function currently only deals with alpha helix
    for i, f in enumerate(record.features):
        if f.type == "TRANSMEM" and "Helix" in f.qualifiers["description"]:
            topology_list.append(["TM", location_to_number(f.location.start)])
        elif f.type == "TOPO_DOM":
            topology_list.append([uni_subcellular_location(
                f.qualifiers["description"].split(".")[0]), location_to_number(f.location.start)])
        else:
            pass

    ordered_list = sorted(topology_list, key=lambda x: x[1])

    for n, i in enumerate(ordered_list):
        if i[0] != "TM":
            io = i[0]
        elif i[0] == "TM":
            # If the next value is a TM, we insert a pseudo topo_dom. This is the result of short loops.
            if ordered_list[n + 1][0] == "TM":
                topology_list.insert([io_flip(io)])
    #print(ordered_list)
    return(odd_even_io(ordered_list))

def integrity_check(tmh_list):
    corrected_tmh_list=tmh_list
    for ref_tmh_order_number, ref_tmh_info in enumerate(tmh_list):
        if ref_tmh_info[2] != len(tmh_list):
            print("Disrepency between tmh number and number of TMHs retreived.")
            corrected_tmh_list[ref_tmh_order_number][2]=len(tmh_list)
        for comp_tmh_order_number, comp_tmh_info in enumerate(tmh_list):
            if ref_tmh_order_number == comp_tmh_order_number:
                pass
            else:
                if ref_tmh_info[1] > comp_tmh_info[1] and ref_tmh_info[3] < comp_tmh_info[3]:
                    print("Missmatch in TMH order. Check manually.")

    return(corrected_tmh_list)


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
                                                            tmh_start = int(tmh_details["Begin"])
                                                            tmh_stop = int(tmh_details["End"])
                                                            n_ter_seq = sequence[tmh_start - flank_length:tmh_start]
                                                            tmh_sequence = sequence[tmh_start:tmh_stop]
                                                            c_ter_seq = sequence[tmh_stop:tmh_stop + flank_length]

                                                            # preparing any non established variables for standard tmh recording.
                                                            full_sequence = sequence
                                                            tm_type = "Helix"
                                                            tmh_list.append([query_id, tmh_number, total_tmh_number, tmh_start + 1, tmh_stop, tmh_topology,
                                                                             evidence_type, membrane_location, n_ter_seq, tmh_sequence, c_ter_seq, evidence_type, full_sequence, tm_type])
                                                        # Although it is about as elegant as a sledgehammer,
                                                        # this catches the previous non tmh environment.
                                                        tmh_topology = tmh_details["Loc"]
                                    tmh_list=integrity_check(tmh_list)
                                    tmh_list=clash_correction(tmh_list)
                                    tmh_to_database(tmh_list)
                                    return(tmh_list)


def Sort(sub_li):
    '''
    # Python code to sort the tuples using second element
    # of sublist Inplace way to sort using sort()
    '''

    # reverse = None (Sorts in Ascending order)
    # key is set to sort using second element of
    # sublist lambda has been used
    sub_li.sort(key = lambda x: x[1])
    return(sub_li)

def clash_correction(tmh_list):
    '''
    Takes a series of tmh_info lists (see below) and cuts any flanks that might overlap.
    The list is returned in the original format.
    0           1           2                   3       4           5           6               7                   8       9               10          11          12              13
    [query_id, tmh_number, total_tmh_number, tmh_start, tmh_stop, tmh_topology,evidence_type, membrane_location, n_ter_seq, tmh_sequence, c_ter_seq, evidence_type, full_sequence, tm_type]
    '''
    #Sort by TMH info incase any XML stuff comes back in the wrong order.
    sorted_tmh_list=Sort(tmh_list)
    print("SORTED",len(sorted_tmh_list),sorted_tmh_list)
    correct_tmh_list=[]
    for ref_tmh_number, ref_tmh_info in enumerate(sorted_tmh_list):
        for comp_tmh_number, comp_tmh_info in enumerate(sorted_tmh_list):
            # Except single pass!
            if ref_tmh_info[2] == 1 and len(sorted_tmh_list)==1:
                correct_n_flank = ref_tmh_info[8]
                correct_c_flank = ref_tmh_info[10]

            elif ref_tmh_number == comp_tmh_number:
                pass
            else:
                # reference TMH n flank is less than the c flank end residue of a comparison TMH with a lower TMH number
                if int(ref_tmh_info[3]-len(ref_tmh_info[8])) < int(comp_tmh_info[4]+len(comp_tmh_info[10])) and ref_tmh_number>comp_tmh_number:
                    correct_n_flank=ref_tmh_info[8][len(ref_tmh_info[8])-int(abs(ref_tmh_info[3]-comp_tmh_info[4])/2)::]
                else:
                    correct_n_flank=ref_tmh_info[8]

                if int(ref_tmh_info[4]+len(ref_tmh_info[10])) > int(comp_tmh_info[3]-len(comp_tmh_info[8])) and ref_tmh_number<comp_tmh_number:
                    #The C terminal flank is greater than the the comparison TMH
                    correct_c_flank=ref_tmh_info[10][0:int(abs(ref_tmh_info[4]-comp_tmh_info[3])/2)]
                else:
                    correct_c_flank=ref_tmh_info[10]

        correct_tmh_info=ref_tmh_info
        correct_tmh_info[8]=correct_n_flank
        correct_tmh_info[10]=correct_c_flank
        correct_tmh_list.append(correct_tmh_info)

    return(correct_tmh_list)

def keyword_to_database(keyword, uniprot_id):
    print("Mapping keyword to", uniprot_id, ":", keyword)
    keyword_for_database, created = Keyword.objects.get_or_create(
        keyword=keyword)
    target_protein = Protein.objects.get(uniprot_id=uniprot_id)
    keyword_for_database.proteins.add(target_protein)


def add_n_flank(tmh_unique_id, n_ter_seq, tmh_topology, current_tmh):

    # Add the N terminal to the database
    current_tmh = Tmh.objects.get(tmh_id=tmh_unique_id)

    record_for_database, created = Flank.objects.update_or_create(
        tmh=current_tmh,
        n_or_c="N",
        defaults={
            "flank_sequence": n_ter_seq,
            "inside_or_outside": tmh_topology[0]
        }
    )


def add_c_flank(tmh_unique_id, c_ter_seq, tmh_topology, current_tmh):
    # Add the C terminal to the database
    if tmh_topology == "Inside":
        c_terminal_inside = "Outside"
    elif tmh_topology == "Outside":
        c_terminal_inside = "Inside"
    else:
        c_terminal_inside = tmh_topology
        print("N terminal was not inside or outside in", tmh_unique_id,
              "... setting C inside to whatever N inside is:", c_terminal_inside)

    record_for_database, created = Flank.objects.update_or_create(
        tmh=current_tmh,
        n_or_c="C",
        defaults={
            "flank_sequence": c_ter_seq,
            "inside_or_outside": c_terminal_inside[0]
        }
    )


def add_a_tmh_to_database(query_id, tmh_number, tmh_total_number, tmh_start, tmh_stop, tmh_topology, evidence_type, membrane_location, n_ter_seq, tmh_sequence, c_ter_seq, evidence, full_sequence, tm_type):

    tmh_protein = Protein.objects.get(uniprot_id=query_id)
    print("Generating tmh id key from\nQuery id:", query_id,
          "\nTMH number:", tmh_number, "\nEvidence:", evidence)
    tmh_unique_id = str(query_id + "." + str(tmh_number) + "." + evidence)

    print(tmh_unique_id)

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
            "n_terminal_inside": tmh_topology,
            "tm_type": tm_type
        }
    )

    # Now we run the calculations on anything that works at the TM level.
    tmsoc(tmh_unique_id, full_sequence, tmh_sequence, tmh_start, tmh_stop)
    deltag(tmh_unique_id, tmh_sequence)
    hydrophobicity(tmh_unique_id, full_sequence,
                   tmh_sequence, tmh_start, tmh_stop)

    transmembrane_helix = Tmh.objects.get(tmh_id=tmh_unique_id)
    add_n_flank(tmh_unique_id, n_ter_seq, tmh_topology, transmembrane_helix)
    add_c_flank(tmh_unique_id, c_ter_seq, tmh_topology, transmembrane_helix)

    # Now we will add residues to the TM residues table

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
            if tmh_topology == "Inside":
                feature_location = "Inside flank"
            elif tmh_topology == "Outside":
                feature_location = "Outside flank"
            else:
                feature_location = "Unknown"
        elif sequence_position > tmh_stop:
            #"C flank"
            if tmh_topology == "Inside":
                feature_location = "Outside flank"
            elif tmh_topology == "Outside":
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

        if "Inside" in tmh_topology:
            amino_acid_location_in_to_out = amino_acid_location_n_to_c
        elif "Outside" in tmh_topology:
            amino_acid_location_in_to_out = 0 - amino_acid_location_n_to_c
        elif "Unknown" in tmh_topology:
            amino_acid_location_in_to_out = None
        #print(tmh_topology)

        # specific_residue = Residue.objects.get(
        #    unique_together=[tmh_protein, sequence_position])
        print("Adding TM residue at position", sequence_position,
              "in ",  query_id, "from", tmh_unique_id)

        specific_residue = Residue.objects.get(
            protein=tmh_protein, sequence_position=int(sequence_position))

        if feature_location == "TMH":
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

        elif feature_location == "Inside flank":
            flank_n_or_c = None
            if tmh_topology in "Inside":
                flank_n_or_c = "N"
            elif tmh_topology in "Outside":
                flank_n_or_c = "C"
            this_flank = Flank.objects.get(
                tmh=transmembrane_helix, n_or_c=flank_n_or_c)
            record_for_database, created = Flank_residue.objects.update_or_create(
                residue=specific_residue,
                flank=this_flank,
                defaults={
                    "amino_acid_type": a_residue,
                    "evidence": evidence,
                    "amino_acid_location_n_to_c": amino_acid_location_n_to_c,
                    "amino_acid_location_in_to_out": amino_acid_location_in_to_out,
                    # inside flank, outside flank. inside flank, outside flank are ONLY flanking TMHs.
                    "feature_location": feature_location
                })

        elif feature_location == "Outside flank":
            flank_n_or_c = None
            if tmh_topology == "Inside":
                flank_n_or_c = "C"
            elif tmh_topology == "Outside":
                flank_n_or_c = "N"
            print(transmembrane_helix, flank_n_or_c)
            this_flank = Flank.objects.get(
                tmh=transmembrane_helix, n_or_c=flank_n_or_c)
            record_for_database, created = Flank_residue.objects.update_or_create(
                residue=specific_residue,
                flank=this_flank,
                defaults={
                    "amino_acid_type": a_residue,
                    "evidence": evidence,
                    "amino_acid_location_n_to_c": amino_acid_location_n_to_c,
                    "amino_acid_location_in_to_out": amino_acid_location_in_to_out,
                    # inside flank, outside flank. inside flank, outside flank are ONLY flanking TMHs.
                    "feature_location": feature_location
                })


def tmh_to_database(tmh_list):
    '''
    This takes the input standardised by the other functions of a TMH and adds them to the database.
    The function also has some integrity checks.
    '''
    # Now we have a complete list of the TMHs.

    for tmh_number_iteration, a_tmh in enumerate(tmh_list):
        print("TMH info:", a_tmh)
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
        evidence = a_tmh[11].replace("\n", "")
        full_sequence = a_tmh[12].replace("\n", "")
        tm_type = a_tmh[13].replace("\n", "")

        add_a_tmh_to_database(query_id, tmh_number, tmh_total_number, tmh_start, tmh_stop, tmh_topology,
                              evidence_type, membrane_location, n_ter_seq, tmh_sequence, c_ter_seq, evidence, full_sequence, tm_type)


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
            #sequence = result[0]
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
    try:
        deltag_result = deltag_result[8].split(" ")
        deltag_result = deltag_result[1]
    except(IndexError):
        pass
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
    topdb_url = "http://topdb.enzim.hu/?m=download&file=topdb_all.xml"
    topdb_file = 'scripts/external_datasets/topdb_all.xml'
    try:
        topdb = ET.parse(topdb_file)
    except FileNotFoundError:
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
        print("Checking tmhs in...")
        print("UniProt...")
        uniprot_tm_check(a_query)
        print("Checking tmhs in...")
        print("TopDB...")
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
        print("Adding UniProt record", a_query, " to table,", query_number + 1, "of", len(input_query), "records...")
        uniprot_table(a_query)

    ### TMH Tables ###
    # tmh_input(input_query)
