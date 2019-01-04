from __future__ import division
import requests
from Bio import SeqIO
import numpy as np
import os
import subprocess
import re
import sys
import xml.etree.ElementTree as ET

# Print clean list


def print_list(a_list):
    '''
    Prints a human readable csv line from a python list.
    '''
    # print(a_list)
    if a_list is not None and str(a_list) != str("[]"):
        output = str(a_list).replace("], [", "\n")
        output = str(output).replace("'", "")
        output = str(output).replace("[", "")
        output = str(output).replace("]", "")
        print(output)

### Database query functions ###

# Check MPtopo


def mptopo_check(query_id):
    '''
    Checks the MPTOPO xml file for transmemembrane regions mapped to a UniProt ID.
    '''

    # Sequences don't exactly match UniProt
    evidence_type = str("MPTopo")

    # for item in root.findall("item"):
    #   ElementTree.dump(item)

    for node in mptopo.findall('.//mptopoProtein'):
        features = node.getchildren()
        for feature in features:

            if str(feature.tag) == str("uniprotNumber") and str(feature.text) == str(query_id):
                # print("Matches query...")
                for tm_find in features:
                    # print(str(tm_find))
                    if str(tm_find.tag) == str("nTerminal"):
                        # Frustratingly, the database only includes the first topology. Re-entrant helices will therefor be incorrect.
                        starting_topology = tm_find.text
                        # print(str(starting_topology))
                    if str(tm_find.tag) == str("tmSegments"):
                        # print("Checking for tmsegments")
                        tmhs = tm_find.getchildren()
                        # print(tmhs)

                        tmh_list = []
                        for tm_number, tmh_segment in enumerate(tmhs):
                            tmh_locations = tmh_segment.getchildren()
                            for tmh_location in tmh_locations:
                                if str(tmh_location.tag) == str("beginIndex"):
                                    tmh_start = int(tmh_location.text)
                                elif str(tmh_location.tag) == str("endIndex"):
                                    tmh_stop = int(tmh_location.text)
                                else:
                                    pass

                            # get around 0 base counting
                            tmh_number = tm_number + 1

                            # %2==0 checks if number is even.
                            # if tmh number is even and N terminal is inside
                            if tmh_number % 2 == 0 and str(starting_topology) == str("in"):
                                tmh_topology = "Outside"
                            # if tmh number is even and N terminal is outside
                            elif tmh_number % 2 == 0 and str(starting_topology) == str("out"):
                                tmh_topology = "Inside"
                            # if tmh number is odd and N terminal is inside
                            elif tmh_number % 2 != 0 and str(starting_topology) == str("in"):
                                tmh_topology = "Inside"
                            # if tmh number is odd and N terminal is inside
                            elif tmh_number % 2 != 0 and str(starting_topology) == str("out"):
                                tmh_topology = "Outside"
                            else:
                                tmh_topology = "None"

                            tmh_list.append(
                                [query_id, tmh_start, tmh_stop, tmh_topology, evidence_type])
                        return(tmh_list)

# Check TOPDB

def topdb_check(query_id):
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

                                                    if str(tmh_details["Loc"]) == str("Membrane"):
                                                        tmh_start = int(
                                                            tmh_details["Begin"])
                                                        tmh_stop = int(
                                                            tmh_details["End"])
                                                        ### NO FLANK CLASH CHECKS! NUMBERS WILL BE WRONG!!! ###
                                                        n_ter_seq = sequence[tmh_start -
                                                                             5:tmh_start]
                                                        tmh_sequence = sequence[tmh_start:tmh_stop]
                                                        c_ter_seq = sequence[tmh_stop:tmh_stop + 5]

                                                        tmh_list.append(
                                                            [query_id, tmh_start, tmh_stop, tmh_topology, evidence_type, membrane_location, n_ter_seq, tmh_sequence, c_ter_seq])

                                                    # Although it is about as elegant as a sledgehammer,
                                                    # this catches the previous non tmh environment.
                                                    tmh_topology = tmh_details["Loc"]
                                return(tmh_list)
        sequence = None
        membrane_location = None

# Check UniProt function

def uniprot_check(query_id):
    '''
    This fetches the uniprot id from either a local bin or the internet and
    checks the annotation for TRANSMEM regions.
    '''

    evidence_type = str("UniProt")
    tmh_list = []

    # The UniProt bin contains lots of uniprot files.
    try:
        filename = str("uniprot_bin/" + query_id + ".txt")
        file = open(filename, "r")
        file.readlines
    # If the file is not found, an attempt is made to grab the file from the internet.
    except(FileNotFoundError):
        uniprot_url = str(
            'https://www.uniprot.org/uniprot/%s.txt' % (query_id))
        r = requests.get(uniprot_url)

        with open(str("uniprot_bin/" + query_id + ".txt"), 'wb') as f:
            f.write(r.content)

    # These are the parameters used by the biopython Seq.IO module

    filename = str("uniprot_bin/" + query_id + ".txt")
    input_format = "swiss"
    feature_type = "TRANSMEM"
    subcellular_location = "TOPO_DOM"
    #subcellular_location = "TOPO_DOM"
    #avoid_features = ["TRANSMEM", "INTRAMEM"]

    # We iterate through each record, parsed by biopython.
    # First we need a list of all the TMH stop start positions
    for record in SeqIO.parse(filename, input_format):
        list_of_tmhs = []
        for i, f in enumerate(record.features):
            if f.type == feature_type:
                list_of_tmhs.append(int(f.location.start))
                list_of_tmhs.append(int(f.location.end))

    # Now we can go through the record and write each TMH and any info to a file (or just print it!)
    for record in SeqIO.parse(filename, input_format):
        new_record = True
        tmd_count = 0
        for i, f in enumerate(record.features):
            if f.type == feature_type:
                n_terminal_start = "none"
                record_present = True
                full_sequence = str(record.seq)
                tmh_start = int(f.location.start)
                tmh_stop = int(f.location.end)
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
                if int(f.location.start) - 5 <= 0 and n_clash == False:
                    n_ter_seq = str(record.seq[0:(f.location.start)])
                elif int(f.location.start) - 5 > 0 and n_clash == False:
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

                if int(f.location.end) + 5 <= len(record.seq) and c_clash == False:
                    c_ter_seq = str(
                        record.seq[(f.location.end):(f.location.end + 5)])
                elif int(f.location.end) + 5 > len(record.seq) and c_clash == False:
                    c_ter_seq = str(
                        record.seq[(f.location.end):(len(record.seq))])

                # A list of common locations. These need sorting into inside/outside locations
                locations = ["Chloroplast intermembrane", "Cytoplasmic", "Extracellular", "Intravacuolar", "Intravirion", "Lumenal", "Lumenal, thylakoid", "Lumenal, vesicle", "Mitochondrial intermembrane",
                             "Mitochondrial matrix", "Periplasmic", "Peroxisomal", "Peroxisomal matrix", "Nuclear", "Perinuclear space", "Stromal", "Vacuolar", "Vesicular", "Virion surface"]
                if n_terminal_start == "none" and tmh_start > 1:
                    previous_feautre_location = tmh_start - 1
                    for index, a_features in enumerate(record.features):
                        tmh_topology = None
                        membrane_location = None
                        if 'UnknownPosition' in str(a_features.location.start) or 'UnknownPosition' in str(a_features.location.end):
                            pass
                        else:
                            if a_features.type == subcellular_location and a_features.location.start < previous_feautre_location and a_features.location.end > previous_feautre_location:
                                inside_locations = [
                                    "Cytoplasmic", "Mitochondrial matrix"]
                                outside_locations = [
                                    "Extracellular", "Lumenal", "Periplasmic", "Mitochondrial intermembrane"]
                                for location in inside_locations:
                                    if location in str(a_features.qualifiers):
                                        tmh_topology = "Inside"
                                        membrane_location = location
                                for location in outside_locations:
                                    if location in str(a_features.qualifiers):
                                        tmh_topology = "Outside"
                                        membrane_location = location

                                tmh_list.append([query_id, tmh_start, tmh_stop, tmh_topology,
                                                 evidence_type, membrane_location, n_ter_seq, tmh_sequence, c_ter_seq])
        return(tmh_list)


def clean_query(query):
    '''
    This aims to generate a clean ascii query of a viable UniProt ID.
    '''
    illegal_characters = ["!", "\n", " ", "@"]
    for char in illegal_characters:
        query = query.replace(char, "")
    a_clean_query = query
    #print("Clean query result:", a_clean_query)
    return(a_clean_query)


### Canonical script starts here ###

# Parse the xml static files since this is the slowest part
topdb = ET.parse('topdb_all.xml')
mptopo = ET.parse('mptopoTblXml.xml')

# Grab the input list
user_file = str(sys.argv[1])
input_file = open(user_file, 'r')
input_query = input_file.readlines()


# Run the searches on each query in the input list
print("UniProt ID, TMH start position, TMH end position, N-terminal starting side, Database source, N-terminal starting subcellular location, n-flank-seq, tmh-seq, c-flank-seq")
for a_query in input_query:
    a_query = clean_query(a_query)
    # print(clean_query(a_query))
    ### OPM needs adding to here also. ###
    # print_list(mptopo_check(a_query))
    print_list(uniprot_check(a_query))
    print_list(topdb_check(a_query))
