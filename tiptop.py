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
    output = str(a_list).replace("'", "")
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


#    for item in root.findall("item"):
#        ElementTree.dump(item)

    for node in mptopo.findall('.//mptopoProtein'):
        record_present = False
        features = node.getchildren()
        for feature in features:

            if str(feature.tag) == str("uniprotNumber"):
                # print("Found a uniprot id")
                # print(str(feature.text))
                if str(feature.text) == str(query_id):
                    # print("Matches query...")
                    for tm_find in features:
                        print(str(tm_find))
                        if str(tm_find.tag) == str("tmSegments"):
                            record_present = True
                            # print("Checking for tmsegments")
                            tmhs = tm_find.getchildren()
                            # print(tmhs)
                            for tmh_segment in tmhs:

                                tmh_locations = tmh_segment.getchildren()

                                for tmh_location in tmh_locations:
                                    if str(tmh_location.tag) == str("beginIndex"):
                                        tmh_start = int(tmh_location.text)
                                    elif str(tmh_location.tag) == str("endIndex"):
                                        tmh_stop = int(tmh_location.text)
                                    else:
                                        pass
                                print_list(
                                    [query_id, tmh_start, tmh_stop, evidence_type])
    return(record_present)

# Check topdb


def topdb_check(query_id):
    '''
    Checks the TOPDB xml file for transmem regions mapped to the UniProt search ID.
    '''

    record_present = False

    evidence_type = str("TOPDB")


#    for item in root.findall("item"):
#        ElementTree.dump(item)

    for node in topdb.findall('.//TOPDB'):
        records = node.getchildren()
        for features in records:
            if str(features.tag) == str("CrossRef"):
                id_types = features.getchildren()
                for id_type in id_types:
                    if id_type.tag == str("UniProt"):
                        acs = id_type.getchildren()
                        for ids in acs:
                            if str(ids.text) == query_id:
                                for feature in records:
                                    if str(feature.tag) == str("Topology"):
                                        topology = feature.getchildren()
                                        for item in topology:
                                            if str(item.tag) == str("Regions"):
                                                tmhs = item.getchildren()
                                                for tmh in tmhs:
                                                    tmh_details = tmh.attrib
                                                    if str(tmh_details["Loc"]) == str("Membrane"):
                                                        record_present = True
                                                        tmh_start = tmh_details["Begin"]
                                                        tmh_stop = tmh_details["End"]
                                                        print_list(
                                                            [query_id, tmh_start, tmh_stop, evidence_type])

    return(record_present)

# Check UniProt function


def uniprot_check(query_id):
    '''
    This fetches the uniprot id from either a local bin or the internet and checks the annotation for TRANSMEM regions.
    '''
    record_present = False
    evidence_type = str("UniProt")

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
    #subcellular_location = "TOPO_DOM"
    #avoid_features = ["TRANSMEM", "INTRAMEM"]

    # We iterate through each record, parsed by biopython.
    for record in SeqIO.parse(filename, input_format):
        new_record = True
        tmd_count = 0
        for i, f in enumerate(record.features):
            if f.type == feature_type:
                record_present = True
                full_sequence = str(record.seq)
                tmh_start = int(f.location.start)
                tmh_stop = int(f.location.end)
                tmh_sequence = str(
                    record.seq[(f.location.start):(f.location.end)])
                print_list([query_id, tmh_start, tmh_stop, evidence_type])
    return(record_present)


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


# run the searches on each query in the input list
print("UniProt ID, TMH start position, TMH end position, Database source")
for a_query in input_query:
    a_query = clean_query(a_query)
    # print(clean_query(a_query))
    mptopo_check(a_query)
    uniprot_check(a_query)
    topdb_check(a_query)
