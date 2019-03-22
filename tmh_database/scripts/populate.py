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
from tmh_db.models import Tmh
from tmh_db.models import Residue
from tmh_db.models import Tmh_residue
from tmh_db.models import Variant

print("Usage:\npython3 manage.py runscript populate --traceback")


def uniprot_bin(query_id):
    try:
        "Already in cache."
        filename = str(
            "scripts/external_datasets/uniprot_bin/" + query_id + ".txt")
        file = open(filename, "r")
        file.readlines
    # If the file is not found, an attempt is made to grab the file from the internet.
    except(FileNotFoundError):
        uniprot_url = str(f'https://www.uniprot.org/uniprot/{query_id}.txt')
        r = requests.get(uniprot_url)

        with open(str("scripts/external_datasets/uniprot_bin/" + query_id + ".txt"), 'wb') as f:
            f.write(r.content)


def uniprot_table(query_id):
    filename = str("scripts/external_datasets/uniprot_bin/"
                   + query_id + ".txt")
    input_format = "swiss"
    feature_type = "TRANSMEM"
    tm_protein = False
    for record in SeqIO.parse(filename, input_format):
        list_of_tmhs = []
        # features locations is a bit annoying as the start location needs +1 to match the sequence IO, but end is the correct sequence value.
        for i, f in enumerate(record.features):
            if f.type == feature_type:
                if "UnknownPosition" in str(f.location.start) or "UnknownPosition" in str(f.location.end):
                    pass
                    print(record.id, "Unknown position for TMH in record")
                    # This doesn't mean that there are coordinates, just that there is a TM somewhere in the protein.
                    tm_protein = True
                else:
                    list_of_tmhs.append(int(f.location.start) + 1)
                    list_of_tmhs.append(int(f.location.end))
                    tm_protein = True
        sequence = record.seq
    # Add the protein to the protein table if it is a TMP
    if tm_protein == True:
        record_for_database = Protein(
            uniprot_id=query_id, full_sequence=str(sequence))
        record_for_database, created = Protein.objects.update_or_create(
            uniprot_id=query_id,
            defaults={
                "full_sequence": str(sequence),
            }
        )


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

    filename = str("scripts/external_datasets/uniprot_bin/"
                   + query_id + ".txt")
    input_format = "swiss"
    feature_type = "TRANSMEM"
    subcellular_location = "TOPO_DOM"
    #subcellular_location = "TOPO_DOM"
    #avoid_features = ["TRANSMEM", "INTRAMEM"]

    # We iterate through each record, parsed by biopython.
    # First we need a list of all the TMH stop start positions
    for record in SeqIO.parse(filename, input_format):
        list_of_tmhs = []
        total_tmh_number = 0
        # features locations is a bit annoying as the start location needs +1 to match the sequence IO, but end is the correct sequence value.
        for i, f in enumerate(record.features):
            if f.type == feature_type:
                if "UnknownPosition" in str(f.location.start) or "UnknownPosition" in str(f.location.end):
                    pass
                    print(record.id, "Unknown position for TMH in record")
                    # this stops unknown tmhs masking polytopic as single pass
                    total_tmh_number = total_tmh_number + 1
                else:
                    list_of_tmhs.append(int(f.location.start) + 1)
                    list_of_tmhs.append(int(f.location.end))
                    total_tmh_number = total_tmh_number + 1

        if len(list_of_tmhs) > 0:  # Checks if it is a TM protein.
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

        new_record = True
        tmd_count = 0
        for i, f in enumerate(record.features):

            if f.type == feature_type:

                if "UnknownPosition" in str(f.location.start) or "UnknownPosition" in str(f.location.end):
                    # Let's count the tmhs even though they do not have a position and will not be in the database.
                    tmd_count = tmd_count + 1
                    pass
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
                            print("C-clash detected")

                    if int(f.location.end) + 5 <= len(record.seq) and c_clash == False:
                        c_ter_seq = str(
                            record.seq[(f.location.end):(f.location.end + 5)])
                    elif int(f.location.end) + 5 > len(record.seq) and c_clash == False:
                        c_ter_seq = str(
                            record.seq[(f.location.end):(len(record.seq))])

                    # A list of common locations. These need sorting into inside/outside locations
                    locations = ["Chloroplast intermembrane", "Cytoplasmic", "Extracellular", "Intravacuolar", "Intravirion", "Lumenal", "Lumenal, thylakoid", "Lumenal, vesicle", "Mitochondrial intermembrane",
                                 "Mitochondrial matrix", "Periplasmic", "Peroxisomal", "Peroxisomal matrix", "Nuclear", "Perinuclear space", "Stromal", "Vacuolar", "Vesicular", "Virion surface"]
                    tmh_topology="Unknown"
                    membrane_location="Unknown"
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
                                     evidence_type, membrane_location, n_ter_seq, tmh_sequence, c_ter_seq, evidence_type])

        tmh_to_database(tmh_list)

        return(tmh_list)


def tmh_to_database(tmh_list):
    '''
    This takes the input standardised by the other functions of a TMH and adds them to the database.
    The function also has some integrity checks.
    '''
    print(tmh_list)
    # Now we have a complete list of the TMHs.
    for tmh_number_iteration, a_tmh in enumerate(tmh_list):
        print(a_tmh)
        query_id = a_tmh[0]
        tmh_number = a_tmh[1]
        tmh_total_number = a_tmh[2]
        tmh_start = a_tmh[3]
        tmh_stop = a_tmh[4]
        tmh_topology = a_tmh[5]
        evidence_type = a_tmh[6]
        # At this point we have all the membrane locations, and some may be dead. Integrity needs to be run to ensure this makes sense, at least in an IO sense.
        membrane_location = a_tmh[7]
        n_ter_seq = a_tmh[8]
        tmh_sequence = a_tmh[9]
        c_ter_seq = a_tmh[10]
        evidence = a_tmh[11]

        if tmh_topology == None:
            tmh_topology="Unknown"

        if tmh_topology in "Inside":
            n_terminal_inside="Inside"
        elif tmh_topology in "Outside":
            n_terminal_inside="Outside"
        else:
            n_terminal_inside="Unknown"
        tmh_protein = Protein.objects.get(uniprot_id=query_id)
        tmh_unique_id = str(query_id + "." + str(tmh_number) + "." + evidence)
        print(tmh_unique_id)
        # The TMH for the database
        record_for_database = Tmh(protein=tmh_protein, tmh_total_number=tmh_total_number, tmh_id=tmh_unique_id, tmh_sequence=tmh_sequence,
                                  tmh_start=tmh_start, tmh_stop=tmh_stop, tmh_evidence=evidence,  membrane_type=membrane_location, tmh_number=tmh_number,n_terminal_inside=n_terminal_inside)
        record_for_database, created = Tmh.objects.update_or_create(
            protein=tmh_protein,
            tmh_id= tmh_unique_id,
            defaults={
                "tmh_sequence": tmh_sequence,
                "tmh_start": tmh_start,
                "tmh_stop": tmh_stop,
                "tmh_evidence": evidence,
                "membrane_type": membrane_location,
                "tmh_number": tmh_number,
                "tmh_total_number": tmh_total_number,
                "n_terminal_inside":n_terminal_inside
            }
        )


def get_uniprot():
    '''
    Downloads UniProt IDs from Human transmembrane proteins from UniProt.org.
    '''
    # Grab the input list
    print("Fetching UniProt TM protein IDs")
    uniprot_list = "https://www.uniprot.org/uniprot/?query=reviewed%3Ayes+AND+organism%3A%22Homo+sapiens+%28Human%29+%5B9606%5D%22+AND+annotation%3A%28type%3Atransmem%29&sort=score&columns=id,&format=tab"
    #"https://www.uniprot.org/uniprot/?query=reviewed%3Ayes+AND+annotation%3A(type%3Atransmem)&sort=score&columns=id,&format=tab"
    #uniprot_list = 'https://www.uniprot.org/uniprot/?query=reviewed%3Ayes+AND+organism%3A"Homo+sapiens+(Human)+[9606]"+AND+annotation%3A(type%3Atransmem)&sort=score&columns=id,&format=tab'

    uniprot_request = urllib.request.urlretrieve(str(uniprot_list))
    # This saves the request to a file for reasons beyond me.
    # So we now need to open the file to recover the items as a list.
    with open(uniprot_request[0]) as f:
        # Somehow this has already made a list.
        lines = f
        #lines = f.read().splitlines()

        input_query = list(lines)
        input_query = input_query[1:] # Entry is the first line, which will break later code as it is not a valid uniprot id!
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
    #print("Clean query result:", a_clean_query)
    return(a_clean_query)


def run():
    '''
    This is what django runs. This is effectively the canonical script,
    even though django forces it to be in a function.
    $ python3 manage.py runscript populate --traceback
    '''

    ### Canonical script starts here ###

    # In full scale mode it will take forever. Here we will just use a watered down list of tricky proteins.
    input_query=get_uniprot()
    #input_query = ["P32897", "Q9NR77", "P31644", "Q96E22", "P47869", "P28472", "P18507", "P05187"]

    # Parse the xml static files since this is the slowest part.
    # Ignore this for now -  we need to sort out uniprot before anything else!
    #topdb = ET.parse('topdb_all.xml')
    #mptopo = ET.parse('mptopoTblXml.xml')

    # Also, parse the variant files which can be massive.
    # humsavar table

    print(input_query)
    print("Starting TMH database population script...")
    os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'tmh_database.settings')

    ### Download uniprot files ###

    for query_number, a_query in enumerate(input_query):
        a_query = clean_query(a_query)
        print("Checking cache/downloading", a_query, ",",
              query_number + 1, "of", len(input_query), "records...")
        uniprot_bin(a_query)

    ### If UniProt says it is a TMH, add it to the protein table ###

    for query_number, a_query in enumerate(input_query):
        a_query = clean_query(a_query)
        uniprot_table(a_query)

    ### TMH Tables ###

    print("Extracting TMH bounadries from...")
    for query_number, a_query in enumerate(input_query):
        a_query = clean_query(a_query)
        print("Extracting TMH boundaries for", a_query, ",",
              query_number + 1, "of", len(input_query), "records.")
        # print(clean_query(a_query))
        ### OPM needs adding to here also. ###
        # mptopo_tm_check(a_query)
        uniprot_tm_check(a_query)
        # topdb_check(a_query)

    # Now get all TM information from these and build a residue table flat file.
    # Check residues in TMHs are consistent. Ditch anything that does not match uniprot and print to log.
    # now generate flat files for VarSite.

    ### Redundancy tables ###
    #print("Finding closest funfams for records...")
    # for a_query in input_query:
    #    # Convert the UniProt binned file to a fasta.
    #    with open(str(f"./scripts/external_datasets/uniprot_bin/{a_query}.txt"), "rU") as input_handle:
    #    with open(str(f"./scripts/external_datasets/fasta_bin/{a_query}.fasta"), "w") as output_handle:
    #        sequences = SeqIO.parse(input_handle, input_format)
    #        count = SeqIO.write(sequences, output_handle, "fasta")
    #    sequence = record.seq
    #    funfam = subprocess.Popen(["perl", "./scripts/external_scripts/cath-tools-seqscan-master/script/cath-tools-seqscan.pl", str(f"--in=./scripts/external_datasets/fasta_bin/{a_query}.fasta"), "--max_hits=1", "--max_aln=3"], stdout=subprocess.PIPE)
    #    print(funfam)

    ### Variant tables ###

    #print("Reading the variant tables...")
    #humsavar_list = []
    # with open('humsavar.txt') as f:
    #    lines = f.read().splitlines()
    #    for i in lines:
    #        i = i.replace('  ', ' ')
    #        humsavar_list.append(i.split())
    # print(humsavar_list)

    # Load the  varsite tsv file from snip clip.
    #clinvar_results = []
    # with open("./scripts/external_datasets/clinvar_snipclipa13_02_2019.tsv") as inputfile:
    #    for line_number, line in enumerate(inputfile):
    #        if line_number == 0:
    #            pass
    #        else:
    #            for a_query in input_query: # This can massively reduce memory use and time.
    #                if clean_query(a_query) in str(line):
    #                    clinvar_results.append(line.strip().split('\t'))
    # print(var_results)

    # Load the clinvar summary file
    #clinvar_summary_lines = []
    # with open("variant_summary.txt") as inputfile:
    #    for line_number, line in enumerate(inputfile):
    #        if line_number == 0:
    #            pass
    #        else:
    #            for variant in var_results:
    #                # print(str(variant[26])) # This should be simply the file ID. This filter kind of speeds it up if there are not many items.
    #                if clean_query(variant[26]) in str(line):
    #                    clinvar_summary_lines.append(line.strip().split('\t'))

    # populate()
