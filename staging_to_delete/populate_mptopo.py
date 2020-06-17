from __future__ import division

import gzip
import json
import os
import re
import shutil
import sys
import time
import urllib
from datetime import date
from datetime import datetime
from datetime import timedelta
from subprocess import check_output

import Bio
import defusedxml.ElementTree as ET
import numpy as np
import pytz
import requests
from Bio import SeqIO
from Bio import SwissProt
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from django.db import models
from django.utils import timezone
from requests import get

from scripts.populate_general_functions import *
from scripts.populate_tmh import clash_correction
from scripts.populate_tmh import integrity_check
from scripts.populate_tmh import tmh_to_database
from scripts.populate_tmh import add_a_tmh_to_database



def mptopo_file_check():

    mptopo_url="https://blanco.biomol.uci.edu/mpstruc/mptopo/mptopoAlphaHlxTblXml"
    mptopo_file="scripts/external_datasets/mptopo_alpha.xml"

    try:
        mptopo_file_parse= ET.parse(mptopo_file)
    except FileNotFoundError:
        download(mptopo_url, mptopo_file)
        mptopo_file_parse= ET.parse(mptopo_file)
    return(mptopo_file_parse)




def check_uniprot_sequence(mptopo_uniprot_id, sequence):
    print(f"Looking for {mptopo_uniprot_id} in the database")
    #MPTOPO is not restricted to humans so a get command will throw an error.
    # Rather than making an error exception, we will filter the results, check there is only one result and then act on that.
    db_seq=None
    database_record=Protein.objects.filter(uniprot_id=mptopo_uniprot_id)
    if database_record.count()==1:
        for i in database_record:
            db_seq=i.full_sequence
    elif database_record.count()>1:
        print("more than one database record match this query. Fix immediately")
    if str(db_seq)==str(sequence):
        return(True)
    elif database_record.count()!=1:
        return(False)
    else:
        print(f"Sequence mismatch in {mptopo_uniprot_id}:")
        print(f"Uniprot:{str(db_seq)}")
        print(f"MPTOPO:{sequence}")
        return(False)



def mptopo_check():
    '''
    Checks the MPTOPO xml file for transmemembrane regions mapped to a UniProt ID.
    '''
    mptopo=mptopo_file_check()
    # Sequences don't exactly match UniProt
    evidence_type = str("MPTOPO")

    # for item in root.findall("item"):
    #   ElementTree.dump(item)

    for node in mptopo.findall('.//mptopoProtein'):

        features = node.getchildren()
        for feature in features:
            if str(feature.tag) == str("sequence"):
                mptopo_sequence=feature.text

        for feature in features:
            if str(feature.tag) == str("uniprotNumber"):
                query_id=str(feature.text)
                # print("Matches query...")
                sequence_match_pass=check_uniprot_sequence(query_id, mptopo_sequence)
                if sequence_match_pass == True:
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



                            #Find the total tmhs before we parse more thoroughly
                            tmh_number_count=0
                            for tm_number, tmh_segment in enumerate(tmhs):
                                tmh_locations = tmh_segment.getchildren()
                                tmh_number_count=tmh_number_count+1
                                # get around 0 base counting
                            total_tmh_number=tmh_number_count



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

                                #tmh_list.append([query_id, tmh_start, tmh_stop, tmh_topology, evidence_type])
                                #as far as I can tell there is no standard membrane location in the xml
                                membrane_location="Unknown"
                                # lets get the slices from the sequence
                                n_ter_seq=mptopo_sequence[int(tmh_start)-5:int(tmh_start)]
                                tmh_sequence=mptopo_sequence[int(tmh_start):int(tmh_stop)]
                                c_ter_seq=mptopo_sequence[int(tmh_stop):int(tmh_stop)+5]
                                full_sequence=mptopo_sequence
                                # This is a dirty assumption that the alpha list from MPTOPO already has the filter
                                tm_type="Helix"
                                tmh_list.append([query_id, tmh_number, total_tmh_number, tmh_start , tmh_stop, tmh_topology, evidence_type, membrane_location, n_ter_seq, tmh_sequence, c_ter_seq, evidence_type, full_sequence, tm_type])
                                tmh_list=integrity_check(tmh_list)
                                tmh_list=clash_correction(tmh_list)
                                print(tmh_list)
                                tmh_to_database(tmh_list)
                                #return(tmh_list)
    return()




# This gubbins is needed to make the file run.
def run():
    mptopo_check()
