from __future__ import division
import requests
import urllib
import shutil
import numpy as np
import os
import collections
import gzip
import json
import re
import defusedxml.ElementTree as ET
from defusedxml.lxml import fromstring
from urllib.error import URLError
import Bio
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
# env LDFLAGS="-I/usr/local/opt/openssl/include -L/usr/local/opt/openssl/lib" pip install psycopg2
from django.conf import settings
from django.db import models
from tmh_db.models import Database_Metadata, Uniref, Go, Structure, Structural_residue, Funfam_residue, Funfamstatus, Protein, Residue, Tmh, Tmh_deltag, Tmh_hydrophobicity, Tmh_residue, Tmh_tmsoc, Variant, Keyword, Binding_residue
from django.utils import timezone
from datetime import date
from scripts.populate_general_functions import *

print("Usage:\npython manage.py runscript populate --traceback")

# How many days should be allowed to not enforce updates
time_threshold = 7
today = date.today()
todaysdate = today.strftime("%d_%m_%Y")



def sifts_mapping(a_query):
    protein_object = Protein.objects.get(uniprot_id=a_query)

    # Download via the API
    print("Fetching sifts information")
    sifts_file = str(
        "scripts/external_datasets/sifts_mapping/" + a_query + ".json")
    sifts_url = f"https://www.ebi.ac.uk/pdbe/api/mappings/all_isoforms/{a_query}"
    download(sifts_url, sifts_file)
    # PARSE and add to database
    with open(sifts_file, 'r') as file:
        sifts_json = json.load(file)

    for record in sifts_json:
        print(record)
        pdb_codes = []
        for pdb_code in sifts_json[a_query]['PDB'].keys():
            pdb_codes.append(pdb_code)
        print(a_query, "maps to", pdb_codes)
        for pdb_code in pdb_codes:

            print("Processing and mapping into database", pdb_code)
            pdb_download_location = f"ftp://ftp.ebi.ac.uk/pub/databases/pdb/data/structures/divided/pdb/{pdb_code[1:3]}/pdb{pdb_code}.ent.gz"
            pdb_file_location = f"./scripts/external_datasets/pdb/{pdb_code}.pdb"

            try:
                pdb_str = gzip.open(urllib.request.urlopen(
                pdb_download_location)).read()

                record_for_database, created = Structure.objects.get_or_create(pdb_id=pdb_code)
                record_for_database.uniprot_protein_id.add(protein_object)

                structure_object = Structure.objects.get(uniprot_protein_id=protein_object, pdb_id=pdb_code)

                with open(pdb_file_location, 'w') as pdb_file:
                    pdb_file.write(pdb_str.decode("utf-8"))

                    structure_sequence_map = get_sequence_resid_chains_dict(pdb_code)  # [1]
                    # ('Q95460', 36): {'A': [15, 14, 'Asp'], 'C': [15, 14, 'Asp']}

                    residue_list = list(Residue.objects.filter(protein=protein_object).values())
                    for residue in residue_list:

                        try:
                            #print("Trying to find ", (a_query, residue_details["sequence_position"]),  "in", structure_sequence_map)
                            structural_residues_to_map = structure_sequence_map[(a_query, residue["sequence_position"])]

                            seq_residue = Residue.objects.get(protein=protein_object, sequence_position=residue["sequence_position"])
                            #print("Residues to map:", structural_residues_to_map)
                            for chain, positions in structural_residues_to_map.items():
                                pdb_chain = chain
                                pdb_position = positions[0]
                                author_position = positions[1]
                                structure_aa = Bio.SeqUtils.IUPACData.protein_letters_3to1[positions[2]]
                                #print("Mapping:,", pdb_chain, pdb_position, author_position)

                                record_for_database, created = Structural_residue.objects.update_or_create(
                                    structure=structure_object,
                                    # residue=seq_residue,
                                    pdb_position=pdb_position,
                                    pdb_chain=pdb_chain,
                                    author_position=author_position,
                                    structure_aa=structure_aa,
                                    uniprot_position=residue["sequence_position"],
                                    memprotmd_headgroups=None,
                                    memprotmd_tail=None
                                )

                                record_for_database.residue.add(seq_residue)


                        except KeyError:
                            pass
            except URLError:
                print("pdb code", pdb_code, "has no pdb file. Potentially there is an mmcif.")


def get_sequence_resid_chains_dict(pdb_code):
    """Returns a dict where keys are Uniprot res and values are pdb res chains"""

    url = ("ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/split_xml/"
           "{0}/{1}.xml.gz".format(pdb_code[1:3], pdb_code))
    xml_str = gzip.open(urllib.request.urlopen(url)).read()
    xml_str = re.sub(b'\sxmlns="[^"]+"', b'', xml_str, count=1)
    root = fromstring(xml_str)
    sequence_chain_dict = collections.defaultdict(dict)
    chain_resid_to_auth_dict = {}
    for entity in root.findall(".//entity"):
        for residue in entity.findall(".//residue"):
            uniprot_cross_rf = residue.find(
                "./crossRefDb[@dbSource='UniProt']")
            if uniprot_cross_rf is None:
                continue
            pdbe_resid = int(residue.attrib['dbResNum'])
            pdb_resid = residue.find(
                "./crossRefDb[@dbSource='PDB']").get('dbResNum')
            pdb_resid = int(re.sub('[^0-9]', '', pdb_resid)
                            ) if pdb_resid != 'null' else None
            pdb_resname = residue.find(
                "./crossRefDb[@dbSource='PDB']").get('dbResName').capitalize()
            u_resid = int(uniprot_cross_rf.attrib['dbResNum'])
            uniprot_id = uniprot_cross_rf.attrib['dbAccessionId']
            sequence_chain_dict[(uniprot_id, u_resid)][entity.attrib['entityId']] = [
                pdbe_resid, pdb_resid, pdb_resname]

    # print(chain_resid_to_auth_dict)
    return sequence_chain_dict


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

    inputs = input_query_process(input_query)
    input_queries = inputs[0]

    # Populate structures
    for a_query in input_query:
        a_query = clean_query(a_query)
        print("Mapping", a_query, "via SIFTS...")
        sifts_mapping(a_query)

# sftp://ebi-login.ebi.ac.uk/nfs/public2/release/thornton-www/software/html/PoreWalker PORE WALKER FILES
