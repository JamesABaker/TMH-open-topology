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

# To run this file, run:
# python3 manage.py runscript populate --traceback
# Could just do a bash script? wget is more reliable than the python modules
print("Usage:\npython3 manage.py runscript populate --traceback")

# Print clean list


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
    benign = ["Unclassified", "Polymorphism", "Affects", "association", "Benign", "Benign/Likely_benign", "Conflicting_interpretations_of_pathogenicity",
              "drug_response", "no_interpretation_for_the_single_variant", "not_provided",  "other", "protective", "risk_factor", "Uncertain_significance"]
    if str(disease_type) in disease:
        pathogenicity = "d"
    elif str(disease_type) in benign:
        pathogenicity = "n"
    else:
        pathogenicity = "u"
        print("Unknown pathogenicity:", str(disease_type))
    return(pathogenicity)


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
## Currently not in use


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
#Currently not in use


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
                print("Membrane known: ", membrane_location)
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
    print("Checking", query_id, "in UniProt.")
    evidence_type = str("UniProt")
    tmh_list = []

    # The UniProt bin contains lots of uniprot files.
    # This bin should either routinely be flushed, or have some sort of timestamp.
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
        # features locations is a bit annoying as the start location needs +1 to match the sequence IO, but end is the correct sequence value.
        for i, f in enumerate(record.features):
            if f.type == feature_type:
                if "UnknownPosition" in str(f.location.start) or "UnknownPosition" in str(f.location.end):
                    pass
                    print(record.id, "Unknown position for TMH in record")
                else:
                    list_of_tmhs.append(int(f.location.start) + 1)
                    list_of_tmhs.append(int(f.location.end))
                    print("TMH in record")
        if len(list_of_tmhs) > 0: #Checks if it is a TM protein.
            # The protein for the database
            # This feels like it could go somewhere else...
            print(record.location)
            membrane_location = record.location


            record_for_database = Protein(uniprot_id=query_id, full_sequence=str(
                record.seq), membrane_type=membrane_location)
            record_for_database, created = Protein.objects.update_or_create(
                uniprot_id=query_id,
                defaults={
                    "full_sequence": str(record.seq),
                    "membrane_type": membrane_location
                }
            )

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

    # Now we can go through the record and write each TMH and any info to a file (or just print it!)
    for record in SeqIO.parse(filename, input_format):
        new_record = True
        tmd_count = 0
        for i, f in enumerate(record.features):
            if f.type == feature_type:
                if "UnknownPosition" in str(f.location.start) or "UnknownPosition" in str(f.location.end):
                    pass
                else:
                    tmd_count=tmd_count +1
                    tmh_number=tmd_count
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
                                    for location in inside_locations:
                                        if location in str(a_features.qualifiers):
                                            tmh_topology = "Inside"
                                            membrane_location = location
                                    for location in outside_locations:
                                        if location in str(a_features.qualifiers):
                                            tmh_topology = "Outside"
                                            membrane_location = location


                                    # The TMH for the database
                                    tmh_record_id = Protein.objects.get(uniprot_id=query_id)
                                    record_for_database = Tmh(tmh_id=str(
                                        query_id + "." + tmh_number), tmh_sequence=tmh_sequence, tmh_start=tmh_start, tmh_stop=tmh_stop, tmh_evidence='UniProt')
                                    record_for_database, created = Tmh.objects.update_or_create(
                                        uniprot_id=tmh_record_id,
                                        defaults={
                                            "tmh_id": str(query_id + "." + tmh_number),
                                            "tmh_sequence": tmh_sequence,
                                            "tmh_start": tmh_start,
                                            "tmh_stop": tmh_stop,
                                            "tmh_evidence":'UniProt'
                                        }
                                    )

                                    # record_for_database.save()

                                    # tmh_list.append([query_id, tmh_start, tmh_stop, tmh_topology,
                                    #                 evidence_type, membrane_location, n_ter_seq, tmh_sequence, c_ter_seq])
        return(tmh_list)
#uniprot_id = models.CharField(max_length=20, unique=True)
#full_sequence = models.TextField()
#membrane_type = models.CharField(max_length=100, default='')
##total_tmh_number = models.IntegerField(default=None)
#created_date = models.DateTimeField(default=timezone.now)


def clinvar_check(tmh_info):
    '''
    Checks if a tmh has any variants in the variant file and spews out a list of
    variants and their position in the tmh.
    '''
    var_source = "ClinVar"
    variants_in_tmh = []
    for var_database_entry in var_results:
        # print(var_database_entry)
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

        except(IndexError):
            #print("Not enough datapoints in line.")
            pass
            # This list should get bigger as scores etc are added.

        # This is a really messy way to get the structure of the tmh data.
        # But heck, I'd rather be human readable than have less code!
        # tmh_info=[query_id, tmh_start, tmh_stop, tmh_topology, , membrane_location, n_ter_seq, tmh_sequence, c_ter_seq, tmh_number, len(list_of_tmhs)]

        query_id = str(tmh_info[0])
        tmh_start = int(tmh_info[1])
        tmh_stop = int(tmh_info[2])
        tmh_topology = str(tmh_info[3]).strip()
        evidence_type = str(tmh_info[4])
        n_location = str(tmh_info[5])
        n_ter_seq = str(tmh_info[6])
        tmh_seq = str(tmh_info[7])
        c_ter_seq = str(tmh_info[8])
        tmh_number = int(tmh_info[9])
        tmd_total = int(tmh_info[10])

        var_record_location = str(SEQ_NO)
        var_record_id = USER_ID
        AA_CHANGE

        '''
        This could be worth investigating if isoforms are an issue
        # ISOFORMS!!!!
        # This bit is fiddly since there are isoforms. First, we need to establish if the record is the right line to save some time.
        id_match = False
        if query_id == UNIPROT_ACCESSION:

            var_record_id = UNIPROT_ACCESSION
            var_record_location = SEQ_NO
            id_match = True

        elif query_id in ALL_TRANSCRIPTS:
            id_match = True

            # ' / ' deliniates isoforms. ',' deliniates items in isoforms.
            # Example:
            #   ENST00000379389,-,P05161,21,S/N,*,ISG15,0,0.73,Missense variant / ENST00000458555,-,P05161,-,-,*,ISG15,0,0.73,Upstream gene variant

            list_of_transcripts = str(ALL_TRANSCRIPTS).split(" / ")
            for isoform in list_of_transcripts:
                this_isoform = isoform.split(",")
                if query_id == this_isoform[2]:
                    var_record_id = this_isoform[2]
                    var_record_location = this_isoform[3]
                else:
                    pass
        else:
            pass
        '''
        if query_id == UNIPROT_ACCESSION:
            variant_type = "Unknown"
            variant_review = "Unknown"

            # Is the variant disease causing?

            for i in clinvar_lines:
                #print("Is", int(i[-1]), "equal to", int(USER_ID), "?" )
                if int(i[-1]) == int(var_record_id):  #  (variant id is last column in summary)
                    # print("clinvar summary and snipclip finally found a hit for variant ",int(var_record_id))

                    variant_type = i[6]
                    variant_review = i[24]
                    #print(variant_type, variant_review)

                    # AlleleID
                    # Type
                    # Name
                    # GeneID
                    # GeneSymbol
                    # HGNC_ID
                    # ClinicalSignificance
                    # ClinSigSimple
                    # LastEvaluated
                    # RS# (dbSNP)
                    #nsv/esv (dbVar)
                    # RCVaccession
                    # PhenotypeIDS
                    # PhenotypeList
                    # Origin
                    # OriginSimple
                    # Assembly
                    # ChromosomeAccession
                    # Chromosome
                    # Start
                    # Stop
                    # ReferenceAllele
                    # AlternateAllele
                    # Cytogenetic
                    # ReviewStatus
                    # NumberSubmitters
                    # Guidelines
                    # TestedInGTR
                    # OtherIDs
                    # SubmitterCategories
                    # VariationID

                    # n_ter_seq_all =  str(n_ter_seq + n_ter_seq_all)
                    # tmh_seq_all = str(tmh_seq + tmh_seq_all)
                    # c_ter_seq_all = str(c_ter_seq + c_ter_seq_all)

                    # Lets sort the flanks into inside outside
                    if tmh_topology == str("Inside"):
                        in_seq = str(n_ter_seq)
                        out_seq = str(c_ter_seq)
                    elif tmh_topology == str("Outside"):
                        in_seq = str(c_ter_seq)
                        out_seq = str(n_ter_seq)
                    else:
                        # print("Topology missing.")
                        pass

                    if "-" not in str(AA_CHANGE):
                        aa_change = AA_CHANGE.split("/")
                        starting_residue = str(aa_change[0])
                        mutation_residue = str(aa_change[1])
                    else:
                        starting_residue = "-"
                        mutation_residue = "-"
                    # We want as much information to be passed onto the next table.

                    variant_class = str(disease_class(variant_type))
                    # Z coordinate is the absolute distance from the central TMH residue to the variant position. I wonder if direction information is useful at this point.
                    if var_record_location == "-":  # Lets do some exceotions without "try"
                        pass
                    else:
                        z_coord_n_c = int(
                            int(tmh_start + (abs(tmh_stop - tmh_start) / 2)) - int(var_record_location))
                        if tmh_topology == str("Inside"):
                            z_coord_in_out = z_coord_n_c
                        elif tmh_topology == str("Outside"):
                            z_coord_in_out = 0 - z_coord_n_c

                        #print (int(var_record_location), tmh_start, int(var_record_location), tmh_stop)
                        if int(var_record_location) > tmh_start and int(var_record_location) <= tmh_stop:
                            variant = [query_id, var_record_location, starting_residue, mutation_residue,
                                       "TM", z_coord_n_c, z_coord_in_out, tmh_number, tmd_total, var_source, variant_class, variant_type, variant_review]
                            variants_in_tmh.append(variant)
                            print_list(variant)
                        elif int(var_record_location) >= tmh_start - len(n_ter_seq) and int(var_record_location) <= tmh_start:
                            if tmh_topology == str("Inside"):
                                variant = [query_id, var_record_location, starting_residue, mutation_residue,
                                           "N-terminal inside flank", z_coord_n_c, z_coord_in_out, tmh_number, tmd_total, var_source, variant_class, variant_type, variant_review]
                                variants_in_tmh.append(variant)
                                print_list(variant)
                            elif tmh_topology == str("Outside"):
                                variant = [query_id, var_record_location, starting_residue, mutation_residue,
                                           "N-terminal outside flank", z_coord_n_c, z_coord_in_out, tmh_number, tmd_total, var_source, variant_class, variant_type, variant_review]
                                variants_in_tmh.append(variant)
                                print_list(variant)
                        elif int(var_record_location) >= tmh_stop and int(var_record_location) <= tmh_stop + len(c_ter_seq):
                            if tmh_topology == str("Inside"):
                                variant = [query_id, var_record_location, starting_residue, mutation_residue,
                                           "C-terminal outside flank", z_coord_n_c, z_coord_in_out, tmh_number, tmd_total, var_source, variant_class, variant_type, variant_review]
                                variants_in_tmh.append(variant)
                                print_list(variant)
                            elif tmh_topology == str("Outside"):
                                variant = [query_id, var_record_location, starting_residue, mutation_residue,
                                           "C-terminal inside flank", z_coord_n_c, z_coord_in_out, tmh_number, tmd_total, var_source, variant_class, variant_type, variant_review]
                                variants_in_tmh.append(variant)
                                print_list(variant)

    return(variants_in_tmh)


def humsavar_var_check(tmh_info):
    var_source = "Humsavar"
    variants_in_tmh = []
    query_id = str(tmh_info[0])
    tmh_start = int(tmh_info[1])
    tmh_stop = int(tmh_info[2])
    tmh_topology = str(tmh_info[3]).strip()
    evidence_type = str(tmh_info[4])
    n_location = str(tmh_info[5])
    n_ter_seq = str(tmh_info[6])
    tmh_seq = str(tmh_info[7])
    c_ter_seq = str(tmh_info[8])
    tmh_number = int(tmh_info[9])
    tmd_total = int(tmh_info[10])

    filename = str("uniprot_bin/" + query_id + ".txt")
    input_format = "swiss"
    feature_type = "TRANSMEM"
    subcellular_location = "TOPO_DOM"

    # print("Checking ", query_id, "in humsavar.txt.")
    for record in SeqIO.parse(filename, input_format):
        tmhs = 0
        for i, feature in enumerate(record.features):
            if feature.type == 'VARIANT':
                for entry in humsavar_list:
                    # print(entry)
                    if str(entry[2]) == str(feature.id):
                        #variant_types=[str('Disease'), str('Polymorphism'), str('Unclassified')]
                        variant_type = str(entry[4])
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
                                    original_residue = str(feature.qualifiers)[
                                        char_num - 2]
                                    variant_residue = str(feature.qualifiers)[
                                        char_num + 3]
                        record_for_database = Variant(aa_wt = original_residue, aa_mut = variant_residue, residue = query_id+".residue."+var_record_location)
                        record_for_database, created = Variant.objects.update_or_create(
                            uniprot_id=query_id,
                            defaults={
                                "aa_wt": str(query_id + "." + tmh_number),
                                "aa_mut": tmh_sequence,
                            }
                        )

                        if tmh_start <= feature.location.end and tmh_stop >= feature.location.start:
                            # print("Transmembrane variant!")
                            # print(original_residue,"to", variant_residue)
                            # print(record.id)
                            tm_variant = True

                            #variant_record(original_residue, variant_residue)

                            var_record_location = str(feature.location.end)



                            # We want as much information to be passed onto the next table.
                            variant_class = str(disease_class(variant_type))
                            # Z coordinate is the absolute distance from the central TMH residue to the variant position. I wonder if direction information is useful at this point.
                            z_coord_n_c = int(
                                int(tmh_start + (abs(tmh_stop - tmh_start) / 2)) - int(var_record_location))
                            if tmh_topology == str("Inside"):
                                z_coord_in_out = z_coord_n_c
                            elif tmh_topology == str("Outside"):
                                z_coord_in_out = 0 - z_coord_n_c

                            if int(var_record_location) > tmh_start and int(var_record_location) <= tmh_stop:
                                variant = [query_id, var_record_location, original_residue, variant_residue, "TM", z_coord_n_c,
                                           z_coord_in_out, tmh_number, tmd_total, var_source, variant_class, variant_type, variant_review]
                                variants_in_tmh.append(variant)
                                print_list(variant)
                            elif int(var_record_location) >= tmh_start - len(n_ter_seq) and int(var_record_location) <= tmh_start:
                                if tmh_topology == str("Inside"):
                                    variant = [query_id, var_record_location, original_residue, variant_residue, "N-terminal inside flank",
                                               z_coord_n_c, z_coord_in_out, tmh_number, tmd_total, var_source, variant_class, variant_type, variant_review]
                                    variants_in_tmh.append(variant)
                                    print_list(variant)
                                elif tmh_topology == str("Outside"):
                                    variant = [query_id, var_record_location, original_residue, variant_residue, "N-terminal outside flank",
                                               z_coord_n_c, z_coord_in_out, tmh_number, tmd_total, var_source, variant_class, variant_type, variant_review]
                                    variants_in_tmh.append(variant)
                                    print_list(variant)
                            elif int(var_record_location) >= tmh_stop and int(var_record_location) <= tmh_stop + len(c_ter_seq):
                                if tmh_topology == str("Inside"):
                                    variant = [query_id, var_record_location, original_residue, variant_residue, "C-terminal outside flank",
                                               z_coord_n_c, z_coord_in_out, tmh_number, tmd_total, var_source, variant_class, variant_type, variant_review]
                                    variants_in_tmh.append(variant)
                                    print_list(variant)
                                elif tmh_topology == str("Outside"):
                                    variant = [query_id, var_record_location, original_residue, variant_residue, "C-terminal inside flank",
                                               z_coord_n_c, z_coord_in_out, tmh_number, tmd_total, var_source, variant_class, variant_type, variant_review]
                                    variants_in_tmh.append(variant)
                                    print_list(variant)


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


def get_uniprot():
    # Grab the input list
    print("Fetching UniProt TM protein IDs")
    #uniprot_list = "https://www.uniprot.org/uniprot/?query=reviewed%3Ayes+AND+annotation%3A(type%3Atransmem)&sort=score&columns=id,&format=tab"
    uniprot_list = 'https://www.uniprot.org/uniprot/?query=reviewed%3Ayes+AND+organism%3A"Homo+sapiens+(Human)+[9606]"+AND+annotation%3A(type%3Atransmem)&sort=score&columns=id,&format=tab'

    uniprot_request = urllib.request.urlretrieve(str(uniprot_list))
    # This saves the request to a file for reasons beyond me.
    # So we now need to open the file to recover the items as a list.
    with open(uniprot_request[0]) as f:
        # Somehow this has already made a list.
        lines = f
        #lines = f.read().splitlines()

        input_query = list(lines)
        input_query = input_query[1:]
    return(input_query)


def run():
    ### Canonical script starts here ###

    # In full scale mode it will take forever. Here we will just use a watered down list of tricky proteins.
    # input_query=get_uniprot()
    input_query = ["P32897", "Q9NR77", "P31644", "P47869", "P28472", "P18507", "P05187"]

    # Parse the xml static files since this is the slowest part.
    # Ignore this for now -  we need to sort out uniprot before anything else!
    #topdb = ET.parse('topdb_all.xml')
    #mptopo = ET.parse('mptopoTblXml.xml')

    # Also, parse the variant files which can be massive.
    # humsavar table
    humsavar_list = []
    with open('humsavar.txt') as f:
        lines = f.read().splitlines()
        for i in lines:
            i = i.replace('  ', ' ')
            humsavar_list.append(i.split())
    # print(humsavar_list)

    # Load the  varsite tsv file from snip clip.
    var_results = []
    with open("dev_snipclip.tsv") as inputfile:
        for line_number, line in enumerate(inputfile):
            if line_number == 0:
                pass
            else:
                for a_query in input_query:
                    if clean_query(a_query) in str(line):
                        var_results.append(line.strip().split('\t'))
    # print(var_results)

    # Load the clinvar summary file
    clinvar_lines = []
    with open("variant_summary.txt") as inputfile:
        for line_number, line in enumerate(inputfile):
            if line_number == 0:
                pass
            else:
                for variant in var_results:
                    # print(str(variant[26])) # This should be simply the file ID. This filter kind of speeds it up if there are not many items.
                    if clean_query(variant[26]) in str(line):
                        clinvar_lines.append(line.strip().split('\t'))

    print(input_query)
    print("Starting TMH database population script...")
    os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'tmh_database.settings')
    for a_query in input_query:
        a_query = clean_query(a_query)
        # print(clean_query(a_query))
        ### OPM needs adding to here also. ###
        # print_list(mptopo_check(a_query))
        uniprot_check(a_query)
        # print_list(topdb_check(a_query))

    # populate()
