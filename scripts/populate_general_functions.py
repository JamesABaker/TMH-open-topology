import collections
from django.utils import timezone
import os.path
import numpy as np
from requests import get
from time import sleep
from tmh_db.models import Disease, Database_Metadata, Flank, Flank_residue, Funfam, FunfamResidue, Go, Keyword, Non_tmh_helix, Non_tmh_helix_residue, Protein, Residue, Signal_peptide, Signal_residue, Structural_residue, Structure, SubcellularLocation, Tmh, Tmh_deltag, Tmh_hydrophobicity, Tmh_residue, Tmh_tmsoc, Uniref, Variant
import time
# These functions are used repeatedly throughout the population process.
# To keep things consistent and easy to manage, they have all been bundled into one file.

time_threshold = 7
todaysdate =time.time() 

test_query_list = ["Q5VTH2", "Q8NHU3", "P01850", "P22760", "P18507", "Q5K4L6",
                   "Q7Z5H4", "O14925", "Q9NR77", "P31644", "Q9NS61", "P02748", "A0A075B6J2"]


def aa_baezo_order():
    '''
    Returns the amino acids in the order according to Baezo-Delgado paper.
    '''
    return(['K', 'R', 'E', 'D', 'Q', 'H', 'N', 'P', 'Y', 'W', 'C', 'M', 'T', 'S', 'G', 'V', 'F', 'A', 'I', 'L'])


def remove_duplicate_variants(list_of_variants):
    remove_duplicate_list_of_variants = set(list_of_variants)
    remove_duplicate_list_of_variants = list(remove_duplicate_list_of_variants)
    truncate_list=[]
    for variant in remove_duplicate_list_of_variants:
        truncate_list.append((variant[0], variant[1]))
    return(truncate_list)


def impossible_subs():
    aa_impossible_subs_dict = {
        "A": ['K', 'R', 'Q', 'H', 'N', 'Y', 'W', 'C', 'M', 'F', 'I', 'L'],
        "R": ['E', 'D', 'N', 'Y', 'V', 'F', 'A'],
        "N": ['R', 'E', 'Q', 'P', 'W', 'C', 'M', 'G', 'V', 'F', 'A', 'L'],
        "D": ['K', 'R', 'Q', 'P', 'W', 'C', 'M', 'T', 'S', 'F', 'I', 'L'],
        "C": ['K', 'E', 'D', 'Q', 'H', 'N', 'P', 'M', 'T', 'V', 'A', 'I', 'L'],
        "Q": ['D', 'N', 'Y', 'W', 'C', 'M', 'T', 'S', 'G', 'V', 'F', 'A', 'I'],
        "E": ['R', 'H', 'N', 'P', 'Y', 'W', 'C', 'M', 'T', 'S', 'F', 'I', 'L'],
        "G": ['K', 'Q', 'H', 'N', 'P', 'Y', 'M', 'T', 'F', 'I', 'L'],
        "H": ['K', 'E', 'W', 'C', 'M', 'T', 'S', 'G', 'V', 'F', 'A', 'I'],
        "I": ['E', 'D', 'Q', 'H', 'P', 'Y', 'W', 'C', 'G', 'A'],
        "L": ['K', 'E', 'D', 'N', 'Y', 'C', 'T', 'G', 'A'],
        "K": ['D', 'H', 'P', 'Y', 'W', 'C', 'S', 'G', 'V', 'F', 'A', 'L'],
        "M": ['E', 'D', 'Q', 'H', 'N', 'P', 'Y', 'W', 'C', 'S', 'G', 'F', 'A'],
        "F": ['K', 'R', 'E', 'D', 'Q', 'H', 'N', 'P', 'W', 'M', 'T', 'G', 'A'],
        "P": ['K', 'E', 'D', 'N', 'Y', 'W', 'C', 'M', 'G', 'V', 'F', 'I'],
        "S": ['K', 'E', 'D', 'Q', 'H', 'M', 'V'],
        "T": ['E', 'D', 'Q', 'H', 'Y', 'W', 'C', 'G', 'V', 'F', 'L'],
        "W": ['K', 'E', 'D', 'Q', 'H', 'N', 'P', 'Y', 'M', 'T', 'V', 'F', 'A', 'I'],
        "Y": ['K', 'R', 'E', 'Q', 'P', 'W', 'M', 'T', 'G', 'V', 'A', 'I', 'L'],
        "V": ['K', 'R', 'Q', 'H', 'N', 'P', 'Y', 'W', 'C', 'T', 'S'],
        "U": []
        }

    return(aa_impossible_subs_dict)


def open_uniprot(uniprot_list_file):
    with open(uniprot_list_file) as f:
        # Somehow this has already made a list.
        lines = f

        input_query = list(lines)
        # Entry is the first line, which will break later code as it is not a valid uniprot id!
        input_query = input_query[1:]
        return(input_query)


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


def get_uniprot():
    '''
    Downloads UniProt IDs from Human transmembrane proteins from UniProt.org.
    '''
    # Grab the input list
    print("Fetching UniProt TM protein IDs")
    # All human proteins
    uniprot_list_url = "https://www.uniprot.org/uniprot/?query=reviewed%3Ayes+AND+organism%3A%22Homo+sapiens+%28Human%29+%5B9606%5D%22&sort=score&columns=id,&format=tab"
    # All human transmembrane proteins
    #uniprot_list_url = "https://www.uniprot.org/uniprot/?query=reviewed%3Ayes+AND+organism%3A%22Homo+sapiens+%28Human%29+%5B9606%5D%22+AND+annotation%3A%28type%3Atransmem%29&sort=score&columns=id,&format=tab"
    # "https://www.uniprot.org/uniprot/?query=reviewed%3Ayes+AND+annotation%3A(type%3Atransmem)&sort=score&columns=id,&format=tab"
    # uniprot_list = 'https://www.uniprot.org/uniprot/?query=reviewed%3Ayes+AND+organism%3A"Homo+sapiens+(Human)+[9606]"+AND+annotation%3A(type%3Atransmem)&sort=score&columns=id,&format=tab'

    uniprot_list_file = "scripts/external_datasets/uniprot_bin/uniprot_list" + \
        todaysdate + ".txt"
    try:
        input_query = open_uniprot(uniprot_list_file)
    except:
        download(uniprot_list_url, uniprot_list_file)
        input_query = open_uniprot(uniprot_list_file)
        # This saves the request to a file for reasons beyond me.
        # So we now need to open the file to recover the items as a list.

    return(input_query)


def output(line_for_output):
    '''
    Adds a time stamped string to a log.txt file
    '''
    output_file = open("log.txt", "a")
    line_for_output = todaysdate+":"+str(line_for_output)
    #line_for_output = line_for_output.translate(None, "()',")
    line_for_output = line_for_output + "\n"
    output_file.write(line_for_output)
    output_file.close()


def input_query_get():
    '''
    Returns a list of uniprot ids.
    '''
    # In full scale mode it will take a long time which may not be suitable for development.
    input_query_list = get_uniprot()
    # Here we will just use a watered down list of tricky proteins. Uncomment this line for testing the whole list.
    # input_query_list = test_query_list
    # This protein currently throws an error in biopython parsing.
    blacklist = []
    with open('scripts/external_datasets/exclusion_list.txt') as f:
        blacklist_lines = f.read().splitlines()
        for i in blacklist_lines:
            blacklist.append(clean_query(i))

    input_set = []
    for i in input_query_list:
        input_set.append(clean_query(i))

    blacklist = set(blacklist)

    input_set = set(input_set)

    new_input_set = input_set - blacklist
    new_input_set = list(new_input_set)
    #print("Test", new_input_set)
    return(new_input_set)

def check_local_file(file):
    '''
    Checks if a local file exists.
    Returns True if the file exists, or False if it does not.
    '''
    exists = os.path.isfile(file)
    return(exists)

def download(url, file_name, pause=0):
    '''
    Downloads the content of a url to a local file.
    '''
    sleep(pause)
    # open in binary mode
    with open(file_name, "wb") as file:
        # get request
        response = None
        while response is None:
            try:
                #print("Donwloading", url, "to", file_name, "...")
                # connect
                response = get(url)
            except ConnectionError:
                print("Connection dropped during download.")

        # write to file
        file.write(response.content)


def clean_query(query):
    '''
    This aims to generate a clean ascii query of a viable UniProt ID from a
     dirty input like a user input.
    '''

    illegal_characters = ["!", "\n", " ", "@",
                          "'", ")", ",", "(", "[", "]", " "]
    for char in illegal_characters:
        query = query.replace(char, "")
    #This gets rid of the isoform dashes and defaults to canonical.
    a_clean_query = query.split('-')[0]
    
    print("Clean query result:", a_clean_query)
    return(a_clean_query)


def list_to_csv(a_list):
    '''
    This takes a list and prints as CSV
    '''
    a_list = str(a_list)
    illegal_characters = ["[", "]", "'", '"']
    for char in illegal_characters:
        a_list = a_list.replace(char, "")
    a_csv = a_list
    # print("Clean query result:", a_clean_query)
    return(a_csv)


def input_query_process(input_query):
    '''
    This returns an input list (such as an opened list.txt file from uniprot) and returns a list and a set of the ids.
    '''

    input_queries = []
    for query_number, a_query in enumerate(input_query):
        a_query = clean_query(a_query)
        #print("Checking cache/downloading", a_query, ",",query_number + 1, "of", len(input_query), "records...")

        input_queries.append(clean_query(a_query))

    input_query_set = set(input_queries)

    return([input_queries, input_query_set])


def heatmap_array(var_freq_dict, aa_order):
    var_freq = collections.Counter(var_freq_dict)
    large_array = []
    for aa_mut in aa_order:
        aa_array = []
        for aa_wt in aa_order:
            # This query is counter intuitive. The aa_wt is first in the tuple, the aa_mut is second. The aa_mut is first in the loop to make sure it is on the y axis.
            aa_array.append(float(var_freq[(aa_wt, aa_mut)]))
        large_array.append(aa_array)
    # print(np.array(large_array))
    return(np.array(large_array))


def uniref_to_uniprot(uniref_id):
    '''
    Returns the uniprot part of a uniref id.
    '''
    uniref_parts = uniref_id.split("_")
    if len(uniref_parts) == 2:
        uniprot_id = uniref_parts[1]
    return(uniprot_id)
