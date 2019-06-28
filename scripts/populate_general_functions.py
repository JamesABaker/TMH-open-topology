from requests import get
from requests.exceptions import ConnectionError
from datetime import date

"""
These functions are used repeatedly throughout the population process.
To keep things consistent and easy to manage, they have all been bundled into one file.
"""

time_threshold = 7
today = date.today()
todaysdate = today.strftime("%d_%m_%Y")

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

def input_query_get():
    '''
    Returns a list of uniprot ids.
    '''
    # In full scale mode it will take a long time which may not be suitable for development.
    input_query_list = get_uniprot()
    # Here we will just use a watered down list of tricky proteins. Uncomment this line for testing the whole list.
    #input_query_list = ["P01850", "P22760", "Q5K4L6","Q7Z5H4", "P32897", "Q9NR77", "P31644", "Q9NS61"]

    return(input_query_list)


def download(url, file_name):
    '''
    Downloads the content of a url to a local file.
    '''
    # open in binary mode
    with open(file_name, "wb") as file:
        # get request
        response = None
        while response is None:
            try:
                print("Donwloading", url, "to", file_name, "...")
                # connect
                response = get(url)
            except (ConnectionError, urllib.error.HTTPError) as e:
                print("Connection dropped during download.")

        # write to file
        file.write(response.content)


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


def input_query_process(input_query):
    '''
    This returns an input list (such as an opened list.txt file from uniprot) and returns a list and a set of the ids.
    '''

    input_queries = []
    for query_number, a_query in enumerate(input_query):
        a_query = clean_query(a_query)
        print("Checking cache/downloading", a_query, ",",
              query_number + 1, "of", len(input_query), "records...")

        input_queries.append(clean_query(a_query))

    input_query_set = set(input_queries)

    return([input_queries, input_query_set])
