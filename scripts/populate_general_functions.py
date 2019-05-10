
def download(url, file_name):
    '''
    Downloads the content of a url to a local file.
    '''
    # open in binary mode
    with open(file_name, "wb") as file:
        # get request
        response = get(url)
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
