import sys

def varmap_header_dict(varmap_file):
    with open(varmap_file) as inputfile:
        for line_number, var_database_entry in enumerate(inputfile):
            # This is the header line
            if line_number == 0:
                varmap_headers = varmap_columns_and_keys(var_database_entry)

    return(varmap_headers)

varmap_headers
