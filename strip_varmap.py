import sys

filename = sys.argv[1]

aas = [
    "K",
    "R",
    "E",
    "D",
    "Q",
    "H",
    "N",
    "P",
    "Y",
    "W",
    "C",
    "M",
    "T",
    "S",
    "G",
    "V",
    "F",
    "A",
    "I",
    "L",
]


def varmap_columns_and_keys(column_headers):
    column_headers = column_headers.split()
    varmap_col_dictionary = {}
    for column_number, column_title in enumerate(column_headers):
        varmap_col_dictionary[column_title] = column_number
    # print(varmap_col_dictionary)
    return varmap_col_dictionary


def varmap_header_dict(varmap_file):
    with open(varmap_file, encoding="ISO-8859-1") as inputfile:
        for line_number, line in enumerate(inputfile):
            # This is the header line
            if line_number == 0:
                varmap_headers = varmap_columns_and_keys(line)
    return varmap_headers


headers = varmap_header_dict(filename)
print(headers)
with open(filename, encoding="ISO-8859-1") as file:
    for line_number, line in enumerate(file):
        line = line.split()
        if line[headers["UNIPROT_AA"]] in aas:
            uni_acc=line[headers["UNIPROT_ACCESSION"]]
            seq_n=line[headers["SEQ_NO"]]
            seq_aa=line[headers["UNIPROT_AA"]]
            pdb_acc=line[headers["CLOSEST_PDB_CODE"]]
            pdb_n=line[headers["RES_NUM"]]
            pdb_aa=line[headers["RES_NAME"]]

            print(
                f'uniprot: {uni_acc} \n sequence {seq_n} \n pdb {pdb_acc} \n sequence {pdb_n}'
            )
