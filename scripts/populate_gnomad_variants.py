from scripts.populate_general_functions import *

# env LDFLAGS="-I/usr/local/opt/openssl/include -L/usr/local/opt/openssl/lib" pip install psycopg2
bulk_variants_to_add=[]

def varmap_columns_and_keys(column_headers):
    column_headers = column_headers.split('\t')
    varmap_col_dictionary = {}
    for column_number, column_title in enumerate(column_headers):
        varmap_col_dictionary[column_title] = column_number
    # #print(varmap_col_dictionary)
    return varmap_col_dictionary



def gnomad_process(varmap_file, input_query_set, version_name):

    varmap_results = []
    proteins=[]
    with open(varmap_file, encoding="ISO-8859-1") as inputfile:
        varmap_line = inputfile.readline()
        cnt = 0
        #for line_number, var_database_entry in enumerate(line):
        while varmap_line:

            # This is the header line
            cnt=cnt+1
            if cnt == 1:
                varmap_headers = varmap_columns_and_keys(varmap_line)
                uniprot_accession_index = varmap_headers["UNIPROT_ACCESSION"]
                varmap_user_id_index = varmap_headers["USER_ID"]
            else:
                #
                varmap_line = varmap_line.strip().split('\t')
                uniprot_accession =clean_query( str(varmap_line[uniprot_accession_index]))
                print("gnomAD in:", uniprot_accession)
                if str(uniprot_accession) in input_query_set:
                    #varmap_results.append(varmap_line)
                    ##print("Line {}: {}".format(cnt, varmap_line))

                    #now lets do the database stuff
                    variant_source = version_name
                    var_record_location = varmap_line[varmap_headers["SEQ_NO"]]
                    var_record_id = varmap_line[varmap_headers["USER_ID"]]
                    uniprot_record = clean_query(varmap_line[varmap_headers["UNIPROT_ACCESSION"]])
                    variant_review = "Unknown"
                    aa_wt = varmap_line[varmap_headers["UNIPROT_AA"]]
                    if len(varmap_line[varmap_headers["AA_CHANGE"]]) == 3 and "/" in varmap_line[varmap_headers["AA_CHANGE"]]:
                        aa_mut = varmap_line[varmap_headers["AA_CHANGE"]].split("/")[1]
                    else:
                        aa_mut = varmap_line[varmap_headers["AA_CHANGE"]]
                    disease_status = ""
                    disease_comments = ""
                    user_id = varmap_line[varmap_headers["USER_ID"]]
                    parsed_id=id_parse(user_id)
                    parsed_af=parsed_id["allele_frequency"]
                    if parsed_af == ".":
                        parsed_af=None
                    parsed_qc=parsed_id["qc"]
                    # Is the variant disease causing?
                    # This only applies to clinvar
                    print("Trying to store variant", var_record_id, "for", uniprot_accession, "to memory.")

                    var_to_database(uniprot_record, var_record_location, aa_wt, aa_mut,
                                    disease_status, disease_comments, variant_source, user_id, variant_maf=parsed_af, variant_qc=parsed_qc)
                    #print(cnt, ":")
                    proteins.append(uniprot_record)
            varmap_line = inputfile.readline()
    print(f'{len(proteins)} protein variants added, {len(set(proteins))} unique proteins')
    return()


def id_parse(varmap_id):
    varmap_id_dict={}
    varmap_id=varmap_id.split('zZz')
    print(varmap_id)
    varmap_order={
        0:"rsid",
        1:"qc",
        2:"allele_count",
        3:"allele_frequency"}
    for n, i in enumerate(varmap_id):
        if n in varmap_order:
            varmap_id_dict[varmap_order[n]]=i
    return(varmap_id_dict)


def var_to_database(uniprot_record, var_record_location, aa_wt, aa_mut, disease_status, disease_comments, variant_source, variant_source_id, variant_maf=None, variant_qc=None):
    '''
    Adds the variant from various external databases and flat files to the database in a standardised way.
    '''
    global bulk_variants_to_add
    if var_record_location == "-":
        #print("Unkown sequence location. Possibly intron: ", uniprot_record, var_record_location, aa_wt, "->", aa_mut, disease_status, disease_comments, variant_source, variant_source_id)
        pass
    elif aa_wt == "-":
        #print("Wildtype amino acid not defined. Assuming this is not an SNP: ", uniprot_record, var_record_location, aa_wt, "->", aa_mut, disease_status, disease_comments, variant_source, variant_source_id)
        pass
    elif len(aa_wt) > 1 or len(aa_mut) > 1:
        #print("More than a single residue changed. Assuming this is not an SNP: ", uniprot_record, var_record_location, aa_wt, "->", aa_mut, disease_status, disease_comments, variant_source, variant_source_id)
        pass
    elif "*" in str(aa_mut):
        #print("Stop codon introduced. This will change more than one residue:", uniprot_record,var_record_location, aa_wt, "->", aa_mut, disease_status, disease_comments, variant_source, variant_source_id)
        pass
    else:
        try:
            protein = Protein.objects.get(uniprot_id=uniprot_record)

            if int(var_record_location) <= len(str(protein.full_sequence)):

                residue_variant = Residue.objects.get(protein=protein, sequence_position=var_record_location)
                if str(residue_variant.amino_acid_type) == str(aa_wt):
                    #print("Adding ", uniprot_record, var_record_location, aa_wt, "->", aa_mut, disease_status, disease_comments, variant_source, "to database variant table.")
                    bulk_variants_to_add.append(
                        Variant(
                        residue=residue_variant,
                        aa_wt=aa_wt,
                        aa_mut=aa_mut,

                        disease_status=disease_status,
                        disease_comments=disease_comments,
                        variant_source=variant_source,
                        variant_source_id=variant_source_id,
                        maf=variant_maf,
                        qc=variant_qc
                    ))
                else:
                    #print("Mismatch between wild-type amino acids. UniProt:", str(residue_variant.amino_acid_type), str(variant_source), ":", str(aa_wt), "for record", uniprot_record, var_record_location, aa_wt, "->", aa_mut, disease_status, disease_comments, variant_source)
                    pass

            else:
                print("Variant position exceeds the length of the protein. Protein length:", len(str(protein.full_sequence)), "Variant position:", var_record_location, "for record", uniprot_record, var_record_location, aa_wt, "->", aa_mut, disease_status, disease_comments, variant_source)
                pass

        except(ObjectDoesNotExist):
            pass
    if len(bulk_variants_to_add)>9999:
        Variant.objects.bulk_create(bulk_variants_to_add)
        bulk_variants_to_add=[]

def run():
    input_query = input_query_get()
    # Also, parse the variant files which can be massive.
    #print("Starting TMH database population script...")
    os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'tmh_database.settings')

    ### Download uniprot files ###
    inputs = input_query_process(input_query)
    #input_queries = inputs[0]
    input_query_set = inputs[1]

    ### VarMap ###

    gnomad3_variant_varmap_file= "scripts/external_datasets/gnomad3_ALL.tsv"

    gnomad_process(gnomad3_variant_varmap_file, input_query_set, "gnomAD3")

    #gnomad2_variant_varmap_file= "scripts/external_datasets/gnomad_coding_regions2.tsv"

    #gnomad_process(gnomad2_variant_varmap_file, input_query_set, "gnomAD2")
