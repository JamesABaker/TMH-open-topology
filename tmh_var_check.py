import sys

var_file = str(sys.argv[1])
tmh_file = str(sys.argv[2])

# Load the  varsite tsv file then the first table (tmh_reference_table).
var_results = []
with open(var_file) as inputfile:
    for line in inputfile:
        var_results.append(line.strip().split('\t'))

tmh_results = []
with open(tmh_file) as inputfile:
    for line in inputfile:
        tmh_results.append(line.strip().split(','))

for var_database_entry in var_results:
    # print(var_database_entry)
    # Don't read header line
    if var_database_entry == var_results[0]:
        pass
    else:
        try:
            # CHROMOSOME = str(var_database_entry[0])
            # COORDS = str(var_database_entry[1])
            USER_BASE = str(var_database_entry[2])
            USER_VARIANT = str(var_database_entry[3])
            # ENSEMBL_BASE = str(var_database_entry[4])
            # VEP_CODING_BASE = str(var_database_entry[5])
            # GENE = str(var_database_entry[6])
            # GENE_ACC = str(var_database_entry[7])
            # REFSEQ_GENE_ACC = str(var_database_entry[8])
            # TRANSCRIPT = str(var_database_entry[9])
            # REFSEQ_TRANSCRIPT = str(var_database_entry[10])
            # STRAND_DIR = str(var_database_entry[11])
            # CODON_CHANGE = str(var_database_entry[12])
            # VEP_AA = str(var_database_entry[13])
            # UNIPROT_AA = str(var_database_entry[14])
            # AA_CHANGE = str(var_database_entry[15])
            # POLYPHEN_SCORE = str(var_database_entry[16])
            # SIFTS_SCORE = str(var_database_entry[17])
            UNIPROT_ACCESSION = str(var_database_entry[18])
            # PROTEIN_NAME = str(var_database_entry[19])
            SEQ_NO = str(var_database_entry[20])
            # CHANGE_TYPE = str(var_database_entry[21])
            # ALL_TRANSCRIPTS = str(var_database_entry[22])
            # NOTE = str(var_database_entry[23])
            # GNOMAD_AF = str(var_database_entry[24])
            # NEGATIVE = str(var_database_entry[25])
            # SYNONYMOUS = str(var_database_entry[26])
            # HAVE_PDB = str(var_database_entry[27])
            # PDB_UNIPROT_MATCH = str(var_database_entry[28])
            # CLOSEST_PDB_CODE = str(var_database_entry[29])
            # PDB_CHAIN = str(var_database_entry[30])
            # PDB_PROTEIN_NAME = str(var_database_entry[31])
            # PDB_EXPT_TYPE = str(var_database_entry[32])
            # PDB_RESOLUTION = str(var_database_entry[33])
            # PDB_RFACT = str(var_database_entry[34])
            # PDB_UNIPROT_ACC = str(var_database_entry[35])
            # PDB_IDENTITY = str(var_database_entry[36])
            # PDB_SW_SCORE = str(var_database_entry[37])
            # PDB_E_VALUE = str(var_database_entry[38])
            # RES_NAME = str(var_database_entry[39])
            # RES_NUM = str(var_database_entry[40])
            # SST = str(var_database_entry[41])
            # CAT_RES = str(var_database_entry[42])
            # DISULPHIDE = str(var_database_entry[43])
            # NTO_DNA = str(var_database_entry[44])
            # NTO_LIGAND = str(var_database_entry[45])
            # NTO_METAL = str(var_database_entry[46])
            # NTO_PROTEINNPDB_RES = str(var_database_entry[47])
            # LIGANDS = str(var_database_entry[48])
            # METALS = str(var_database_entry[49])
            # PFAM_DOMAIN = str(var_database_entry[50])
            # PFAM_NAME = str(var_database_entry[51])
            # CATH_DOMAIN = str(var_database_entry[52])
            # CATH_NAME = str(var_database_entry[53])
            # RES_CONSERVATION = str(var_database_entry[54])
            # NCONS_SEQS = str(var_database_entry[55])
            # DISEASES = str(var_database_entry[56])
            # DISEASE_VARIANTS = str(var_database_entry[57])
            # NVARIANTS = str(var_database_entry[58])
            # NAT_VARIANTS = str(var_database_entry[59])
        except(IndexError):
            #print("Not enough datapoints in line.")
            pass
        for tmh_database_entry in tmh_results:
            # print(tmh_database_entry)
            # Don't read header line
            if tmh_database_entry == tmh_results[0]:
                pass
            else:
                # This list should get bigger as scores etc are added.
                query_id = str(tmh_database_entry[0])
                tmh_start = int(tmh_database_entry[1])
                tmh_stop = int(tmh_database_entry[2])
                tmh_topology = str(tmh_database_entry[3]).strip()
                evidence_type = str(tmh_database_entry[4])
                n_location = str(tmh_database_entry[5])
                n_ter_seq = str(tmh_database_entry[6])
                tmh_seq = str(tmh_database_entry[7])
                c_ter_seq = str(tmh_database_entry[8])

                #n_ter_seq_all =  str(n_ter_seq + n_ter_seq_all)
                #tmh_seq_all = str(tmh_seq + tmh_seq_all)
                #c_ter_seq_all = str(c_ter_seq + c_ter_seq_all)

                # Lets sort the flanks into inside outside
                if tmh_topology == str("Inside"):
                    in_seq = str(n_ter_seq)
                    out_seq = str(c_ter_seq)
                elif tmh_topology == str("Outside"):
                    in_seq = str(c_ter_seq)
                    out_seq = str(n_ter_seq)
                else:
                    #print("Topology missing.")
                    pass

                try:
                    # We want as much information to be passed onto the next table.

                    # Z coordinate is the absolute distance from the central TMH residue to the variant position. I wonder if direction information is useful at this point.
                    z_coord_n_c = int(int(tmh_start + (abs(tmh_stop-tmh_start)/2)) - int(SEQ_NO))
                    if tmh_topology == str("Inside"):
                        z_coord_in_out = z_coord_n_c
                    elif tmh_topology == str("Outside"):
                        z_coord_in_out = 0-z_coord_n_c

                    if UNIPROT_ACCESSION == query_id and int(SEQ_NO) >= tmh_start and int(SEQ_NO) <= tmh_stop:
                        print(query_id, ",", SEQ_NO, ",", USER_BASE, ",",
                              USER_VARIANT, ",", "TM", ",",z_coord_n_c, ",", z_coord_in_out, ",", evidence_type)
                    elif UNIPROT_ACCESSION == query_id and int(SEQ_NO) >= tmh_start - len(n_ter_seq) and int(SEQ_NO) <= tmh_start:
                        if tmh_topology == str("Inside"):
                            print(query_id, ",", SEQ_NO, ",", USER_BASE, ",", USER_VARIANT, ",",
                                  "N-terminal inside flank", ",",z_coord_n_c, ",", z_coord_in_out, ",", evidence_type)
                        elif tmh_topology == str("Outside"):
                            print(query_id, ",", SEQ_NO, ",", USER_BASE, ",", USER_VARIANT, ",",
                                  "N-terminal outside flank", ",",z_coord_n_c, ",", z_coord_in_out, ",", evidence_type)
                    elif UNIPROT_ACCESSION == query_id and int(SEQ_NO) >= tmh_stop and int(SEQ_NO) <= tmh_stop + len(c_ter_seq):
                        if tmh_topology == str("Inside"):
                            print(query_id, ",", SEQ_NO, ",", USER_BASE, ",", USER_VARIANT, ",",
                                  "C-terminal outside flank", ",",z_coord_n_c, ",", z_coord_in_out, ",", evidence_type)
                        elif tmh_topology == str("Outside"):
                            print(query_id, ",", SEQ_NO, ",", USER_BASE, ",", USER_VARIANT, ",",
                                  "C-terminal inside flank", ",",z_coord_n_c, ",", z_coord_in_out, ",", evidence_type)

                except(ValueError):
                    pass
