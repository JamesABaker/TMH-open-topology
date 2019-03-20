
fname='ian files list.txt'
with open(fname) as f:
    content = f.readlines()
# you may also want to remove whitespace characters like `\n` at the end of each line
content = [x.strip() for x in content]
funfam_files=content

tmh_clinvar='clinvar_test_tmh.csv'

results=[]
with open(tmh_clinvar) as inputfile:
    for line in inputfile:
        results.append(line.strip().split(','))


for entry in results:
    # print(entry)
    # Don't read header line
    if entry == results[0]:
        pass
    else:
        # query_id, ",", SEQ_NO, ",", USER_BASE, ",", USER_VARIANT, ",", "N-terminal outside flank", ",",z_coord_n_c, ",", z_coord_in_out, ",", evidence_type
        # This list should get bigger as scores etc are added.
        query_id = str(entry[0])
        seq_no = int(entry[1])
        USER_BASE = str(entry[2])
        USER_VARIANT = str(entry[3]).strip()
        N_terminal_outside_flank = str(entry[4])
        z_coord_n_c = str(entry[5])
        z_coord_in_out = str(entry[6])
        evidence_type = str(entry[7])
        for i in funfam_files:
            funfam_file=str("ian_funfam_alignments/sc-alignments.2018_09_28/"+i)
            with open(funfam_file) as inputfile:
                for line in inputfile:
                    if str(query_id) in str(line):
                        print(query_id, "in", funfam_file)
