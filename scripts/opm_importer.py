import os


def thickness(pdb_content=""):
    '''
    Returns the bilayer thickness as a float.
    '''
    bilayer_thickness = []
    for line in pdb_content:
        if "REMARK" in line and "bilayer thickness" in line:
            bilayer_thickness.append(line.split()[-1])
    if len(bilayer_thickness) == 1:
        return(float(bilayer_thickness[0]))
    else:
        print("Bilayer weirdness:", bilayer_thickness)
        return(False)


def membrane_check(z_positions=[], membrane_cutoff=None):
    '''
    Checks if all the atoms, some of the atoms, or none of the atoms fall within the membrane cutoffs
    '''
    membrane = []
    non_membrane = []
    for a_z in z_positions:
        if abs(a_z) <= membrane_cutoff:
            membrane.append(a_z)
        else:
            non_membrane.append(a_z)
    if len(membrane) > 0 and len(non_membrane) == 0:
        return("membrane")
    elif len(membrane) == 0 and len(non_membrane) > 0:
        return("globular")
    elif len(membrane) > 0 and len(non_membrane) > 0:
        return("interface")


def parse(pdb_id="", pdb_filename=""):
    '''
    parses the opm pdb files into a list of membrane residues.
    '''
    with open(pdb_filename) as f:
        content = f.readlines()
    residue_atom_list = []
    # Start at one so the first residue atom is not immediately
    # considered as a residue.
    # This is a very not elegant way of doing this.
    previous_residue_number=1
    for line in content:  # each line is an atom
        line_content = line.split()
        #print(line_content)

        if line_content != []:  # annoying empty line exceptions.
            bilayer_cutoff = thickness(pdb_content=content)
            if bilayer_cutoff == False:
                return(False)
            if line_content[0] == 'ATOM':
                atom_number = line_content[1]
                atom_type = line_content[2]
                aa_type = line_content[3]
                chain = line_content[4]
                residue_number = line_content[5]
                #x = line_content[6]
                #y = line_content[7]
                try:
                    z = float(line_content[8])  # this is the TMH description
                    residue_atom_list.append(z)

                # A bug in which PDB format cause columns to bleed into one another.
                except(ValueError):
                    print("Borked PDB columns. Could not identify z axis value.")
                    return(False)
                if residue_number != previous_residue_number:
                    position=membrane_check(z_positions=residue_atom_list, membrane_cutoff=bilayer_cutoff)
                    print(pdb_id, chain, residue_number, position)
                    residue_atom_list = []
                previous_residue_number = residue_number

    return()


directory = r'scripts/external_datasets/opm/'
for filename in os.listdir(directory):
    if filename.endswith(".pdb"):
        # print("Opening...")
        # print(os.path.join(directory, filename))
        parse(pdb_filename=os.path.join(directory, filename), pdb_id=filename)

    else:
        continue
