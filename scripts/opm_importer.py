import os
from csv import DictReader
from scripts.populate_general_functions import clean_query


def thickness(pdb_content=""):
    '''
    Returns the bilayer thickness as a float.
    '''
    bilayer_thickness = []
    for line in pdb_content:
        if "REMARK" in line and "bilayer thickness" in line:
            bilayer_thickness.append(line.split()[-1])
    if len(bilayer_thickness) >= 1:
        return(float(bilayer_thickness[0]))

    else:
        print("Bilayer weirdness:", bilayer_thickness)
        return(False)


def membrane_check(z_positions=[], membrane_cutoff=None, error=0):
    '''
    Checks if all the atoms, some of the atoms,
    or none of the atoms fall within the membrane cutoffs
    '''
    membrane = []
    non_membrane = []
    interface = []
    for a_z in z_positions:
        if abs(a_z) <= membrane_cutoff-error:
            membrane.append(a_z)
        # This catches non membrane,
        # but within the error of the membrane boundary.
        elif abs(a_z) <= membrane_cutoff+error:
            interface.append(a_z)
        else:
            non_membrane.append(a_z)
    if len(membrane) > 0 and len(non_membrane) == 0:
        return("membrane")
    elif len(membrane) == 0 and len(non_membrane) > 0:
        return("globular")
    elif len(membrane) > 0 and len(non_membrane) > 0:
        return("interface")
    elif len(interface) > 0:
        return("interface")


def find_error(pdb_id_clean=None):    # open file in read mode
    '''
    Scans a lookup CSV for error values of membrane thickness in OPM.
    '''
    opm_csv_filename = "scripts/external_datasets/opm/proteins-2021-04-06.csv"
    with open(opm_csv_filename, 'r') as read_obj:
        # pass the file object to DictReader() to get the DictReader object
        csv_dict_reader = DictReader(read_obj)
        # iterate over each line as a ordered dictionary
        for row in csv_dict_reader:
            if clean_query(str(row["pdbid"])) == clean_query(pdb_id_clean):
                error = float(row["thicknesserror"])
                print(error)
                return(error)
    print("Failed to find error, falling back on REMARK")
    return(0)


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
    previous_residue_number = 1
    membrane_error = find_error(pdb_id_clean=os.path.splitext(pdb_id)[0])
    for line in content:  # each line is an atom
        line_content = line.split()
        if line_content != []:  # annoying empty line exceptions.
            bilayer_cutoff = thickness(pdb_content=content)
            if bilayer_cutoff is False:
                return(False)
            if line_content[0] == 'ATOM':
                atom_number = line_content[1]
                atom_type = line_content[2]
                aa_type = line_content[3]
                chain = line_content[4]
                residue_number = line_content[5]
                # x = line_content[6]
                # y = line_content[7]
                try:
                    z = float(line_content[8])  # this is the relative TMH pos
                    residue_atom_list.append(z)

                # A bug in which PDB format causes
                # columns to bleed into one another.
                except(ValueError):
                    print("Borked PDB columns; no clear z axis value.")
                    return(False)
                if residue_number != previous_residue_number:
                    position = membrane_check(
                        z_positions=residue_atom_list,
                        membrane_cutoff=bilayer_cutoff, error=membrane_error)
                    print(pdb_id, chain, residue_number, position)
                    residue_atom_list = []
                previous_residue_number = residue_number
    return()


def run():
    directory = r'scripts/external_datasets/opm/'
    for filename in os.listdir(directory):
        if filename.endswith(".pdb"):
            # print("Opening...")
            # print(os.path.join(directory, filename))
            parse(pdb_filename=os.path.join(
                directory, filename), pdb_id=filename)

        else:
            continue
