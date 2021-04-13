import os
import Bio
from Bio.PDB import PDBParser
from Bio.PDB import PDBIO
from csv import DictReader
from scripts.populate_general_functions import clean_query
from tmh_db.models import (
    Disease,
    Database_Metadata,
    Flank,
    Flank_residue,
    Funfam,
    FunfamResidue,
    Go,
    Keyword,
    Non_tmh_helix,
    Non_tmh_helix_residue,
    Protein,
    Residue,
    Signal_peptide,
    Signal_residue,
    Structural_residue,
    Structure,
    SubcellularLocation,
    Tmh,
    Tmh_deltag,
    Tmh_hydrophobicity,
    Tmh_residue,
    Tmh_tmsoc,
    Uniref,
    Variant,
)


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


def prot_database_check(pdb_id_to_check=""):
    '''
    Checks that the protein is in the database.
    '''
    matching_pdbs = Structure.objects.filter(pdb_id=pdb_id_to_check)
    if len(matching_pdbs) > 0:
        return(True)
    else:
        return(False)


def parse(pdb_id="", pdb_filename=""):
    '''
    parses the opm pdb files into a list of membrane residues.
    '''
    with open(pdb_filename) as f:
        content = f.readlines()

    # Start at one so the first residue atom is not immediately
    # considered as a residue.

    clean_id = pdb_id_clean = os.path.splitext(pdb_id)[0]

    bilayer_cutoff = thickness(pdb_content=content)
    if bilayer_cutoff is False:
        return(False)

    if prot_database_check(clean_id) is True:
        membrane_error = find_error(clean_id)

        parser = PDBParser()  # recruts the parsing class
        structure = parser.get_structure(pdb_filename, pdb_filename)
        for model in structure:   # X-Ray generally only have 1 model, while more in NMR
            for chain in model:
                chain_id = chain.get_id()
                print(chain_id)
                for residue in chain:
                    residue_number = residue.get_id()[1]
                    
                    #Lets ignore the DUM membrane placeholder HETATM
                    if str(residue.get_resname()) != str("DUM"):
                        for atom in residue:
                            aa_type = residue.get_resname()

                            #print(atom.get_coord())
                            # this is the relative TMH pos
                            z = float(atom.get_coord()[2])
                            residue_atom_list.append(z)

                        position = membrane_check(
                                z_positions=residue_atom_list,
                                membrane_cutoff=bilayer_cutoff, error=membrane_error)
                        print(clean_id, chain_id, residue_number)
                        Structural_residue.objects.filter(
                            structure__pdb_id=clean_id, pdb_chain=chain_id, pdb_position=residue_number).update(opm_status=position)
                        residue_atom_list = []
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
