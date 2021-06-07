import os
import Bio
from Bio.PDB import PDBParser
from Bio.PDB import PDBIO
from csv import DictReader
from scripts.populate_general_functions import clean_query, download
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

list_of_vartmh_proteins = Structure.objects.all().values_list('pdb_id',
                                                              flat=True).distinct()


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


def local_memprot_tail(pdb_id=""):
    memprotmd_url = f'http://memprotmd.bioch.ox.ac.uk/data/memprotmd/simulations/{pdb_id}_default_dppc/files/structures/group_tail.contacts.pdb'
    memprot_file = f'scripts/external_datasets/memprotmd/{pdb_id}_tail.pdb'
    download(memprotmd_url, memprot_file, pause=5)
    return(str(memprot_file))


def local_memprot_head(pdb_id=""):
    memprotmd_url = f'http://memprotmd.bioch.ox.ac.uk/data/memprotmd/simulations/{pdb_id}_default_dppc/files/structures/group_head.contacts.pdb'
    memprot_file = f'scripts/external_datasets/memprotmd/{pdb_id}_head.pdb'
    download(memprotmd_url, memprot_file, pause=5)
    return(str(memprot_file))


def parse(pdb_id="", pdb_filename=""):
    '''
    parses the opm pdb files into a list of membrane residues.
    '''
    with open(pdb_filename) as f:
        content = f.readlines()

    # Start at one so the first residue atom is not immediately
    # considered as a residue.

    clean_id = os.path.splitext(pdb_id)[0]
    if clean_id in list_of_vartmh_proteins:

        tail_file = local_memprot_tail(pdb_id=clean_id)
        head_file = local_memprot_head(pdb_id=clean_id)

        parser = PDBParser()  # recruts the parsing class
        tail_structure = parser.get_structure(tail_file, tail_file)
        head_structure = parser.get_structure(head_file, head_file)

        for model in tail_structure:   # X-Ray generally only have 1 model, while more in NMR
            for chain in model:
                chain_id = chain.get_id()
                print(chain_id)
                for residue in chain:
                    residue_number = residue.get_id()[1]
                    bfactors = []
                    for atom in residue:
                        if atom.bfactor > 0:
                            bfactors.append(atom.bfactor)
                    if len(bfactors) > 0:
                        Structural_residue.objects.filter(
                            structure__pdb_id=clean_id, pdb_chain=chain_id, pdb_position=residue_number).update(memprotmd_tail=True)
                    else:
                        Structural_residue.objects.filter(
                            structure__pdb_id=clean_id, pdb_chain=chain_id, pdb_position=residue_number).update(memprotmd_tail=False)
        for model in head_structure:   # X-Ray generally only have 1 model, while more in NMR
            for chain in model:
                chain_id = chain.get_id()
                print(chain_id)
                for residue in chain:
                    residue_number = residue.get_id()[1]
                    bfactors = []
                    for atom in residue:
                        if atom.bfactor > 0:
                            bfactors.append(atom.bfactor)
                    if len(bfactors) > 0:
                        Structural_residue.objects.filter(
                            structure__pdb_id=clean_id, pdb_chain=chain_id, pdb_position=residue_number).update(memprotmd_head=True)
                    else:
                        Structural_residue.objects.filter(
                            structure__pdb_id=clean_id, pdb_chain=chain_id, pdb_position=residue_number).update(memprotmd_head=False)
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
