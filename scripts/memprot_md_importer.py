import os
import Bio
from Bio.PDB import PDBParser
from Bio.PDB import PDBIO
import csv
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
from requests.exceptions import ConnectionError


def thickness(pdb_content=""):
    """
    Returns the bilayer thickness as a float.
    """
    bilayer_thickness = []
    for line in pdb_content:
        if "REMARK" in line and "bilayer thickness" in line:
            bilayer_thickness.append(line.split()[-1])
    if len(bilayer_thickness) >= 1:
        return float(bilayer_thickness[0])

    else:
        print("Bilayer weirdness:", bilayer_thickness)
        return False


def membrane_check(z_positions=[], membrane_cutoff=None, error=0):
    """
    Checks if all the atoms, some of the atoms,
    or none of the atoms fall within the membrane cutoffs
    """
    membrane = []
    non_membrane = []
    interface = []
    for a_z in z_positions:
        if abs(a_z) <= membrane_cutoff - error:
            membrane.append(a_z)
        # This catches non membrane,
        # but within the error of the membrane boundary.
        elif abs(a_z) <= membrane_cutoff + error:
            interface.append(a_z)
        else:
            non_membrane.append(a_z)
    if len(membrane) > 0 and len(non_membrane) == 0:
        return "membrane"
    elif len(membrane) == 0 and len(non_membrane) > 0:
        return "globular"
    elif len(membrane) > 0 and len(non_membrane) > 0:
        return "interface"
    elif len(interface) > 0:
        return "interface"


def find_error(pdb_id_clean=None):  # open file in read mode
    """
    Scans a lookup CSV for error values of membrane thickness in OPM.
    """
    opm_csv_filename = "scripts/external_datasets/opm/proteins-2021-04-06.csv"
    with open(opm_csv_filename, "r") as read_obj:
        # pass the file object to DictReader() to get the DictReader object
        csv_dict_reader = DictReader(read_obj)
        # iterate over each line as a ordered dictionary
        for row in csv_dict_reader:
            if clean_query(str(row["pdbid"])) == clean_query(pdb_id_clean):
                error = float(row["thicknesserror"])
                print(error)
                return error
    print("Failed to find error, falling back on REMARK")
    return 0


def local_memprot_tail(pdb_id=""):
    memprotmd_url = f"http://memprotmd.bioch.ox.ac.uk/data/memprotmd/simulations/{pdb_id}_default_dppc/files/structures/group_tail.contacts.pdb"
    memprotmd_file = f"scripts/external_datasets/memprotmd/{pdb_id}_tail.pdb"
    download(memprotmd_url, memprotmd_file, pause=5)
    return str(memprotmd_file)


def local_mapping(pdb_id=""):
    memprotmd_url = f"http://memprotmd.bioch.ox.ac.uk/data/memprotmd/simulations/{pdb_id}_default_dppc/files/contacts/by_resid_postprocess.csv"
    memprotmd_file = f"scripts/external_datasets/memprotmd/{pdb_id}_mapping.csv"
    download(memprotmd_url, memprotmd_file, pause=5)
    return str(memprotmd_file)


def local_memprot_head(pdb_id=""):
    memprotmd_url = f"http://memprotmd.bioch.ox.ac.uk/data/memprotmd/simulations/{pdb_id}_default_dppc/files/structures/group_head.contacts.pdb"
    memprotmd_file = f"scripts/external_datasets/memprotmd/{pdb_id}_head.pdb"
    download(memprotmd_url, memprotmd_file, pause=5)
    return str(memprotmd_file)


def residue_mapping(pdb_code=None):
    clean_id = os.path.splitext(pdb_code)[0]
    mapping_csv = local_mapping(clean_id)
    map = {}
    print(clean_id)
    with open(mapping_csv) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=",")
        line_count = 0
        try:
            for n, row in enumerate(csv_reader):
                #Skip header line

                if n > 0:
                    try:
                        memprot_md_number=row[0]
                        chain_number = row[3]
                        residue_number_corrected = row[7]
                        map[str(memprot_md_number)] = (float(residue_number_corrected), chain_number)
                    except(ValueError):
                        pass
                elif str("No file is located at that path.") in str(row):
                    return(False)

                        # return(residue_number_corrected, chain_number)
        except(IndexError):
            #URL not found return False
            return(False)

    return(map)


def parse(pdb_id="", pdb_filename=""):
    """
    parses the opm pdb files into a list of membrane residues.
    """
    with open(pdb_filename) as f:
        content = f.readlines()

    # Start at one so the first residue atom is not immediately
    # considered as a residue.
    #print(pdb_id)
    clean_id = os.path.splitext(pdb_id)[0]

    memprotmd_to_pdb = residue_mapping(pdb_code=pdb_id)

    if memprotmd_to_pdb is False:
        print(f"Could not make sense of MemProtMD mapping file for {pdb_id}")
        return(False)
    else:
        tail_file = local_memprot_tail(pdb_id=clean_id)
        head_file = local_memprot_head(pdb_id=clean_id)

        parser = PDBParser()  # recruts the parsing class
        tail_structure = parser.get_structure(tail_file, tail_file)
        head_structure = parser.get_structure(head_file, head_file)
        structural_residues=Structural_residue.objects.filter(structure__pdb_id=clean_id)
        for (
            model
        ) in tail_structure:  # X-Ray generally only have 1 model, while more in NMR
            for chain in model:
                chain_id = chain.get_id()
                print(chain_id)
                for residue in chain:
                    residue_number = residue.get_id()[1]

                    bfactors = []
                    for atom in residue:
                        if atom.bfactor > 0:
                            bfactors.append(atom.bfactor)
                    try:
                        real_residue_number, chain_id = memprotmd_to_pdb[str(residue_number)]
                        if len(bfactors) > 0:
                            # This should not be here, it should be below the definition for residue_number and the parser should return a dictionary, but I am lazy and in a rush

                            structural_residues.filter(
                                structure__pdb_id=clean_id,
                                pdb_chain=chain_id,
                                pdb_position=real_residue_number,
                            ).update(memprotmd_tail=True)
                            print(
                                f"Residue {real_residue_number}, and memprotmd residue number {residue_number}"
                            )
                        else:
                            structural_residues.filter(
                                structure__pdb_id=clean_id,
                                pdb_chain=chain_id,
                                pdb_position=real_residue_number,
                            ).update(memprotmd_tail=False)
                    except(KeyError):
                        pass

        for (
            model
        ) in head_structure:  # X-Ray generally only have 1 model, while more in NMR
            for chain in model:
                chain_id = chain.get_id()
                for residue in chain:
                    residue_number = residue.get_id()[1]

                    print(
                        f"Residue {real_residue_number}, and memprotmd residue number {residue_number}"
                    )
                    bfactors = []
                    for atom in residue:
                        if atom.bfactor > 0:
                            bfactors.append(atom.bfactor)
                    try:
                        real_residue_number, chain_id = memprotmd_to_pdb[str(residue_number)]
                        if len(bfactors) > 0:

                            structural_residues.filter(
                                structure__pdb_id=clean_id,
                                pdb_chain=chain_id,
                                pdb_position=real_residue_number,
                            ).update(memprotmd_head=True)
                        else:
                            structural_residues.filter(
                                structure__pdb_id=clean_id,
                                pdb_chain=chain_id,
                                pdb_position=real_residue_number,
                            ).update(memprotmd_head=False)
                    except(KeyError):
                        pass
    return(True)


def run():
    tmh_structures=Structure.objects.filter(uniprot_protein_id__total_tmh_number__gte=1)
    directory = r"scripts/external_datasets/opm/"
    for filename in os.listdir(directory):
        if filename.endswith(".pdb"):
            # print("Opening...")
            # print(os.path.join(directory, filename))
            pdb_id_to_check = os.path.splitext(filename)[0]
            try:
                if len(tmh_structures.filter(pdb_id=pdb_id_to_check)) > 0:
                    parse(pdb_filename=os.path.join(directory, filename), pdb_id=filename)
            except ConnectionError as e:
                if len(tmh_structures.filter(pdb_id=pdb_id_to_check)) > 0:
                    parse(pdb_filename=os.path.join(directory, filename), pdb_id=filename)
        else:
            continue
