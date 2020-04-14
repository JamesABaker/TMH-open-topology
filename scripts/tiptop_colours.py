# Shell Plus Model Imports
from django.conf import settings
from django.contrib.admin.models import LogEntry
from django.contrib.auth import get_user_model
from django.contrib.auth.models import Group
from django.contrib.auth.models import Permission
from django.contrib.auth.models import User
from django.contrib.contenttypes.models import ContentType
from django.contrib.sessions.models import Session
from django.core.cache import cache
from django.db import transaction
from django.db.models import Avg
from django.db.models import Case
from django.db.models import Count
from django.db.models import Exists
from django.db.models import F
from django.db.models import Max
from django.db.models import Min
from django.db.models import OuterRef
from django.db.models import Prefetch
from django.db.models import Q
from django.db.models import Subquery
from django.db.models import Sum
from django.db.models import When
from django.urls import reverse
from django.utils import timezone

from scripts.populate_general_functions import *
from tmh_db.models import Database_Metadata
from tmh_db.models import Flank
from tmh_db.models import Flank_residue
from tmh_db.models import Funfam
from tmh_db.models import Funfam_residue
from tmh_db.models import Funfamstatus
from tmh_db.models import Go
from tmh_db.models import Keyword
from tmh_db.models import Non_tmh_helix
from tmh_db.models import Non_tmh_helix_residue
from tmh_db.models import Protein
from tmh_db.models import Residue
from tmh_db.models import Signal_peptide
from tmh_db.models import Signal_residue
from tmh_db.models import Structural_residue
from tmh_db.models import Structure
from tmh_db.models import Subcellular_location
from tmh_db.models import Tmh
from tmh_db.models import Tmh_deltag
from tmh_db.models import Tmh_hydrophobicity
from tmh_db.models import Tmh_residue
from tmh_db.models import Tmh_tmsoc
from tmh_db.models import Uniref
from tmh_db.models import Variant
# Shell Plus Django Imports

colors = {
    "disease": "ruby",
    "gnomAD": "forest",
    "featureless": "white",
    "tmh_spontaneous": "cyan",
    "tmh_non_spontaneous": "orange",
    "pore": "nitrogen"
}


def delta_g_tmh_color(structural_residue_object):
    score = structural_residue_object.residue.get().tmh_residue_set.filter(
        tmh_id__meta_tmh=True).first().tmh_id.tmh_deltag_set.get().test_score
    if score > 0:
        return(colors["tmh_spontaneous"])
    else:
        return(colors["tmh_non_spontaneous"])


def residue_exists(chain, residue):
    '''
    Checks to make sure the intended atoms are present.
    '''
    if cmd.count_atoms("chain " + chain) > 0:
        return True
    elif cmd.count_atoms("chain " + residue) == 0:
        print("...chain not found.")
        return False


def backbone_object_to_cmd(pdb_id="", chain="", residue_position="", color=colors["featureless"]):
    '''
    Outputs the pml style for a general cartoon feature.
    '''
    pml_output(str("color " + color + ", " + pdb_id + " and chain " +
                   chain + " and resi " + str(residue_position)), pdb_id)


def stick_object_to_cmd(pdb_id=None, chain=None, residue_position=None, color=colors["featureless"], transparency=0):
    '''
    Outputs the pml style for a variant.
    '''
    pml_output(str("show stick, " + pdb_id + " and chain " +
                   chain + " and resi " + str(residue_position)), pdb_id)
    pml_output(str("set stick_color, " + color + ", " + pdb_id +
                   " and chain " + chain + " and resi " + str(residue_position)), pdb_id)
    pml_output(str("set_bond stick_transparency, " + str(transparency) + ", " +
                   pdb_id + " and chain " + chain + " and resi " + str(residue_position)), pdb_id)


def query_to_stick_cmd(structural_residue_query, color):
    '''
    Turns the Django query into a list of cmds for a pml file.
    '''
    q = structural_residue_query.values(
        "structure__pdb_id", "pdb_chain", "pdb_position")
    for i in q:
        # return(clean_query(str(i[1])), clean_query(str(i[0])), color)
        if color == colors["gnomAD"]:
            res_transparency = 0.8
        else:
            res_transparency = 0
        stick_object_to_cmd(pdb_id=str(i["structure__pdb_id"]), chain=str(
            i["pdb_chain"]), residue_position=str(i["pdb_position"]), color=color, transparency=res_transparency)


def pml_output(line, pml_pdb):
    '''
    Puts a string as a line in the pml file.
    '''
    pml_file = str(pml_pdb + "_varTMH.pml")
    with open(pml_file, "a") as f:
        f.write(line)
        f.write("\n")


def color_structure(pdb):
    '''
    Performs a django query and sends to the results to the pml file by calling other functions.
    '''

    tmh = Structural_residue.objects.filter(
        structure__pdb_id=pdb, residue__tmh_residue__tmh_id__meta_tmh=True).distinct('pk')

    for i in tmh:
        tmh_color = delta_g_tmh_color(i)
        backbone_object_to_cmd(pdb_id=pdb, chain=i.pdb_chain,
                               residue_position=i.pdb_position, color=tmh_color)

    pore = Structural_residue.objects.filter(
        structure__pdb_id=pdb, pore_residue=True).distinct('pk')
    for i in pore:
        backbone_object_to_cmd(pdb_id=pdb, chain=i.pdb_chain,
                               residue_position=i.pdb_position, color=colors["pore"])

    gnomad = Structural_residue.objects.filter(
        structure__pdb_id=pdb, residue__variant__variant_source__contains="gnomAD").distinct('pk')
    query_to_stick_cmd(gnomad, colors["gnomAD"])

    disease = Structural_residue.objects.filter(
        structure__pdb_id=pdb, residue__variant__disease_status="d").distinct('pk')
    query_to_stick_cmd(disease, colors["disease"])

    print(f"{pdb},{tmh.count()},{pore.count()},{gnomad.count()},{disease.count()}")


def run():
    '''
    This nonsense is needed by django
    '''
    with open('pore_residue_tmp_structures.txt') as f:
        lines = f.read().splitlines()
    print("PDB, TMHs, Pores, gnomAD, Disease")
    for i in lines:
        structure_pdb_id = clean_query(i)
        structure = Structure.objects.get(pdb_id=structure_pdb_id)

        fetch_command = str("fetch " + structure_pdb_id + ', type=pdb')
        pml_output(fetch_command, structure_pdb_id)
        pml_output(str("hide all"), structure_pdb_id)
        pml_output(str("show cartoon, " + structure_pdb_id), structure_pdb_id)
        pml_output(
            str("color " + colors["featureless"] + ", " + structure_pdb_id), structure_pdb_id)
        #cmd.color("white", "resi " + i[0] + " and chain "+ i[1])

        color_structure(structure_pdb_id)


#a=Structural_residue.objects.filter(structure__pdb_id="2zw3", residue__variant__variant_source__contains="gnomAD").values_list("pdb_position","pdb_chain").distinct('pk')
