# Shell Plus Model Imports
from django.contrib.admin.models import LogEntry
from django.contrib.auth.models import Group, Permission, User
from django.contrib.contenttypes.models import ContentType
from django.contrib.sessions.models import Session
from tmh_db.models import Database_Metadata, Flank, Flank_residue, Funfam, Funfam_residue, Funfamstatus, Go, Keyword, Non_tmh_helix, Non_tmh_helix_residue, Protein, Residue, Signal_peptide, Signal_residue, Structural_residue, Structure, Subcellular_location, Tmh, Tmh_deltag, Tmh_hydrophobicity, Tmh_residue, Tmh_tmsoc, Uniref, Variant
# Shell Plus Django Imports
from django.core.cache import cache
from django.conf import settings
from django.contrib.auth import get_user_model
from django.db import transaction
from django.db.models import Avg, Case, Count, F, Max, Min, Prefetch, Q, Sum, When, Exists, OuterRef, Subquery
from django.utils import timezone
from django.urls import reverse
import pymol
from pymol import cmd

def residue_exists(chain, residue):
    if cmd.count_atoms("chain " + chain) > 0:
        return True
    elif cmd.count_atoms("chain " + residue) == 0:
        print("...chain not found.")
        return False

def f(query, color):
    q=query.values_list("pdb_chain", "pdb_position")
    for i in q:
        print(clean_query(str(i[1])), clean_query(str(i[0])), color)




def run():

    pymol.finish_launching(['pymol', '-q'])

    structure=Structure.objects.get(pdb_id="1rx0").pdb_id

    cmd.fetch(structure, type="pdb")

    with open("resi_list.txt") as f:
        lines = f.read().splitlines()
        residues=[]
        for i in lines:
            a=i.split()
            residues.append(a)

    #cmd.color("white", "resi " + i[0] + " and chain "+ i[1])

    coloured=0
    not_coloured=0
    for i in residues:
        print("Attempting to colour residue", i[0], "on chain", i[1], "...")

        if residue_exists(i[1], i[0]) == True:
            cmd.color(i[2], "chain " + i[1]+ " resi " + i[0])
            print("...coloured")
            coloured=coloured+1
        else:
            not_coloured=not_coloured+1

    print("Coloured:", coloured, "\nNot found:", not_coloured)



#a=Structural_residue.objects.filter(structure__pdb_id="2zw3", residue__variant__variant_source__contains="gnomAD").values_list("pdb_position","pdb_chain").distinct('pk')
