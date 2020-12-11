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


def structure_residue_uniprot_residue_check():
    a = (
        Structural_residue.objects.all()
        .distinct("pk")
        .values_list(
            "structure__pdb_id",
            "residue__protein__uniprot_id",
            "structure_aa",
            "residue__amino_acid_type",
        )
    )
    test_pass = True
    mismatches = []
    for i in a:
        if i[2] != i[3]:
            test_pass = False
            mismatches.append(i)

    return (test_pass, mismatches)


def run():
    residue_structure_mismatches = structure_residue_uniprot_residue_check()
    if residue_structure_mismatches[0] == True:
        print("All structures match their uniprot residue aa types.")
    elif residue_structure_mismatches[0] == False:
        print(
            len(residue_structure_mismatches[1]),
            "structural residues did not match between their structure and their sequence.",
        )
