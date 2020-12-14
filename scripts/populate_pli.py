# Shell Plus Model Imports
from django.contrib.admin.models import LogEntry
from django.contrib.auth.models import Group, Permission, User
from django.contrib.contenttypes.models import ContentType
from django.contrib.sessions.models import Session
from tmh_db.models import (
    Database_Metadata,
    Disease,
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

# Shell Plus Django Imports
from django.core.cache import cache
from django.conf import settings
from django.contrib.auth import get_user_model
from django.db import transaction
from django.db.models import (
    Avg,
    Case,
    Count,
    F,
    Max,
    Min,
    Prefetch,
    Q,
    Sum,
    When,
    Exists,
    OuterRef,
    Subquery,
)
from django.utils import timezone
from django.urls import reverse


def pli_parser(pli_file):
    with open(pli_file) as inputfile:
        for line_number, pli_entry in enumerate(inputfile):
            if line_number > 0:
                column_entries = pli_entry.split()
                uniprot = column_entries[3].split("-")[0]
                pli = column_entries[5]
                missense = column_entries[8]
                print(f"{uniprot}, {pli}, {missense}")
                if "NA" not in pli and "NA" not in missense:
                    database_update(
                        marcia_uniprot_id=uniprot,
                        pli_score=pli,
                        missense_score=missense,
                    )


def database_update(marcia_uniprot_id=None, pli_score=None, missense_score=None):
    Protein.objects.filter(uniprot_id=marcia_uniprot_id).update(
        oe_mis_upper_gn=missense_score, pLI_gn=pli_score
    )


def run():
    pli_parser(
        "scripts/external_datasets/ensembl_uniprot_MANE_metrics-ALL-07-10-2020_simple_uniprotSPC"
    )
