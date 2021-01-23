# Shell Plus Model Imports
from django.conf import settings
from django.contrib.admin.models import LogEntry
from django.contrib.auth import get_user_model
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

from scripts.graphs import *
from scripts.populate_general_functions import *
# Shell Plus Django Imports

residues = [
    "A",
    "R",
    "N",
    "D",
    "C",
    "Q",
    "E",
    "G",
    "H",
    "I",
    "L",
    "K",
    "M",
    "F",
    "P",
    "S",
    "T",
    "W",
    "Y",
    "V",
]


def basic_num():
    print(
        "STATS OVERVIEW\n",
    )
    # Tmh.objects.exclude(residue__tmh_residue=None).filter(disease_status='d').distinct('pk').count()
    protein_num = Protein.objects.distinct("pk").count()
    print("UniProt IDs,", protein_num)
    residue_num = Residue.objects.distinct("pk").count()
    print("Residues,", residue_num)

    for i in residues:
        count = Residue.objects.filter(amino_acid_type=i).distinct("pk").count()
        print(i, ", ", count)

    tmh_residue_num = (
        Residue.objects.filter(tmh_residue__feature_location="TMH")
        .distinct("pk")
        .count()
    )
    print("TMH residues,", tmh_residue_num)

    for i in residues:
        count = (
            Residue.objects.filter(
                amino_acid_type=i, tmh_residue__feature_location="TMH"
            )
            .distinct("pk")
            .count()
        )
        print(i, ", ", count)

    # flanks
    flank_residue_num = (
        Residue.objects.filter(flank_residue__feature_location="Inside flank")
        .distinct("pk")
        .count()
    )
    print("inside flank residues,", flank_residue_num)

    for i in residues:
        count = (
            Residue.objects.filter(
                amino_acid_type=i, flank_residue__feature_location="Inside flank"
            )
            .distinct("pk")
            .count()
        )
        print(i, ", ", count)

    flank_residue_num = (
        Residue.objects.filter(flank_residue__feature_location="Outside flank")
        .distinct("pk")
        .count()
    )
    print("outside flank residues,", flank_residue_num)

    for i in residues:
        count = (
            Residue.objects.filter(
                amino_acid_type=i, flank_residue__feature_location="Outside flank"
            )
            .distinct("pk")
            .count()
        )
        print(i, ", ", count)

    flank_residue_num = (
        Residue.objects.filter(
            Q(flank_residue__feature_location="Inside flank")
            | Q(flank_residue__feature_location="Outside flank")
        )
        .distinct("pk")
        .count()
    )
    print("flank residues,", flank_residue_num)

    print("Helix residues")
    for i in residues:
        count = (
            Residue.objects.filter(amino_acid_type=i, non_tmh_helix_residue=True)
            .distinct("pk")
            .count()
        )
        print(i, ", ", count)

    non_tmh_residue_num = (
        Residue.objects.filter(tmh_residue=None).distinct("pk").count()
    )
    print("Non-TMH residues,", non_tmh_residue_num)
    for i in residues:
        count = (
            Residue.objects.filter(
                amino_acid_type=i, flank_residue=None, tmh_residue=None
            )
            .distinct("pk")
            .count()
        )
        print(i, ", ", count)

    structural_residue_num = (
        Residue.objects.exclude(structural_residue=None).distinct("pk").count()
    )
    print(
        "Residues with a map to at least one residue in a structure,",
        structural_residue_num,
    )
    non_structural_residue_num = (
        Residue.objects.filter(structural_residue=None).distinct("pk").count()
    )
    print(
        "Residues with no map to a residue in a structure,", non_structural_residue_num
    )
    tmh_structural_residue_num = (
        Residue.objects.exclude(tmh_residue=None, structural_residue=None)
        .distinct("pk")
        .count()
    )
    print(
        "TMH residues with a map to at least one residue in a structure,",
        tmh_structural_residue_num,
    )
    tmh_non_structural_residue_num = (
        Residue.objects.exclude(tmh_residue=None)
        .exclude(structural_residue=None)
        .distinct("pk")
        .count()
    )
    print(
        "TMH residues with no map to a residue in a structure,",
        tmh_non_structural_residue_num,
    )

    print("\n\nTMH boundaries\n")
    tmh_boundary_num = Tmh.objects.distinct("pk").count()
    print("TMP TMH boundaries,", tmh_boundary_num)
    uniprot_tmh_boundary = (
        Tmh.objects.filter(tmh_evidence="UniProt").distinct("pk").count()
    )
    print("UniProt TMH boundaries,", uniprot_tmh_boundary)
    topdb_tmh_boundary = Tmh.objects.filter(tmh_evidence="TOPDB").distinct("pk").count()
    print("TopDB TMH boundaries,", topdb_tmh_boundary)
    mptopo_tmh_boundary = (
        Tmh.objects.filter(tmh_evidence="MPTOPO").distinct("pk").count()
    )
    print("MPTOPO TMH boundaries,", mptopo_tmh_boundary)

    print("\n\nVariants\n")
    # Complex query example. How many variants are in the TMH?
    g_variants_num = (
        Variant.objects.filter(variant_source="gnomAD").distinct("pk").count()
    )
    print("gnomAD variants,", g_variants_num)
    tmh_g_variants_num = (
        Variant.objects.filter(
            residue__tmh_residue__feature_location="TMH", variant_source="gnomAD"
        )
        .distinct("pk")
        .count()
    )
    print("TMH gnomAD variants,", tmh_g_variants_num)

    inside_flank_g_variants_num = (
        Variant.objects.filter(residue__flank_residue__feature_location="Inside flank")
        .filter(variant_source="gnomAD")
        .distinct("pk")
        .count()
    )
    print("Inside flank gnomAD variants,", inside_flank_g_variants_num)

    outside_flank_g_variants_num = (
        Variant.objects.filter(residue__flank_residue__feature_location="Outside flank")
        .filter(variant_source="gnomAD")
        .distinct("pk")
        .count()
    )
    print("Outside flank gnomAD variants,", outside_flank_g_variants_num)

    non_tmh_g_variants_num = (
        Variant.objects.filter(
            residue__tmh_residue=None,
            residue__flank_residue=None,
            variant_source="gnomAD",
        )
        .distinct("pk")
        .count()
    )
    print("Non-TMH gnomAD variants,", non_tmh_g_variants_num)

    d_variants_num = Variant.objects.filter(disease_status="d").distinct("pk").count()
    print("Disease variants,", d_variants_num)
    tmh_d_variants_num = (
        Variant.objects.filter(
            residue__tmh_residue__feature_location="TMH", disease_status="d"
        )
        .distinct("pk")
        .count()
    )
    print("TMH disease variants,", tmh_d_variants_num)

    flank_d_variants_num = (
        Variant.objects.exclude(
            residue__tmh_residue=None, residue__tmh_residue__feature_location="TMH"
        )
        .filter(disease_status="d")
        .distinct("pk")
        .count()
    )

    inside_flank_d_variants_num = (
        Variant.objects.filter(residue__flank_residue__feature_location="Inside flank")
        .filter(disease_status="d")
        .distinct("pk")
        .count()
    )
    print("Inside flank disease variants,", inside_flank_d_variants_num)

    outside_flank_d_variants_num = (
        Variant.objects.filter(residue__flank_residue__feature_location="Outside flank")
        .filter(disease_status="d")
        .distinct("pk")
        .count()
    )
    print("Outside flank disease variants,", outside_flank_d_variants_num)

    tmh_d_variants_clinvar_num = (
        Variant.objects.filter(
            residue__tmh_residue__feature_location="TMH",
            disease_status="d",
            variant_source="ClinVar",
        )
        .distinct("pk")
        .count()
    )
    print("TMH ClinVar disease variants,", tmh_d_variants_clinvar_num)
    flank_d_variants_clinvar_num = (
        Variant.objects.exclude(
            residue__tmh_residue=None, residue__tmh_residue__feature_location="TMH"
        )
        .filter(disease_status="d", variant_source="ClinVar")
        .distinct("pk")
        .count()
    )
    print("Flank ClinVar disease variants,", flank_d_variants_clinvar_num)

    flank_d_variants_humsavar_num = (
        Variant.objects.exclude(
            residue__tmh_residue=None, residue__tmh_residue__feature_location="TMH"
        )
        .filter(disease_status="d", variant_source="Humsavar")
        .distinct("pk")
        .count()
    )
    print("Flank Humsavar disease variants,", flank_d_variants_humsavar_num)

    tmh_d_variants_humsavar_num = (
        Variant.objects.exclude(residue__tmh_residue=None)
        .filter(disease_status="d", variant_source="Humsavar")
        .distinct("pk")
        .count()
    )

    d_variants_humsavar_num = (
        Variant.objects.filter(disease_status="d", variant_source="Humsavar")
        .distinct("pk")
        .count()
    )
    print("Humsavar disease variants,", d_variants_humsavar_num)
    print("TMH Humsavar disease variants,", tmh_d_variants_humsavar_num)
    non_tmh_d_variants_num = (
        Variant.objects.filter(residue__tmh_residue=None, residue__flank_residue=None)
        .filter(disease_status="d")
        .distinct("pk")
        .count()
    )
    print("Non-TMH disease variants,", non_tmh_d_variants_num)
    non_tmh_d_variants_clinvar_num = (
        Variant.objects.filter(
            residue__tmh_residue=None,
            residue__flank_residue=None,
            disease_status="d",
            variant_source="ClinVar",
        )
        .distinct("pk")
        .count()
    )
    print("Non-TMH ClinVar disease variants,", non_tmh_d_variants_clinvar_num)
    non_tmh_d_variants_humsavar_num = (
        Variant.objects.filter(
            residue__tmh_residue=None,
            residue__flank_residue=None,
            disease_status="d",
            variant_source="Humsavar",
        )
        .distinct("pk")
        .count()
    )
    print("Non-TMH Humsavar disease variants,", non_tmh_d_variants_humsavar_num)


def run():
    basic_num()
