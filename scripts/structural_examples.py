# Shell Plus Model Imports
from django.contrib.admin.models import LogEntry
from django.contrib.auth.models import Group, Permission, User
from django.contrib.contenttypes.models import ContentType
from django.contrib.sessions.models import Session
from tmh_db.models import Database_Metadata, Disease, Flank, Flank_residue, Funfam, FunfamResidue, Go, Keyword, Non_tmh_helix, Non_tmh_helix_residue, Protein, Residue, Signal_peptide, Signal_residue, Structural_residue, Structure, SubcellularLocation, Tmh, Tmh_deltag, Tmh_hydrophobicity, Tmh_residue, Tmh_tmsoc, Uniref, Variant
# Shell Plus Django Imports
from django.core.cache import cache
from django.conf import settings
from django.contrib.auth import get_user_model
from django.db import transaction
from django.db.models import Avg, Case, Count, F, Max, Min, Prefetch, Q, Sum, When
from django.utils import timezone
from django.urls import reverse
from django.db.models import Exists, OuterRef, Subquery



#Examples



def filter(mut=None, wt=None, variant_query=[], structural_residue_query=[]):
    vars = variant_query.filter(aa_wt=wt, aa_mut=mut)
    residue= structural_residue_query.filter(residue__variant__in=vars)

    print(f"{wt}->{mut}")
    for i in residue:
        structure=Structure.objects.get(structural_residue__pk=i.pk)
        print(structure.pdb_id, i.author_position, i.pdb_chain)


def run():
    alldisvars=Variant.objects.filter(disease_status="d")
    allbenvars=Variant.objects.filter(variant_source="gnomAD3")
    pore = Structural_residue.objects.filter(pore_residue=True, residue__tmh_residue__tmh_id__meta_tmh=True).distinct('pk')
    lipid = Structural_residue.objects.filter(Q(residue__structural_residue__memprotmd_head=True) | Q(residue__structural_residue__memprotmd_tail=True)).filter(residue__tmh_residue__tmh_id__meta_tmh=True).distinct('pk')

    print("Pore disease")
    filter(mut="A", wt="P", variant_query=alldisvars, structural_residue_query=pore)
    print("Lipid Benign")
    filter(mut="A", wt="P", variant_query=allbenvars, structural_residue_query=lipid)

    print("Pore disease")
    filter(mut="L", wt="S", variant_query=alldisvars, structural_residue_query=pore)
    print("Lipid Benign")
    filter(mut="L", wt="S", variant_query=allbenvars, structural_residue_query=lipid)

    print("Pore disease")
    filter(mut="S", wt="L", variant_query=alldisvars, structural_residue_query=pore)
    print("Lipid Benign")
    filter(mut="S", wt="L", variant_query=allbenvars, structural_residue_query=lipid)

    print("Pore benign")
    filter(mut="R", wt="G", variant_query=allbenvars, structural_residue_query=pore)
    print("Lipid disease")
    filter(mut="R", wt="G", variant_query=alldisvars, structural_residue_query=lipid)
