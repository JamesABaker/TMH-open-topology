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
from django.db.models import Avg, Case, Count, F, Max, Min, Prefetch, Q, Sum, When, Exists, OuterRef, Subquery
from django.utils import timezone
from django.urls import reverse
from django.db.models import Prefetch


Variant.objects.all().prefetch_related("residue", "residue__funfamresidue",
                                       "residue__funfamresidue__funfam", "residue__protein")
#variants=Variant.objects.filter(residue__protein__uniprot_id="Q8WZ75").distinct('pk')
variants = Variant.objects.all().distinct('pk')
for n, var_object in enumerate(variants):
    #print(f"{n} of {variants.count()} variants")
    origin_id = var_object.residue.protein.uniprot_id
    origin_position = var_object.residue.sequence_position
    variant_source = var_object.variant_source
    variant_status = var_object.disease_status
    variant_id = var_object.variant_source_id
    funfam_info = []
    for i in FunfamResidue.objects.filter(residue=var_object.residue):
        funfam_id = i.funfam.funfam_id
        superfam_id = i.funfam.superfamily
        funfam_position = i.funfam_position
        funfam_details = [funfam_id, superfam_id, funfam_position]
        funfam_info.append(funfam_details)

    print(f"{origin_id}, {origin_position}, {variant_source}, {variant_status}, {variant_id}, {funfam_info}")
