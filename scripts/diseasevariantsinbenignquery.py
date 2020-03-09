# Shell Plus Model Imports
#diseasevariantsinbenignquery
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

Variant.objects.all().prefetch_related("residue", "residue__protein")

disease=Variant.objects.filter(disease_status="d", variant_source="ClinVar").values_list("residue__protein__uniprot_id", "residue__sequence_position", "aa_wt", "aa_mut")
gnomad2=Variant.objects.filter(variant_source="gnomAD2").values_list("residue__protein__uniprot_id", "residue__sequence_position", "aa_wt", "aa_mut")
gnomad3=Variant.objects.filter(variant_source="gnomAD3").values_list("residue__protein__uniprot_id", "residue__sequence_position", "aa_wt", "aa_mut")

for i in disease:
    for x in gnomad2:
        if i==x:
            print(i, "found in gnomAD2 and disease")

for i in disease:
    for x in gnomad3:
        if i==x:
            print(i, "found in gnomAD3 and disease")
