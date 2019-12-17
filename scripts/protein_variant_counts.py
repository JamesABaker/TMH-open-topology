# Shell Plus Model Imports
from django.contrib.admin.models import LogEntry
from django.contrib.auth.models import Group, Permission, User
from django.contrib.contenttypes.models import ContentType
from django.contrib.sessions.models import Session
from tmh_db.models import Binding_residue, Database_Metadata, Flank, Flank_residue, Funfam, Funfam_residue, Funfamstatus, Go, Keyword, Non_tmh_helix, Non_tmh_helix_residue, Pfam, Pfam_residue, Phmmer_proteins, Phmmer_residues, Protein, Residue, Signal_peptide, Signal_residue, Structural_residue, Structure, Subcellular_location, Tail_anchor, Tmh, Tmh_deltag, Tmh_hydrophobicity, Tmh_residue, Tmh_tmsoc, Uniref, Variant
# Shell Plus Django Imports
from django.core.cache import cache
from django.conf import settings
from django.contrib.auth import get_user_model
from django.db import transaction
from django.db.models import Avg, Case, Count, F, Max, Min, Prefetch, Q, Sum, When, Exists, OuterRef, Subquery
from django.utils import timezone
from django.urls import reverse

a=Protein.objects.all().distinct("pk").values_list("uniprot_id", "full_sequence")
for i in a:
    a_uniprot_id=i[0]
    sequence=len(i[1])
    number_of_disease_variant=Variant.objects.filter(residue__protein__uniprot_id=a_uniprot_id, disease_status="d").distinct("pk").count()
    number_of_benign_variant=Variant.objects.filter(residue__protein__uniprot_id=a_uniprot_id).exclude(disease_status="d").distinct("pk").count()
    print(a_uniprot_id, sequence, number_of_disease_variant, number_of_benign_variant)
