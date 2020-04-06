# Shell Plus Model Imports
#diseasevariantsinbenignquery
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
