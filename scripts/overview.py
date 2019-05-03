from django.contrib.admin.models import LogEntry
from django.contrib.auth.models import Group, Permission, User
from django.contrib.contenttypes.models import ContentType
from django.contrib.sessions.models import Session
from tmh_db.models import Binding_residue, Database_Metadata, Funfam_residue, Funfamstatus, Go, Keyword, Protein, Residue, Structural_residue, Structure, Tmh, Tmh_deltag, Tmh_hydrophobicity, Tmh_residue, Tmh_tmsoc, Variant
from django.core.cache import cache
from django.conf import settings
from django.contrib.auth import get_user_model
from django.db import transaction
from django.db.models import Avg, Case, Count, F, Max, Min, Prefetch, Q, Sum, When, Exists, OuterRef, Subquery
from django.utils import timezone
from django.urls import reverse

print("STATS OVERVIEW\n",)

#General
print("Total UniProt IDs,",)
print("Total residues,",)
print("Total TMH residues,",)
print("Total UniProt TMH residues,",)

#structure
print("Number of UniProt IDs with structures,",)
print("Number of structures,",)
print("Total residues from structures,",)
print("Total TMH residues from structures,",)

#topology
print("Total known topologies TOPDB,",)
print("Total known topologies UniProt,",)

#variants
print("Total variants,",)
print("Total disease variants,",)
print("Total TMH variants,",)
print("Total TMH+flanks variants,",)
print("Total TMH disease variants,",)
print("Total TMH+flanks disease variants,",)
print("Total TMH+flanks disease variants in binding residues,",)
print("Total TMH+flanks disease variants in catalytic residues,",)
