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

disvars=Variant.objects.filter(disease_status="d", aa_wt="I", aa_mut="N")
benvars=Variant.objects.filter(variant_source="gnomAD3", aa_wt="I", aa_mut="N")

struc_res=Structural_residue.objects.filter(pore_residue=True, residue__variant__in=disvars).distinct('pk')
print("I->N Pore disease")
for i in struc_res:

    structure=Structure.objects.get(structural_residue__pk=i.pk)
    print(structure.pdb_id, i.pdb_position, i.pdb_chain)


struc_res=Structural_residue.objects.filter(memprotmd_tail=True, residue__variant__in=benvars).distinct('pk')
print("I->N Lipid contact benign")
for i in struc_res:

    structure=Structure.objects.get(structural_residue__pk=i.pk)
    print(structure.pdb_id, i.pdb_position, i.pdb_chain)




disvars=Variant.objects.filter(disease_status="d", aa_wt="S", aa_mut="L")
benvars=Variant.objects.filter(variant_source="gnomAD3", aa_wt="S", aa_mut="L")

struc_res=Structural_residue.objects.filter(pore_residue=True, residue__variant__in=disvars).distinct('pk')
print("S->L Pore disease")
for i in struc_res:
    structure=Structure.objects.get(structural_residue__pk=i.pk)
    print(structure.pdb_id, i.pdb_position, i.pdb_chain)


struc_res=Structural_residue.objects.filter(memprotmd_head=True, residue__variant__in=benvars).distinct('pk')
print("S->L Lipid contact benign")
for i in struc_res:
    structure=Structure.objects.get(structural_residue__pk=i.pk)
    print(structure.pdb_id, i.pdb_position, i.pdb_chain)






disvars=Variant.objects.filter(disease_status="d", aa_wt="G", aa_mut="R")
benvars=Variant.objects.filter(variant_source="gnomAD3", aa_wt="G", aa_mut="R")

struc_res=Structural_residue.objects.filter(pore_residue=True, residue__variant__in=benvars).distinct('pk')
print("G->R Pore benign")
for i in struc_res:

    structure=Structure.objects.get(structural_residue__pk=i.pk)
    print(structure.pdb_id, i.pdb_position, i.pdb_chain)


struc_res=Structural_residue.objects.filter(memprotmd_tail=True, residue__variant__in=disvars).distinct('pk')
print("G->R Lipid contact disease")
for i in struc_res:
    structure=Structure.objects.get(structural_residue__pk=i.pk)
    print(structure.pdb_id, i.pdb_position, i.pdb_chain)
