# Shell Plus Model Imports
from django.contrib.admin.models import LogEntry
from django.contrib.auth.models import Group, Permission, User
from django.contrib.contenttypes.models import ContentType
from django.contrib.sessions.models import Session
from tmh_db.models import Binding_residue, Database_Metadata, Funfam_residue, Funfamstatus, Go, Keyword, Protein, Residue, Structural_residue, Structure, Tmh, Tmh_deltag, Tmh_hydrophobicity, Tmh_residue, Tmh_tmsoc, Variant
# Shell Plus Django Imports
from django.core.cache import cache
from django.conf import settings
from django.contrib.auth import get_user_model
from django.db import transaction
from django.db.models import Avg, Case, Count, F, Max, Min, Prefetch, Q, Sum, When, Exists, OuterRef, Subquery
from django.utils import timezone
from django.urls import reverse
#Charts
import matplotlib.pyplot as plt
import numpy as np


def barchart(objects, performance, source, state):
    color_dic ={
      "d": "red",
      "b": "green",
      "a": "yellow"
    }
    y_pos = np.arange(len(objects))
    plt.bar(y_pos, performance, color = color_dic[state], align='center', alpha=0.5, edgecolor="grey", width = 0.3)
    plt.xticks(y_pos, objects)
    plt.xlabel('Residue type')
    plt.ylabel('Variants per residue')
    plt.title(source)
    filename = f"images/enrichment_{source}.png"
    plt.savefig(filename)
    plt.clf()

def run():
    # Tmh.objects.exclude(residue__tmh_residue=None).filter(disease_status='d').count()
    print("\n\nResidues\n")
    protein_num = Protein.objects.count()
    print("UniProt IDs,", protein_num)
    residue_num = Residue.objects.count()
    print("Residues,", residue_num)
    tmh_residue_num = Residue.objects.exclude(tmh_residue=None).count()
    print("TMH residues,", tmh_residue_num)
    non_tmh_residue_num = Residue.objects.filter(tmh_residue=None).count()
    print("Non-TMH residues,", non_tmh_residue_num)

    print("\n\nTMH boundaries\n")
    tmh_boundary_num = Tmh.objects.count()
    print("All TMH boundaries,", tmh_boundary_num)
    uniprot_tmh_boundary = Tmh.objects.filter(tmh_evidence='UniProt').count()
    print("UniProt TMH boundaries,", uniprot_tmh_boundary)
    topdb_tmh_boundary = Tmh.objects.filter(tmh_evidence='TOPDB').count()
    print("TopDB TMH boundaries,", topdb_tmh_boundary)
    mptopo_tmh_boundary = Tmh.objects.filter(tmh_evidence='MPTOPO').count()
    print("MPTOPO TMH boundaries,", mptopo_tmh_boundary)


    print("\n\nVariants\n")
    # Complex query example. How many variants are in the TMH?
    d_variants_num = Variant.objects.filter(disease_status='d').count()
    print("Disease variants,", d_variants_num)
    tmh_d_variants_num = Variant.objects.exclude(residue__tmh_residue=None).filter(disease_status='d').count()
    print("TMH disease variants,", tmh_d_variants_num)
    d_variants_clinvar_num = Variant.objects.filter(disease_status='d', variant_source="ClinVar").count()
    print("ClinVar disease variants,", d_variants_clinvar_num)
    tmh_d_variants_clinvar_num = Variant.objects.exclude(residue__tmh_residue=None).filter(disease_status='d', variant_source="ClinVar").count()
    print("TMH ClinVar disease variants,", tmh_d_variants_clinvar_num)
    tmh_d_variants_humsavar_num = Variant.objects.exclude(residue__tmh_residue=None).filter(disease_status='d', variant_source="Humsavar").count()
    d_variants_humsavar_num = Variant.objects.filter(disease_status='d', variant_source="Humsavar").count()
    print("Humsavar disease variants,", d_variants_humsavar_num)
    print("TMH Humsavar disease variants,", tmh_d_variants_humsavar_num)
    non_tmh_d_variants_num = Variant.objects.filter(residue__tmh_residue=None).filter(disease_status='d').count()
    print("Non-TMH disease variants,", non_tmh_d_variants_num)
    non_tmh_d_variants_clinvar_num = Variant.objects.filter(residue__tmh_residue=None).filter(disease_status='d', variant_source="ClinVar").count()
    print("Non-TMH ClinVar disease variants,", non_tmh_d_variants_clinvar_num)
    non_tmh_d_variants_humsavar_num = Variant.objects.filter(residue__tmh_residue=None).filter(disease_status='d', variant_source="Humsavar").count()
    print("Non-TMH Humsavar disease variants,", non_tmh_d_variants_humsavar_num)

    print("\n\nVariant enrichment\n")
    print("Disease variants per residue,", d_variants_num/residue_num)
    print("Disease variants per TMH residue,", tmh_d_variants_num/tmh_residue_num)
    print("Disease variants per non-TMH residue,", non_tmh_d_variants_num/non_tmh_residue_num)

    objects = ("Residues","TMH ±5 residues", "Non-TMH residues")
    performance = [d_variants_num/residue_num, tmh_d_variants_num/tmh_residue_num, non_tmh_d_variants_num/non_tmh_residue_num]
    barchart (objects, performance, "All disease variants", "d")

    objects = ("Residues","TMH ±5 residues", "Non-TMH residues")
    performance = [d_variants_clinvar_num/residue_num, tmh_d_variants_clinvar_num/tmh_residue_num, non_tmh_d_variants_clinvar_num/non_tmh_residue_num]
    barchart (objects, performance, "ClinVar disease variants", "d")

    objects = ("Residues","TMH ±5 residues", "Non-TMH residues")
    performance = [d_variants_humsavar_num/residue_num, tmh_d_variants_humsavar_num/tmh_residue_num, non_tmh_d_variants_humsavar_num/non_tmh_residue_num]
    barchart (objects, performance, "Humsavar disease variants", "d")
