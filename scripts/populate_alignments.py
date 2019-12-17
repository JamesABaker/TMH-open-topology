# Shell Plus Model Imports
from django.contrib.admin.models import LogEntry
from django.contrib.auth.models import Group, Permission, User
from django.contrib.contenttypes.models import ContentType
from django.contrib.sessions.models import Session
# Shell Plus Django Imports
from django.core.cache import cache
from django.conf import settings
from django.contrib.auth import get_user_model
from django.db import transaction
from django.db.models import Avg, Case, Count, F, Max, Min, Prefetch, Q, Sum, When, Exists, OuterRef, Subquery
from django.utils import timezone
from django.urls import reverse
# General scripts
from scripts.populate_general_functions import *

from Bio import AlignIO


protein_query_set=Protein.objects.all().distinct("pk").values_list("uniprot_id")


for ids in protein_query_set
target_seq=Protein.objects.get(uniprot_id=)
ref_seq=

inp =
outp =
cline = MuscleCommandline(input=inp, out=outp)
cline()
