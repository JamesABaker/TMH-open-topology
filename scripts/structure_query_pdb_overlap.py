# Shell Plus Model Imports
from django.contrib.admin.models import LogEntry
from django.contrib.auth.models import Group, Permission, User
from django.contrib.contenttypes.models import ContentType
from django.contrib.sessions.models import Session
from tmh_db.models import Structure 
# Shell Plus Django Imports
from django.core.cache import cache
from django.conf import settings
from django.contrib.auth import get_user_model
from django.db import transaction
from django.db.models import Avg, Case, Count, F, Max, Min, Prefetch, Q, Sum, When
from django.utils import timezone
from django.urls import reverse
from django.db.models import Exists, OuterRef, Subquery


structures=Structure.objects.all().distinct('pk')
matches=0
for i in structures:
    match=False
    #open text file in read mode
    text_file = open("permutated_vars/total_mappable.csv", "r")
 
    #read whole file to a string
    data = text_file.read()
    if i.pdb_id in data:
        match=True 
    #close file
    text_file.close() 
    if match==True:
        matches=matches+1
print(matches)
