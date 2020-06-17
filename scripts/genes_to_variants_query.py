import time
import urllib.parse
import urllib.request

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

from tmh_db.models import Database_Metadata, Flank, Flank_residue, Funfam, Funfam_residue, Go, Keyword, Non_tmh_helix, Non_tmh_helix_residue, Protein, Residue, Signal_peptide, Signal_residue, Structural_residue, Structure, SubcellularLocation, Tmh, Tmh_deltag, Tmh_hydrophobicity, Tmh_residue, Tmh_tmsoc, Uniref, Variant

# Shell Plus Model Imports
# Shell Plus Django Imports


def uniprot_query_to_list(uniprot_result):
    lines = uniprot_result.splitlines()
    uniprot_ids=[]
    for line_number, line_content in enumerate(lines):
        if line_number==0:
            pass
        elif line_number==1:
            pass
        else:
            items=line_content.split("\t")
            #print(items[1])
            if items[1] == None:
                pass
            else:
                uniprot_ids.append(items[1])
    return(uniprot_ids)

def gene_to_unirptot_id(text):
    no_download=True

    while no_download==True:
        try:
            url = 'https://www.uniprot.org/uploadlists/'

            params = {
                #'from': 'GENENAME',
                'from': 'GENENAME',
                'to': 'ACC',
                'format': 'tab',
                'query': text
            }

            data = urllib.parse.urlencode(params)
            data = data.encode('utf-8')
            req = urllib.request.Request(url, data)
            with urllib.request.urlopen(req) as f:
                response = f.read()
            uniprot_id = response.decode('utf-8')
            #print(uniprot_id)
            no_download=False
            return(uniprot_id)
        except:
            print("OOPS we clogged up UniProt!")
            time.sleep(5)


def variant_query(id):
    variants=Variant.objects.filter(residue__protein__uniprot_id=id).exclude(disease_status="d").count()
    d_variants=Variant.objects.filter(residue__protein__uniprot_id=id, disease_status="d").count()
    residues=Residue.objects.filter(protein__uniprot_id=id).count()
    #values_list("aa_wt", "aa_mut", "residue__sequence_position", "residue__protein__uniprot_id", "variant_source", "disease_status")
    return(str(str(variants)+", "+str(d_variants) + ", "+str(residues)))

def run():
    text_list="covid_queries/supp_table_2_pedro_paper.txt"
    with open(text_list, "r") as f:
        lines = f.read().splitlines()
        for line_text in lines:
            if "#" in line_text:
                pass
            else:
                uniprot_result = gene_to_unirptot_id(line_text)
                uniprot_ids=uniprot_query_to_list(uniprot_result)
                #print(uniprot_ids)
                for uniprot_id in uniprot_ids:

                    all_variants = variant_query(uniprot_id)
                    print(line_text, ",", uniprot_id, ",", all_variants)
