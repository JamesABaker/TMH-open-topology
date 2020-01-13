from django.conf import settings
from django.db import models
from datetime import datetime, timedelta
from django.utils import timezone
from datetime import date
import pytz
from scripts.populate_general_functions import *

def overlap(tmh_start, tmh_stop, comparison_start, comparison_stop):
    if tmh_start >= comparison_start and tmh_start <=comparison_stop:
        overlap=True
    elif tmh_stop <= comparison_stop and tmh_stop >= comparison_start:
        overlap=True
    else:
        overlap=False
    return(overlap)



def priority(evidence, comparison_evidence):
    priority={"UniProt":1, "TOPDB":2, "MPTOPO":3, "OPM":4}
    if priority[evidence] > priority[comparison_evidence]:
        return(True)
    else:
        return(False)

# Canonical script starts here
# Fetch TMHs from the database
Tmh.objects.all().prefetch_related("protein")
tmhs=Tmh.objects.filter().distinct("pk")
tmh_list=[]

for i in tmhs:
    tmh_start=i.tmh_start
    tmh_stop=i.tmh_stop
    tmh_protein=i.protein.uniprot_id
    tmh_evidence=i.tmh_evidence

    tmhs_for_comparison=[]
    for x in tmhs:
        comparison_start=x.tmh_start
        comparison_stop=x.tmh_stop
        comparison_protein=x.protein.uniprot_id
        comparison_evidence=x.tmh_evidence
        if tmh_protein == comparison_protein:
            if tmh_evidence == comparison_evidence:
                pass
            else:
                is_overlap = overlap(tmh_start, tmh_stop, comparison_start, comparison_stop)

                if is_overlap == True:
                    tmhs_for_comparison.append(x)

    meta=True
    rep=None
    for comparison_tmh in tmhs_for_comparison:
        comparison_evidence=comparison_tmh.tmh_evidence
        is_priority = priority(tmh_evidence, comparison_tmh)
        if is_priotity == False:
            meta=False
            rep=comparison_tmh
            #i.meta_tmh_rep=comparison_tmh
        else:
            rep=None

        i.meta_tmh=meta

bulk_update(tmhs)
