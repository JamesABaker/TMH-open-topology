from __future__ import division
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from collections import defaultdict
from tmh_db.models import Binding_residue, Database_Metadata, Flank, Funfam, Funfam_residue, Funfamstatus, Go, Keyword, Flank_residue, Pfam, Pfam_residue, Protein, Residue, Structural_residue, Structure, Subcellular_location, Tmh, Tmh_deltag, Tmh_hydrophobicity, Tmh_residue, Tmh_tmsoc, Uniref, Variant
from django.core.cache import cache
from django.conf import settings
from django.contrib.auth import get_user_model
from django.db import transaction
from django.db.models import Avg, Case, Count, F, Max, Min, Prefetch, Q, Sum, When, Exists, OuterRef, Subquery
from django.utils import timezone
from django.urls import reverse
from scripts.graphs import *
from scripts.populate_general_functions import *
import numpy as np

io_dict={"I": "Inside", "O": "Outside", "Helix":"TMH"}
with open("marcia_list.txt", "r") as f:
    content = f.readlines()

for line, row in enumerate(content):
    if line >1: #head line
        row = row.split()
        #print(row)
        this_uniprot_id=clean_query(row[1])
        this_position=clean_query(row[2])
        base_residue = clean_query(row[3])
        #get command throws errors if residue does not match
        db_residue=Residue.objects.get(protein__uniprot_id=this_uniprot_id, sequence_position=this_position, amino_acid_type=base_residue)

        # tmhs=clean_query(str(list(Tmh.objects.filter(protein__uniprot_id=this_uniprot_id, tmh_residue__residue__sequence_position=this_position, tmh_residue__residue__amino_acid_type=base_residue, tmh_evidence="UniProt").values_list("tm_type"))))
        # flanks=clean_query(str(list(Flank.objects.filter(tmh__protein__uniprot_id=this_uniprot_id, flank_residue__residue__sequence_position=this_position, flank_residue__residue__amino_acid_type=base_residue, tmh__tmh_evidence="UniProt").values_list("inside_or_outside"))))
        tmhs=list(Tmh_residue.objects.filter(tmh_id__protein__uniprot_id=this_uniprot_id, residue__sequence_position=this_position, residue__amino_acid_type=base_residue, tmh_id__tmh_evidence="UniProt").values_list("tmh_id__tm_type", "amino_acid_location_in_to_out"))
        flanks=list(Flank_residue.objects.filter(flank__tmh_id__protein__uniprot_id=this_uniprot_id, residue__sequence_position=this_position, residue__amino_acid_type=base_residue, flank__tmh_id__tmh_evidence="UniProt").values_list("flank__inside_or_outside", "amino_acid_location_in_to_out"))


        if len(flanks)>0:
            flanks=list(flanks[0])
            flanks[0]=io_dict[flanks[0]]
            row.append(flanks)
        elif len(tmhs)>0:
            tmhs=list(tmhs[0])
            tmhs[0]=io_dict[tmhs[0]]
            row.append(tmhs)
        else:
            row.append("Non-TMH")

        print(list_to_csv(row))
