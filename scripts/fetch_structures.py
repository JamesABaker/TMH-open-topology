# env LDFLAGS="-I/usr/local/opt/openssl/include -L/usr/local/opt/openssl/lib" pip install_vars psycopg2
from django.db import models

from tmh_db.models import Binding_residue
from tmh_db.models import Database_Metadata
from tmh_db.models import Funfam_residue
from tmh_db.models import Funfamstatus
from tmh_db.models import Go
from tmh_db.models import Keyword
from tmh_db.models import Protein
from tmh_db.models import Residue
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


# This uses Uniprot Keywords to fetch as many pore-containing proteins as possible and returns a list of corresponding pdb ids.
list_of_channel_keyword = [
    "Ion channel",
    "Ion transport",
    "Translocase",
    "Translocation",
    "Sugar transport",
    "Nuclear pore complex",
    "Amino-acid transport",
    "Ammonia transport",
    "Copper transport",
    "Electron transport",
    "ER-Golgi transport",
    "Iron transport",
    "Lipid transport",
    "mRNA transport",
    "Neurotransmitter transport",
    "Peptide transport",
    "Phosphate transport",
    "Protein transport",
    "Sugar transport",
    "Transport",
    "Zinc transport",
]
pdb_ids = (
    Structure.objects.filter(
        uniprot_protein_id__keywords__keyword__in=list_of_channel_keyword,
        uniprot_protein_id__total_tmh_number__gte=1,
    )
    .distinct("pk")
    .values_list("pdb_id")
)

# Just incase of duplicates we list a set of the list.
for i in list(set(pdb_ids)):
    print(i[0])
