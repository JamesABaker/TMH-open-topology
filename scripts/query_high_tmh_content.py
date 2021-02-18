from scripts.populate_general_functions import *
from scripts.graphs import *
import scipy.stats as stats


tmps=Protein.objects.filter(total_tmh_number__gte=1).distinct('pk')

for i in tmps:
  residues=Residue.objects.filter(protein__uniprot_id=i.uniprot_id).distinct('pk')
  tmh_residues=residues.filter(tmh_residue__tmh_id__meta_tmh=True).distinct('pk')
  fraction=tmh_residues.count()/residues.count()
  print(i.uniprot_id, ",", fraction)
