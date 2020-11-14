from scripts.populate_general_functions import *

a=Protein.objects.filter(residue__variant__variant_source='gnomAD3').distinct('pk').values_list('uniprot_id')
for i in a:
    print(i[0])
