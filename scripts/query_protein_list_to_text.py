
from scripts.populate_general_functions import *

# Shell Plus Django Imports

def run():
        proteins=Protein.objects.all().distinct('pk')

        for i in proteins:
                print(i.uniprot_id)
