from scripts.populate_general_functions import *

def run():
        #We have to iterate to ensure memory doesn't clap out!
        proteins=Protein.objects.all()
        for i in proteins:

                vars=Variant.objects.filter(residue__protein__uniprot_id=i.uniprot_id)
                print(f'{i.uniprot_id} has {vars.count()} variants to delete...')
                vars.delete()
                print('delete complete')
