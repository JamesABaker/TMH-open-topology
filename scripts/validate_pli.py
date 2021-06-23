from scripts.populate_general_functions import *
from scipy import stats
import numpy as np
import matplotlib.pyplot as plt




def run():
    '''
    Django needs this function at the end to run properly
    '''

    # TMP query
    proteins=Protein.objects.all()


    print("UniProt ID, pLI score, Missesnse z-score")
    for i in proteins:
        print(f'{i.uniprot_id}, {i.pLI_gn}, {i.total_tmh_number}')
