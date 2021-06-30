from scripts.populate_general_functions import *
from scipy import stats
import numpy as np
import matplotlib.pyplot as plt


def basic(numbers):
    print("Mean", np.mean(numbers))
    print("Median", np.median(numbers))
    print("Standard deviation", np.std(numbers))


def statistics(list_one=[], list_two=[]):
    '''
    compares how many from list one and two scored above and below the
    threshold.
    '''

    odds, pvalue=stats.fisher_exact([list_one, list_two])
    #print(stats.ks_2samp(list_one, list_two))
    #print(stats.kruskal(list_one, list_two))
    return((odds, pvalue))


def pli_score_iterator(queryset):

    scorelist=[]
    for i in queryset:
        if i.pLI_gn != None:
            scorelist.append(i.pLI_gn)
    return(scorelist)



def pli_threshold(pli_scores):
    above=[]
    below=[]
    for i in pli_scores:
        if i >= 0.9:
            above.append(i)
        elif i < 0.9:
            below.append(i)
    return((len(above), len(below)))


def query_to_scores(queryset_of_proteins):
    '''
    Turns the query set into a row of mean, above, below.
    '''
    plis = pli_score_iterator(queryset_of_proteins)
    pli_average=np.mean(plis)
    pli_sd=np.std(plis)
    pli_above_below=pli_threshold(plis)
    return([pli_average, pli_above_below[0], pli_above_below[1], pli_sd])



def run():
    '''
    Django needs this function at the end to run properly
    '''
    print("Proteins, Average pLI, Tolerant, Intolerant, Fisher test against non-TMPs")

    # Non TMP query
    ref_query = Protein.objects.exclude(residue__tmh_residue__tmh_id__meta_tmh=True).distinct('uniprot_id')
    ref_scores=query_to_scores(ref_query)
    print(f'Non-TMPs,{ref_scores[0]},{ref_scores[1]},{ref_scores[2]}, {statistics([ref_scores[1], ref_scores[2]],[ref_scores[1], ref_scores[2]])}, {ref_scores[3]}')


    # Different TMPs

    query=Protein.objects.filter(total_tmh_number__gt=1).distinct('uniprot_id')
    scores=query_to_scores(query)
    print(f'Multipass,{scores[0]},{scores[1]},{scores[2]}, {statistics([scores[1], scores[2]],[ref_scores[1], ref_scores[2]])}, {scores[3]}')

    query=Protein.objects.filter(total_tmh_number=1).distinct('uniprot_id')
    scores=query_to_scores(query)
    print(f'Singlepass,{scores[0]},{scores[1]},{scores[2]}, {statistics([scores[1], scores[2]],[ref_scores[1], ref_scores[2]])}, {scores[3]}')

    query=Protein.objects.filter(total_tmh_number__gte=1, keywords__keyword="G-protein coupled receptor").distinct('uniprot_id')
    scores=query_to_scores(query)
    print(f'GPCRs,{scores[0]},{scores[1]},{scores[2]}, {statistics([scores[1], scores[2]],[ref_scores[1], ref_scores[2]])}, {scores[3]}')

    query=Protein.objects.filter(total_tmh_number__gte=1, keywords__keyword="Ion channel").distinct('uniprot_id')
    scores=query_to_scores(query)
    print(f'Ion channel,{scores[0]},{scores[1]},{scores[2]}, {statistics([scores[1], scores[2]],[ref_scores[1], ref_scores[2]])},  {scores[3]}')


    query=Protein.objects.filter(keywords__keyword="Transport").exclude(total_tmh_number__gte=1, keywords__keyword="Ion channel").distinct('uniprot_id')
    scores=query_to_scores(query)
    print(f'Transport excluding ion channels,{scores[0]},{scores[1]},{scores[2]}, {statistics([scores[1], scores[2]],[ref_scores[1], ref_scores[2]])},  {scores[3]}')

    query=Protein.objects.filter(keywords__keyword="Kinase").distinct('uniprot_id')
    scores=query_to_scores(query)
    print(f'Kinase,{scores[0]},{scores[1]},{scores[2]}, {statistics([scores[1], scores[2]],[ref_scores[1], ref_scores[2]])},  {scores[3]}')
