from scripts.populate_general_functions import *
from scripts.graphs import *
from django.db.models import F
import scipy.stats as stats
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

Variant.objects.all().prefetch_related("residue", "residue__tmh_residue__tmh_id")
Residue.objects.all().prefetch_related("tmh_residue__tmh_id")
Tmh_deltag.objects.all().prefetch_related("tmh_id")


def get_tmhs():
    '''
    Gets all the meta TMHs.
    '''
    tmh_objects=Tmh.objects.filter(meta_tmh=True).distinct('pk')
    return(tmh_objects)

def get_scores(tmh_objects=[]):
    '''
    extracts the number of residues, variants, and delta g from the tmh objects
    '''
    data_returned=[]
    for i in tmh_objects:
        residues=Residue.objects.filter(tmh_residue__tmh_id__pk=i.pk).distinct('pk').count()
        disease_variants=Variant.objects.filter(disease_status="d", residue__tmh_residue__tmh_id__pk=i.pk).distinct('pk').count()
        #benign_variants=Variant.objects.filter(variant_source="gnomAD3", residue__tmh_residue__tmh_id__pk=i.pk).distinct('pk').count()
        spont_data=Tmh_deltag.objects.get(tmh_id__pk=i.pk)
        spont=spont_data.test_score
        data_returned.append([np.divide(disease_variants, residues), spont])
    return(data_returned)

def between(list1,low,high):
    list2 = []
    for i in list1:
        if(i > low and i < high) and i != 0:
            list2.append(i)
    return list2

def plot_data(values=[]):
    '''
    plots the values as a scatter plot separating a list of lists to x and y
    [[1,5],[2,5]] -> x=[1,2], y=[5,5]
    '''
    x=[]
    y=[]
    for i in values:
        if i[0] != 0:
            x.append(i[1])
            y.append(i[0])
    avg_in_limit=[]
    avg_index=[]

    for i in np.linspace(-5, 5, 0.1):
        print(i)
        in_limits=between(x, i-0.05, i+0.05)
        average_in_limit=np.mean(in_limits)
        avg_in_limit.append(average_in_limit)
        avg_index.append(i)

    plt.xlim([-5.5, 5.5])
    plt.scatter(x, y, alpha=0.05)
    plt.scatter( avg_in_limit, avg_index, alpha=1, color="red")
    plt.title('Scatter plot pythonspot.com')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.show()

    plt.clf()



def run():
    '''
    The main loop
    '''
    tmhs=get_tmhs()
    data=get_scores(tmh_objects=tmhs)
    plot_data(values=data)

    print("complete")
