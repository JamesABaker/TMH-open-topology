from scripts.populate_general_functions import *
from scipy import stats
import numpy as np 
import matplotlib.pyplot as plt 

def statistics(list1=[], list2=[]):
    print(stats.ks_2samp(list1, list2))
    print(stats.kruskal(list1, list2))

def basic(numbers):
    print("Mean", np.mean(numbers))
    print("Median", np.median(numbers))
    print("Standard deviation", np.std(numbers))

def violin(list1=[], list1name="List 1", list2=[], list2name="List 2", title="Numbers"):


    # creating figure and axes to 
    # plot the image 
    fig, (ax1, ax2) = plt.subplots(nrows = 1,  
                                   ncols = 2, 
                                   figsize =(9, 4), 
                                   sharey = True) 
      
    # plotting violin plot for 
    # uniform distribution 
    ax1.set_title(list1name) 
    ax1.set_ylabel(title) 
    ax1.hist(list1) 
      
      
    # plotting violin plot for  
    # normal distribution 
    ax2.set_title(list2name)
    ax2.hist(list2) 
      
    # Function to show the plot 
    plt.show() 

def pli_score_iterator(queryset):
    scorelist=[]
    for i in queryset:
        if i.pLI_gn != None:
           scorelist.append(i.pLI_gn) 
    return(scorelist)

def missense_score_iterator(queryset):
    scorelist=[]
    for i in queryset:
        if i.oe_mis_upper_gn != None:
           scorelist.append(i.oe_mis_upper_gn) 
    return(scorelist)

def run():
    ''' Django needs this function at the end to run properly'''
    #TMP query
    tmps=Protein.objects.filter(residue__tmh_residue__tmh_id__meta_tmh=True).distinct('uniprot_id')
    plistmh=pli_score_iterator(tmps)
    mistmp=missense_score_iterator(tmps)
    #Non TMP query
    nontmp=Protein.objects.exclude(residue__tmh_residue__tmh_id__meta_tmh=True).distinct('uniprot_id')
    plis=pli_score_iterator(nontmp)
    mis=missense_score_iterator(nontmp)
    
    print("pLI")
    statistics(plistmh, plis)
    print("TMP numbers")
    basic(plistmh)
    print("Non-tmp numbers")
    basic(plis)  
    violin(list1=plistmh, list1name="TMPs", list2=plis, list2name="Non-TMPs", title="pLI scores")  
    
    print("\nMissense")
    statistics(mistmp, mis)
    print("TMP numbers")
    basic(mistmp)
    print("Non-tmp numbers")
    basic(mis)  
    violin(list1=mistmp, list1name="TMPs", list2=mis, list2name="Non-TMPs", title="Missense upper threshold scores")  
