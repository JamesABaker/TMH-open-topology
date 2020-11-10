from scripts.populate_general_functions import *
import numpy as np 
import matplotlib.pyplot as plt 


def run():
    tmps=Protein.objects.filter(residue__tmh_residue__tmh_id__meta_tmh=True).distinct('uniprot_id')
    plistmh=[]
    for i in tmps:
        if i.pLI_gn != None:
            plistmh.append(i.pLI_gn)

    nontmp=Protein.objects.exclude(residue__tmh_residue__tmh_id__meta_tmh=True).distinct('uniprot_id')
    plis=[]
    for i in tmps:
        if i.pLI_gn != None:
            plis.append(i.pLI_gn)
    
  
      
    # creating figure and axes to 
    # plot the image 
    fig, (ax1, ax2) = plt.subplots(nrows = 1,  
                                   ncols = 2, 
                                   figsize =(9, 4), 
                                   sharey = True) 
      
    # plotting violin plot for 
    # uniform distribution 
    ax1.set_title('TMP') 
    ax1.set_ylabel('pLI scores') 
    ax1.violinplot(plistmh) 
      
      
    # plotting violin plot for  
    # normal distribution 
    ax2.set_title('Non TMP') 
    ax2.violinplot(plis) 
      
    # Function to show the plot 
    plt.show() 

