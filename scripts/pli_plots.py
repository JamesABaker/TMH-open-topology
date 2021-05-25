from scripts.populate_general_functions import *
from scipy import stats
import numpy as np
import matplotlib.pyplot as plt


def basic(numbers):
    print("Mean", np.mean(numbers))
    print("Median", np.median(numbers))
    print("Standard deviation", np.std(numbers))


def statistics(list_one=[], list_two=[], threshold=0):
    '''
    compares how many from list one and two scored above and below the
    threshold.
    '''

    list_one_above = []
    list_one_below = []
    list_two_above = []
    list_two_below = []

    for i in list_one:
        if i >= threshold:
            list_one_above.append(i)
        elif i < threshold:
            list_one_below.append(i)

    for i in list_two:
        if i >= threshold:
            list_two_above.append(i)
        elif i < threshold:
            list_two_below.append(i)

    table = [[len(list_one_above), len(list_one_below)], [len(list_two_above), len(list_two_below)]]

    print(table)
    print(stats.fisher_exact(table))
    print(stats.ks_2samp(list_one, list_two))
    print(stats.kruskal(list_one, list_two))
    return()


def violin(list1=[], list1name="List 1", list2=[], list2name="List 2", list3=[], list3name="List 3", title="Numbers"):

    a = list1
    b = list2
    c = list3

    # creating figure and axes to
    # plot the image
    common_params = dict(bins=20,
                         range = (0, 2),
                         density = True)
    # label=[list1name, list2name,list3name])

    plt.subplots_adjust(hspace = .7)
    plt.subplot(311)
    plt.title('Default')
    plt.hist([a, b, c], stacked = True, color = [
     "blue", "red", "green"],  **common_params)
    plt.subplot(312)
    plt.title('Skinny shift - 3 at a time')
    plt.hist((a, b, c), color = ["blue", "red", "green"], **common_params)
    plt.subplot(313)
    common_params['histtype']='step'
    plt.title('With steps')
    plt.hist(a, color = "blue", **common_params)
    plt.hist(b, color = "red", **common_params)
    plt.hist(c, color = "green", **common_params)

    plt.savefig(f'{title}.png')
    plt.cla()
    plt.clf()


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
    '''
    Django needs this function at the end to run properly
    '''

    # TMP query
    tmps=Protein.objects.filter(
        residue__tmh_residue__tmh_id__meta_tmh=True).distinct('uniprot_id')
    plistmh = pli_score_iterator(tmps)
    mistmp = missense_score_iterator(tmps)

    # High Fraction TMP query
    my_file = open('TMP_residues_to_TMH_gte_0.5.txt')
    all_the_lines = my_file.readlines()
    hf_tmp_ids = []
    for i in all_the_lines:
        hf_tmp_ids.append(clean_query(i))
    #print(hf_tmp_ids)

    hf_tmp = Protein.objects.filter(
        residue__tmh_residue__tmh_id__meta_tmh=True, uniprot_id__in=hf_tmp_ids).distinct('uniprot_id')
    plishftmp = pli_score_iterator(hf_tmp)
    mishftmp = missense_score_iterator(hf_tmp)
    # Non TMP query
    nontmp = Protein.objects.exclude(
        residue__tmh_residue__tmh_id__meta_tmh=True).distinct('uniprot_id')
    plis = pli_score_iterator(nontmp)
    mis = missense_score_iterator(nontmp)

    print("pLI")
    print("All TMPS")
    statistics(list_one=plistmh, list_two=plis, threshold=0.9)
    print("High TMH residue fraction TMPS")
    statistics(list_one=plishftmp, list_two=plis, threshold=0.9)
    print("high fraction versus all TMP")
    statistics(list_one=plishftmp, list_two=plistmh, threshold=0.9)
    print("TMP numbers")
    basic(plistmh)
    print("High fraction TMP numbers")
    basic(plishftmp)
    print("Non-tmp numbers")
    basic(plis)
    violin(list1=plistmh, list1name="TMPs", list2=plis, list2name="Non-TMPs",
           list3=plishftmp, list3name="High TMH fraction TMPs", title="pLI scores")

    print("\nMissense")
    print("All TMPS")
    statistics(list_one=mistmp, list_two=mis, threshold=0.35)
    print("High TMH residue fraction TMPS")
    statistics(list_one=mishftmp, list_two=mis, threshold=0.35)
    print("high fraction versus all TMP")
    statistics(list_one=mishftmp, list_two=mistmp, threshold=0.35)

    print("TMP numbers")
    basic(mistmp)
    print("High fraction TMP numbers")
    basic(mishftmp)
    print("Non-tmp numbers")
    basic(mis)
    violin(list1=mistmp, list1name="TMPs", list2=mis, list2name="Non-TMPs", list3=mishftmp,
           list3name="High TMH fraction TMPs", title="Missense upper threshold scores")
