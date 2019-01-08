from __future__ import division
import sys
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import datetime
from matplotlib.ticker import MaxNLocator
from collections import namedtuple
from matplotlib import cm
from matplotlib import rcParams


file = str(sys.argv[1])
currentDT = datetime.datetime.now()
# referred to a few time in the functions; this is a list of amino acids in decending order of hydrophobicity according to the Kyte and Doolittle scale.
residues = ["I", "V", "L", "F", "C", "M", "A", "G", "T",
            "S", "W", "Y", "P", "H", "E", "Q", "D", "N", "K", "R"]
### Functions ###


def aa_count(seq):

    aa_type = []
    aa_counts = []
    for aa_all_list in residues:
        aa_all_list_count = 0
        for aa_seq in seq:
            if aa_all_list == aa_seq:
                aa_all_list_count = aa_all_list_count + 1
        aa_type.append(aa_all_list)
        aa_counts.append(aa_all_list_count)
    return(aa_type, aa_counts)


def plot(inside_freqs, tmh_freqs, outside_freqs, y_title):
    n_groups = 20
    fig, ax = plt.subplots()

    index = np.arange(n_groups)
    bar_width = 0.2

    opacity = 0.5
    error_config = {'ecolor': '0.3'}

    rects1 = ax.bar(index, inside_freqs, bar_width,
                    alpha=opacity, label='Inside flank')

    rects2 = ax.bar(index + bar_width, tmh_freqs,
                    bar_width, alpha=opacity, label='TMH')

    rects3 = ax.bar(index + 2 * bar_width, outside_freqs,
                    bar_width, alpha=opacity, label='Outside flank')

    ax.set_xlabel('Residue Type')
    ax.set_ylabel(y_title)
    ax.set_xticks(index + bar_width / 2)
    ax.set_xticklabels(residues)
    ax.legend()

    fig.tight_layout()
    # plt.show()
    plt.savefig(file + y_title + str(currentDT) + ".pdf")
    plt.clf()


def normalise_percentage_total(int_list):
    normalise_percentage_total_list = []
    for i in int_list:
        normalise_percentage_total_list.append((100 / sum(int_list)) * i)
    return(normalise_percentage_total_list)


### Code starts here ###
n_ter_seq_all = str()
in_seq_all = str()
tmh_seq_all = str()
c_ter_seq_all = str()
out_seq_all = str()

results = []

with open(file) as inputfile:
    for line in inputfile:
        results.append(line.strip().split(','))


for entry in results:
    # print(entry)
    # Don't read header line
    if entry == results[0]:
        pass
    else:
        # This list should get bigger as scores etc are added.
        query_id = str(entry[0])
        tmh_start = int(entry[1])
        tmh_stop = int(entry[2])
        tmh_topology = str(entry[3]).strip()
        evidence_type = str(entry[4])
        n_location = str(entry[5])
        n_ter_seq = str(entry[6])
        tmh_seq = str(entry[7])
        c_ter_seq = str(entry[8])

        n_ter_seq_all = str(n_ter_seq + n_ter_seq_all)
        tmh_seq_all = str(tmh_seq + tmh_seq_all)
        c_ter_seq_all = str(c_ter_seq + c_ter_seq_all)

        # Lets sort the flanks into inside outside
        if tmh_topology == str("Inside"):
            in_seq_all = str(n_ter_seq + in_seq_all)
            out_seq_all = str(c_ter_seq + out_seq_all)
        elif tmh_topology == str("Outside"):
            in_seq_all = str(c_ter_seq + in_seq_all)
            out_seq_all = str(n_ter_seq + out_seq_all)
        else:
            #print("Topology missing.")
            pass

inside_freqs = (aa_count(in_seq_all))
tmh_freqs = (aa_count(tmh_seq_all))
outside_freqs = (aa_count(out_seq_all))

print("\nInside flank\n",aa_count(in_seq_all))
print("\nTMH\n",aa_count(tmh_seq_all))
print("\nOutside flank\n"aa_count(out_seq_all))


plot(inside_freqs[1], tmh_freqs[1], outside_freqs[1], str("Residue Frequency"))

plot(normalise_percentage_total(inside_freqs[1]), normalise_percentage_total(tmh_freqs[1]), normalise_percentage_total(
    outside_freqs[1]), str("Percentage of Residues of Feature"))
