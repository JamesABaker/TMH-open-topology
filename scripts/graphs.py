from requests import get
from datetime import date
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
import collections
from datetime import datetime
date = datetime.now().strftime('%Y-%m-%d %H:%M:%S')


impossible_snps={"K": ["M", "T", "E","Q","R","N"],
                "P":["T","Q","R","H"],
                "H":["L","Y","D","Q","N","P"],
                "N":["I","Y","T","D","H","K"],
                "G":["W","V","Y","T","E","D"],
                "S":[],
                "R":[],
                "Q":[],
                "E":[],
                "T":[],
                "A":[],
                "Y":[],
                "C":[],
                "V":[],
                "M":[],
                "I":[],
                "L":[],
                "W":[],
                "F":[],
                }


def barchart(objects, performance, source, state, x_label, y_label):
    color_dic = {
        "d": "red",
        "n": "green",
        "a": "grey"
    }
    y_pos = np.arange(len(objects))
    plt.bar(y_pos, performance, color=color_dic[state], align='center', alpha=0.5, edgecolor="grey", width=0.3)
    plt.xticks(y_pos, objects)
    #plt.xlabel('Residue type')
    plt.xlabel(x_label)
    plt.xticks(rotation=45)
    #plt.ylabel('Variants per residue')
    plt.ylabel(y_label)
    plt.title(source)
    filename = f"images/{source}_{date}.png"
    plt.savefig(filename, bbox_inches='tight')
    plt.clf()

def histogram(performance, source, state,x_label, y_label):
    color_dic = {
        "d": "red",
        "n": "green",
        "a": "grey"
    }
    plt.yscale('log', nonposy='clip')
    plt.hist(performance, bins=len(performance))
    plt.gca().set(title=source, ylabel=y_label, xlabel=x_label);
    plt.title(source)
    filename = f"images/{source}_{date}.png"
    plt.savefig(filename, bbox_inches='tight')
    plt.clf()

def heatmap(var_freqs_list, title, aa_list_baezo_order, state, is_it_log):

    color_dic = {
        "d": "Reds",
        "n": "Greens",
        "a": "Greys",
        "s": "afmhot"
    }

    fig, ax = plt.subplots()
    im = ax.imshow(var_freqs_list, cmap=color_dic[state], interpolation='nearest', norm=is_it_log)

    # We want to show all ticks...
    ax.set_xticks(np.arange(len(aa_list_baezo_order)))
    ax.set_yticks(np.arange(len(aa_list_baezo_order)))
    # ... and label them with the respective list entries
    ax.set_xticklabels(aa_list_baezo_order)
    ax.set_yticklabels(aa_list_baezo_order)
    # Add boxes
    ax.grid(which='minor', color='b', linestyle='-', linewidth=2)
    ax.set_xlabel("Wildtype (from)")
    ax.set_ylabel("Variant (to)")
    ax.set_title(title)



    fig.tight_layout()
    filename = f"images/{title}_{date}.png"
    # And now the colorbar
    # --------------------------------------------------------
    fig.colorbar(im)

    for x_coordinate, x_amino_acid in enumerate(aa_list_baezo_order):
        for y_coordinate, y_amino_acid in enumerate(aa_list_baezo_order):
            plt.Circle((x_coordinate, y_coordinate), 0.5, color='black', fill=False)

    plt.savefig(filename)
    plt.clf()