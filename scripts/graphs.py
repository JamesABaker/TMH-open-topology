from requests import get
from requests.exceptions import ConnectionError
from datetime import date
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
import collections
from matplotlib.colors import LogNorm

def barchart(objects, performance, source, state, x_label, y_label):
    color_dic = {
        "d": "red",
        "b": "green",
        "a": "grey"
    }
    y_pos = np.arange(len(objects))
    plt.bar(y_pos, performance, color=color_dic[state], align='center', alpha=0.5, edgecolor="grey", width=0.3)
    plt.xticks(y_pos, objects)
    #plt.xlabel('Residue type')
    plt.xlabel(x_label)
    #plt.ylabel('Variants per residue')
    plt.ylabel(y_label)
    plt.title(source)
    filename = f"images/enrichment_{source}.png"
    plt.savefig(filename)
    plt.clf()

def heatmap(var_freqs_list, title, aa_list_baezo_order, state, is_it_log):

    color_dic = {
        "d": "Reds",
        "b": "Greens",
        "a": "Greys",
        "s": "cubehelix"
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
    filename = f"images/{title}.png"
    # And now the colorbar
    # --------------------------------------------------------
    fig.colorbar(im)

    plt.savefig(filename)
    plt.clf()
