import collections
from datetime import date
from datetime import datetime

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
from requests import get

date = datetime.now().strftime("%Y-%m-%d %H:%M:%S")


def barchart(objects, performance, source, state, x_label, y_label):
    color_dic = {"d": "red", "n": "green", "a": "grey"}
    y_pos = np.arange(len(objects))
    plt.bar(
        y_pos,
        performance,
        color=color_dic[state],
        align="center",
        alpha=0.5,
        edgecolor="grey",
        width=0.3,
    )
    plt.xticks(y_pos, objects)
    # plt.xlabel('Residue type')
    plt.xlabel(x_label)
    plt.xticks(rotation=45)
    # plt.ylabel('Variants per residue')
    plt.ylabel(y_label)
    plt.title(source)
    filename = f"images/{source}_{date}.png"
    plt.savefig(filename, bbox_inches="tight")
    plt.clf()


def histogram(performance, source, state, x_label, y_label):

    plt.yscale("log", nonposy="clip")
    plt.hist(performance, bins=len(performance))
    plt.gca().set(title=source, ylabel=y_label, xlabel=x_label)
    plt.title(source)
    filename = f"images/{source}_{date}.png"
    plt.savefig(filename, bbox_inches="tight")
    plt.clf()


def impossible_cordinates(aa_order):

    impossible_substitutions = scripts.populate_general_functions.impossible_subs()
    impossible_aa = []
    impossible_x_coordinates = []
    impossible_y_coorindates = []
    for x_coordinate, x_amino_acid in enumerate(aa_order):
        for y_coordinate, y_amino_acid in enumerate(aa_order):
            if y_amino_acid in impossible_substitutions[x_amino_acid]:
                impossible_x_coordinates.append(x_coordinate)
                impossible_y_coorindates.append(y_coordinate)
    return (impossible_x_coordinates, impossible_y_coorindates)


def x_histogram(var_freqs_list):
    x_values = []
    for column_number, frequencies in enumerate(var_freqs_list):
        x_coord_total = []
        for row in var_freqs_list:
            x_coord_total.append(row[column_number])
        x_values.append(sum(x_coord_total))
    return x_values


def y_histogram(var_freqs_list):
    y_values = []
    for row_number, frequencies in enumerate(var_freqs_list):
        total_residues_in_row = 0
        for i in var_freqs_list[row_number]:
            total_residues_in_row = total_residues_in_row + i
        y_values.append(total_residues_in_row)
    y_values = y_values[::-1]
    return y_values


def annotate_heatmap(
    im,
    data=None,
    valfmt="{x:.2f}",
    textcolors=["black", "white"],
    threshold=None,
    **textkw,
):
    """
    A function to annotate a heatmap.

    Parameters
    ----------
    im
        The AxesImage to be labeled.
    data
        Data used to annotate.  If None, the image's data is used.  Optional.
    valfmt
        The format of the annotations inside the heatmap.  This should either
        use the string format method, e.g. "$ {x:.2f}", or be a
        `matplotlib.ticker.Formatter`.  Optional.
    textcolors
        A list or array of two color specifications.  The first is used for
        values below a threshold, the second for those above.  Optional.
    threshold
        Value in data units according to which the colors from textcolors are
        applied.  If None (the default) uses the middle of the colormap as
        separation.  Optional.
    **kwargs
        All other arguments are forwarded to each call to `text` used to create
        the text labels.
    """

    if not isinstance(data, (list, np.ndarray)):
        data = im.get_array()

    # Normalize the threshold to the images color range.
    if threshold is not None:
        threshold = im.norm(threshold)
    else:
        threshold = im.norm(data.max()) / 2.0

    # Set default alignment to center, but allow it to be
    # overwritten by textkw.
    kw = dict(horizontalalignment="center", verticalalignment="center")
    kw.update(textkw)

    # Get the formatter in case a string is supplied
    if isinstance(valfmt, str):
        valfmt = matplotlib.ticker.StrMethodFormatter(valfmt)

    # Loop over the data and create a `Text` for each "pixel".
    # Change the text's color depending on the data.
    texts = []
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            kw.update(color=textcolors[int(im.norm(data[i, j]) > threshold)])
            text = im.axes.text(j, i, valfmt(data[i, j], None), **kw)
            texts.append(text)

    return texts


def heatmap(
    var_freqs_list,
    title,
    aa_order,
    color,
    is_it_log,
    annotation_format="{x:n}",
    bars=True,
):
    """
    var_freqs_list
        heatmap array in list of lists
    title
        title of data for the graph
    aa_list_baezo_order
        ordered list of amino acids for the list of lists
    is_it_log
        is the scale logarithmic? Accepts None.
    """

    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    spacing = 0.01

    rect_heatmap = [left, bottom, width, height]
    rect_histx = [left, bottom + height + spacing, width, 0.2]
    rect_histy = [left + width + spacing, bottom, 0.2, height]
    rect_scale = [left, bottom, width, height]

    # fig, ax = plt.subplots()
    plt.figure(figsize=(10, 10), dpi=80)

    ax_heatmap = plt.axes(rect_heatmap)
    ax_heatmap.tick_params(direction="in", top=True, right=True)
    ax_histy = plt.axes(rect_histy, xticklabels=[], yticklabels=[])
    ax_histy.tick_params(direction="in", labelleft=False)
    ax_histx = plt.axes(rect_histx, xticklabels=[], yticklabels=[])
    ax_histx.tick_params(direction="in", labelbottom=False)
    ax_scale = plt.axes(rect_scale)

    # print(var_freqs_list)
    im = ax_heatmap.imshow(
        var_freqs_list, cmap=color, interpolation="nearest", norm=is_it_log
    )
    texts = annotate_heatmap(im, valfmt=annotation_format, fontsize=20)

    # We want to show all ticks...
    ax_heatmap.set_xticks(np.arange(len(aa_order)))
    ax_heatmap.set_yticks(np.arange(len(aa_order)))
    # ... and label them with the respective list entries
    ax_heatmap.set_xticklabels(aa_order, fontsize=20)
    ax_heatmap.set_yticklabels(aa_order, fontsize=20)
    # Add boxes
    ax_heatmap.grid(which="minor", color="b", linestyle="-", linewidth=2)
    ax_heatmap.set_xlabel("Wildtype (from)", fontsize=20)
    ax_heatmap.set_ylabel("Variant (to)", fontsize=20)
    # ax.set_ylim(len(aa_order)-0.5, -0.5) # temp bug workaround: https://github.com/matplotlib/matplotlib/issues/14751
    # ax.set_xlim(len(aa_order)-0, -0.5) # temp bug workaround: https://github.com/matplotlib/matplotlib/issues/14751
    # ax_heatmap.set_title(title)

    flagged_cells = impossible_cordinates(aa_order)
    ax_heatmap.scatter(
        flagged_cells[0], flagged_cells[1], marker="x", color="silver", s=75
    )

    ax_histy.axis("off")
    ax_histx.axis("off")
    if bars == True:
        # histogram on the right
        ax_histx.bar(
            list(range(0, 20)),
            x_histogram(var_freqs_list),
            orientation="vertical",
            color="grey",
        )

        ax_histx.axis("off")
        ax_histx.margins(0.00)
        rects = ax_histx.patches

        # Make some labels.
        labels = x_histogram(var_freqs_list)

        for rect, label in zip(rects, labels):
            height = rect.get_height()
            ax_histx.text(
                rect.get_x() + rect.get_width() / 2,
                height + 5,
                int(label),
                ha="center",
                va="bottom",
            )

        # ax_histx.hist(x_histogram(var_freqs_list), 21, histtype='stepfilled', orientation='vertical', color='grey')

        # histogram in the top
        # print(list(range(0, 20)), y_histogram(var_freqs_list))
        # ax_histy.bar(list(range(0, 20)), y_histogram(var_freqs_list), orientation='horizontal', color='grey')
        y_values_ordered = y_histogram(var_freqs_list).reverse()
        ax_histy.barh(
            list(range(0, 20)),
            y_histogram(var_freqs_list),
            orientation="horizontal",
            color="grey",
        )

        ax_histy.margins(0.00)
        for i, v in enumerate(y_histogram(var_freqs_list)):
            ax_histy.text(v + 3, i, int(v))
        # ax_histy.barh(list(range(0, 20)), y_values_ordered, color='gainsboro')

    # ax_histx.set_xlim(ax_heatmap.get_xlim())
    # ax_histy.set_ylim(ax_heatmap.get_ylim())
    # And now the colorbar
    # --------------------------------------------------------
    ax_scale = plt.colorbar(im)

    # plt.tight_layout()
    filename = f"images/{title}_{date}.png"

    for x_coordinate, x_amino_acid in enumerate(aa_order):
        for y_coordinate, y_amino_acid in enumerate(aa_order):
            plt.Circle((x_coordinate, y_coordinate), 0.5, color="black", fill=False)

    # plt.tight_layout()
    plt.savefig(filename)
    plt.clf()
    plt.close()
