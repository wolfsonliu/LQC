import math
from itertools import accumulate
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use("Agg")

magentas = ["#7a0177", "#ae017e", "#dd3497",
            "#f768a1", "#fa9fb5", "#fcc5c0",
            "#fde0dd", "#fff7f3"]


def get_facet_row_col(number):
    row = int(math.sqrt(number))
    col = row
    while row * col < number:
        col += 1
    return row, col


def determine_figure_size(row, col,
                          base_width = 4,
                          base_height = 3):
    return col * base_width, row * base_height


def plot_readstat_bar(readstat_list,
                      feature,
                      width = None,
                      height = None):
    assert feature in [
        "Read count",
        "Median read length",
        "Mean read length",
        "Insertions per read",
        "Insertions per read per kb",
        "Deletions per read",
        "Deletions per read per kb",
        "Mismatches per read",
        "Mismatches per read per kb",
        "Mean intron number",
        "N50", "L50"
    ]
    if feature == "Read count":
        read_feature = [
            [i, a.label, a.get_read_count()]
            for i, a in enumerate(readstat_list)
        ]
    elif feature == "Median read length":
        read_feature = [
            [i, a.label, a.get_median_length()]
            for i, a in enumerate(readstat_list)
        ]
    elif feature == "Mean read length":
        read_feature = [
            [i, a.label, a.get_mean_length()]
            for i, a in enumerate(readstat_list)
        ]
    elif feature == "Insertions per read":
        read_feature = [
            [i, a.label, a.get_mean_insertions()]
            for i, a in enumerate(readstat_list)
        ]
    elif feature == "Insertions per read per kb":
        read_feature = [
            [i, a.label, a.get_mean_length_normalized_insertions()]
            for i, a in enumerate(readstat_list)
        ]
    elif feature == "Deletions per read":
        read_feature = [
            [i, a.label, a.get_mean_deletions()]
            for i, a in enumerate(readstat_list)
        ]
    elif feature == "Deletions per read per kb":
        read_feature = [
            [i, a.label, a.get_mean_length_normalized_deletions()]
            for i, a in enumerate(readstat_list)
        ]
    elif feature == "Mismatches per read":
        read_feature = [
            [i, a.label, a.get_mean_mismatches()]
            for i, a in enumerate(readstat_list)
        ]
    elif feature == "Mismatches per read per kb":
        read_feature = [
            [i, a.label, a.get_mean_length_normalized_mismatches()]
            for i, a in enumerate(readstat_list)
        ]
    elif feature == "Mean intron number":
        read_feature = [
            [i, a.label, a.get_mean_introns()]
            for i, a in enumerate(readstat_list)
        ]
    elif feature == "N50":
        read_feature = [
            [i, a.label, a.get_length_NL(50)[0]]
            for i, a in enumerate(readstat_list)
        ]
    elif feature == "L50":
        read_feature = [
            [i, a.label, a.get_length_NL(50)[1]]
            for i, a in enumerate(readstat_list)
        ]
    else:
        pass

    if width is None or height is None:
        width = max(5, len(read_feature) * 0.5)
        height = 4
    else:
        pass
    fig, axes = plt.subplots(
        1, 1, figsize = (width, height)
    )
    axes.bar(
        [read_feature[i][0]
         for i in range(len(read_feature))],
        [read_feature[i][2]
         for i in range(len(read_feature))],
        width = 0.66,
        color = ["#7a0177", "#dd3497"],
        fill = True
    )
    axes.set_xticks(
        [read_feature[i][0]
         for i in range(len(read_feature))]
    )
    axes.set_xticklabels(
        [read_feature[i][1]
         for i in range(len(read_feature))]
    )
    axes.set_xlabel("Contig")
    axes.set_ylabel(feature)
    plt.tight_layout()
    return fig


def plot_readstat_bar_mean_element_per_read(readstat_list,
                                            width = None,
                                            height = None):
    read_feature = [
        [a.get_mean_insertions(),
         a.get_mean_deletions(),
         a.get_mean_mismatches(),
         a.get_mean_introns(),
         i,
         a.label]
        for i, a in enumerate(readstat_list)
    ]

    row, col = get_facet_row_col(
        len(read_feature)
    )

    if width is None or height is None:
        width, height = determine_figure_size(
            row, col,
            base_width = 2, base_height = 2
        )
        width = max(width, 5)
        height = max(height, 5)
    else:
        pass

    colors = ["#ef3b2c", "#f16913",
              "#4292c6", "#41ab5d"]
    names = ["Insertion", "Deletion",
             "Mismatch", "Intron"]
    binw = 0.22
    fig, axes = plt.subplots(
        row, col,
        sharex = True, sharey = True,
        figsize = (width, height)
    )
    for ai in range(row * col):
        if ai in list(range(len(read_feature))):
            for ti in [0, 1, 2, 3]:
                fig.axes[ai].bar(
                    [(ti * 2 - 3) * 0.5 * binw],
                    [read_feature[ai][ti]],
                    width = binw - 0.01,
                    color = colors[ti],
                    fill = True,
                    label = names[ti]
                )
                fig.axes[ai].set_xticks([])
                fig.axes[ai].set_title(
                    read_feature[ai][5]
                )
        else:
            fig.axes[ai].set_frame_on(False)
            fig.axes[ai].set_axis_off()

    Bar, Label = fig.axes[0].get_legend_handles_labels()
    fig.legend(
        Bar, Label,
        ncol = len(names),
        loc = 'lower center',
        bbox_to_anchor = (0.5, 0)
    )
    fig.supylabel("Mean count per read")
    plt.tight_layout(rect = [0, 0.05, 1, 1])
    return fig


def plot_readstat_bar_mean_element_per_read_per_kb(
        readstat_list, width = None, height = None
):
    read_feature = [
        [a.get_mean_length_normalized_insertions() * 1000,
         a.get_mean_length_normalized_deletions() * 1000,
         a.get_mean_length_normalized_mismatches() * 1000,
         a.get_mean_length_normalized_introns() * 1000,
         i,
         a.label]
        for i, a in enumerate(readstat_list)
    ]

    row, col = get_facet_row_col(
        len(read_feature)
    )

    if width is None or height is None:
        width, height = determine_figure_size(
            row, col,
            base_width = 2, base_height = 2
        )
        width = max(width, 5)
        height = max(height, 5)
    else:
        pass

    colors = ["#ef3b2c", "#f16913",
              "#4292c6", "#41ab5d"]
    names = ["Insertion", "Deletion",
             "Mismatch", "Intron"]
    binw = 0.22
    fig, axes = plt.subplots(
        row, col,
        sharex = True, sharey = True,
        figsize = (width, height)
    )
    for ai in range(row * col):
        if ai in list(range(len(read_feature))):
            for ti in [0, 1, 2, 3]:
                fig.axes[ai].bar(
                    [(ti * 2 - 3) * 0.5 * binw],
                    [read_feature[ai][ti]],
                    width = binw - 0.01,
                    color = colors[ti],
                    fill = True,
                    label = names[ti]
                )
                fig.axes[ai].set_xticks([])
                fig.axes[ai].set_title(
                    read_feature[ai][5]
                )
        else:
            fig.axes[ai].set_frame_on(False)
            fig.axes[ai].set_axis_off()

    Bar, Label = fig.axes[0].get_legend_handles_labels()
    fig.legend(
        Bar, Label,
        ncol = len(names),
        loc = 'lower center',
        bbox_to_anchor = (0.5, 0)
    )
    fig.supylabel("Mean count per read per kb")
    plt.tight_layout(rect = [0, 0.05, 1, 1])
    return fig


def plot_readstat_bar_ratio_with_element(readstat_list,
                                         width = None,
                                         height = None):
    read_feature = [
        [1 -
         a.get_read_count_with_n_insertions(0) /
         a.get_read_count(),
         1 -
         a.get_read_count_with_n_deletions(0) /
         a.get_read_count(),
         1 -
         a.get_read_count_with_n_mismatches(0) /
         a.get_read_count(),
         i,
         a.label]
        for i, a in enumerate(readstat_list)
    ]

    row, col = get_facet_row_col(
        len(read_feature)
    )

    if width is None or height is None:
        width, height = determine_figure_size(
            row, col,
            base_width = 3, base_height = 3
        )
        width = max(width, 5)
        height = max(height, 5)
    else:
        pass

    colors = ["#ef3b2c", "#f16913",
              "#4292c6"]
    names = ["Insertion", "Deletion",
             "Mismatch"]
    binw = 0.3
    fig, axes = plt.subplots(
        row, col,
        sharex = True, sharey = True,
        figsize = (width, height)
    )
    for ai in range(row * col):
        if ai in list(range(len(read_feature))):
            for ti in [0, 1, 2]:
                fig.axes[ai].bar(
                    [(ti - 1) * binw],
                    [read_feature[ai][ti]],
                    width = binw - 0.01,
                    color = colors[ti],
                    fill = True,
                    label = names[ti]
                )
                ylimit = fig.axes[ai].get_ylim()
                fig.axes[ai].annotate(
                    '{:.2}'.format(
                        read_feature[ai][ti]
                    ),
                    ((ti - 1) * binw,
                     read_feature[ai][ti] - 0.1 * ylimit[1]),
                    color = 'white',
                    ha = "center", va = "top"
                )
                fig.axes[ai].set_xticks([])
                fig.axes[ai].set_title(
                    read_feature[ai][4]
                )
        else:
            fig.axes[ai].set_frame_on(False)
            fig.axes[ai].set_axis_off()

    Bar, Label = fig.axes[0].get_legend_handles_labels()
    fig.legend(
        Bar, Label,
        ncol = len(names),
        loc = 'lower center',
        bbox_to_anchor = (0.5, 0)
    )
    fig.supylabel("Ratio of reads with elements")
    plt.tight_layout(rect = [0, 0.05, 1, 1])
    return fig


def plot_readstat_cumulative_length(readstat_list,
                                    width = None,
                                    height = None):
    row, col = get_facet_row_col(
        len(readstat_list)
    )

    if width is None or height is None:
        width, height = determine_figure_size(
            row, col,
            base_width = 3, base_height = 3
        )
    else:
        pass

    fig, axes = plt.subplots(
        row, col, figsize = (width, height)
    )
    for ai in range(row * col):
        if ai in list(range(len(readstat_list))):
            lengths = sorted(
                readstat_list[ai].get_lengths()
            )
            length_cumsum = list(
                accumulate([a for a in lengths])
            )
            N50, L50 = readstat_list[ai].get_length_NL(50)
            # plot
            fig.axes[ai].plot(
                lengths, length_cumsum,
                color = "#969696",
                linestyle = '-'
            )
            cross_x = lengths[-1 * L50]
            cross_y = length_cumsum[-1 * L50]
            xlimits = fig.axes[ai].get_xlim()
            xspan = xlimits[1] - xlimits[0]
            ylimits = fig.axes[ai].get_ylim()
            yspan = ylimits[1] - ylimits[0]
            fig.axes[ai].axvline(
                cross_x,
                ymin = 0, ymax = 1,
                color = "#c51b8a", linestyle = ':'
            )
            fig.axes[ai].axhline(
                cross_y,
                xmin = 0, xmax = 1,
                color = "#c51b8a", linestyle = ':'
            )
            fig.axes[ai].plot(
                cross_x, cross_y,
                color = "#7a0177", marker = "."
            )
            fig.axes[ai].annotate(
                "\n".join(
                    ["N50: {}".format(N50),
                     "L50: {}".format(L50)]
                ),
                (cross_x, cross_y),
                (cross_x + 0.02 * xspan,
                 cross_y - 0.02 * yspan),
                ha = "left", va = "top"
            )
            fig.axes[ai].set_title(
                readstat_list[ai].label
            )
        else:
            fig.axes[ai].set_frame_on(False)
            fig.axes[ai].set_axis_off()
    fig.supxlabel("Read length")
    fig.supylabel("Cumulative bases")
    plt.tight_layout()
    return fig


def plot_readstat_length_hist(readstat_list,
                              width = None,
                              height = None):

    row, col = get_facet_row_col(
        len(readstat_list)
    )
    if width is None and height is None:
        width, height = determine_figure_size(
            row, col,
            base_width = 3, base_height = 2
        )
        width = max(width, 5)
        height = max(height, 5)
    else:
        pass
    fig, axes = plt.subplots(
        row, col, figsize = (width, height)
    )
    for ai in range(row * col):
        if ai in list(range(len(readstat_list))):
            fig.axes[ai].hist(
                readstat_list[ai].get_lengths(),
                color = "#7a0177",
                label = readstat_list[ai].label
            )
            fig.axes[ai].set_title(
                readstat_list[ai].label
            )
            xlimits = fig.axes[ai].get_xlim()
            ylimits = fig.axes[ai].get_ylim()
            fig.axes[ai].text(
                xlimits[1] * 0.98,
                ylimits[1] * 0.98,
                "Median:\n{}".format(int(
                    readstat_list[ai].get_median_length()
                )),
                ha = "right", va = "top"
            )
        else:
            fig.axes[ai].set_frame_on(False)
            fig.axes[ai].set_axis_off()

    fig.supxlabel("Read length")
    fig.supylabel("Read count")
    plt.tight_layout()
    return fig


def plot_element_total_count(element_list,
                             element_name,
                             width = None,
                             height = None):
    assert element_name in [
        "Insertion", "Deletion", "Mismatch", "Intron"
    ], "element_name should be one of [Insertion, Deletion, Mismatch, Intron]."
    if element_name == "Intron":
        count_list = [
            [i, a.label,
             a.get_total_splice_pair_count()]
            for i, a in enumerate(element_list)
        ]
    else:
        count_list = [
            [i, a.label,
             a.get_total_count()]
            for i, a in enumerate(element_list)
        ]
    if width is None or height is None:
        width = max(5, len(count_list) * 0.5)
        height = 4
    else:
        pass
    fig, axes = plt.subplots(
        1, 1, figsize = (width, height)
    )
    axes.bar(
        [count_list[i][0]
         for i in range(len(count_list))],
        [count_list[i][2]
         for i in range(len(count_list))],
        width = 0.66,
        color = ["#7a0177", "#dd3497"],
        fill = True
    )
    axes.set_xticks(
        [count_list[i][0]
         for i in range(len(count_list))]
    )
    axes.set_xticklabels(
        [count_list[i][1]
         for i in range(len(count_list))]
    )
    axes.set_xlabel("Contig")
    axes.set_ylabel("Count")
    axes.set_title(element_name)
    plt.tight_layout()
    return fig


def plot_indel_hist_length(indel_list,
                           width = None,
                           height = None):

    row, col = get_facet_row_col(
        len(indel_list)
    )
    if width is None and height is None:
        width, height = determine_figure_size(
            row, col,
            base_width = 3, base_height = 2
        )
        width = max(width, 5)
        height = max(height, 5)
    else:
        pass
    fig, axes = plt.subplots(
        row, col, figsize = (width, height)
    )
    for ai in range(row * col):
        if ai in list(range(len(indel_list))):
            fig.axes[ai].hist(
                indel_list[ai].get_lengths(),
                color = "#7a0177",
                label = indel_list[ai].label
            )
            fig.axes[ai].set_title(
                indel_list[ai].label
            )
            xlimits = fig.axes[ai].get_xlim()
            ylimits = fig.axes[ai].get_ylim()
            fig.axes[ai].text(
                xlimits[1] * 0.98,
                ylimits[1] * 0.98,
                "Median:\n{}".format(int(
                    indel_list[ai].get_median_length()
                )),
                ha = "right", va = "top"
            )
        else:
            fig.axes[ai].set_frame_on(False)
            fig.axes[ai].set_axis_off()

    fig.supxlabel("Element length")
    fig.supylabel("Element count")
    plt.tight_layout()
    return fig


def plot_indel_hist_location(indel_list,
                             width = None,
                             height = None):

    row, col = get_facet_row_col(
        len(indel_list)
    )
    if width is None and height is None:
        width, height = determine_figure_size(
            row, col,
            base_width = 3, base_height = 2
        )
        width = max(width, 5)
        height = max(height, 5)
    else:
        pass
    fig, axes = plt.subplots(
        row, col, figsize = (width, height)
    )
    for ai in range(row * col):
        if ai in list(range(len(indel_list))):
            fig.axes[ai].hist(
                indel_list[ai].get_locations(),
                color = "#7a0177",
                label = indel_list[ai].label
            )
            fig.axes[ai].set_title(
                indel_list[ai].label
            )
        else:
            fig.axes[ai].set_frame_on(False)
            fig.axes[ai].set_axis_off()

    fig.supxlabel("Normalized read location")
    fig.supylabel("Element count")
    plt.tight_layout()
    return fig


def plot_mismatch_type_count(mismatch_list,
                             width = None,
                             height = None):
    miscolors = {
        "ac": "#a50f15", "ag": "#ef3b2c",
        "at": "#fc9272", "ca": "#54278f",
        "cg": "#807dba", "ct": "#bcbddc",
        "ga": "#08519c", "gc": "#4292c6",
        "gt": "#9ecae1", "ta": "#006d2c",
        "tc": "#41ab5d", "tg": "#a1d99b"
    }
    mistypes = ["ac", "ag", "at",
                "ca", "cg", "ct",
                "ga", "gc", "gt",
                "ta", "tc", "tg"]
    row, col = get_facet_row_col(
        len(mismatch_list) + 1
    )

    if width is None and height is None:
        width, height = determine_figure_size(
            row, col,
            base_width = 2.5, base_height = 2
        )
        width = max(width, 5)
        height = max(height, 4)
    else:
        pass

    mislist = list()

    for a in mismatch_list:
        alist = list()
        for mt in mistypes:
            alist.append(a.get_type_count()[mt])
        alist.append(a.label)
        mislist.append(alist)

    fig, axes = plt.subplots(
        row, col, figsize = (width, height)
    )
    for ai in range(row * col):
        if ai in list(range(len(mismatch_list))):
            for j in range(len(mistypes)):
                fig.axes[ai].bar(
                    j, mislist[ai][j],
                    color = miscolors[mistypes[j]],
                    label = '{}=>{}'.format(
                        mistypes[j][0],
                        mistypes[j][1]
                    )
                )
            fig.axes[ai].set_xticks([])
            fig.axes[ai].set_title(
                mislist[ai][len(mistypes)]
            )
        else:
            if ai == len(mismatch_list):
                Bar, Label = fig.axes[0].get_legend_handles_labels()
                fig.axes[ai].legend(
                    Bar, Label,
                    ncol = 2,
                    loc = 'center',
                    title = "Mismatch types"
                )
            else:
                pass
            fig.axes[ai].set_frame_on(False)
            fig.axes[ai].set_axis_off()
    fig.supylabel("Count of mismatches")
    plt.tight_layout()
    return fig


def plot_mismatch_hist_location(mismatch_list,
                                width = None,
                                height = None):

    row, col = get_facet_row_col(
        len(mismatch_list)
    )
    if width is None and height is None:
        width, height = determine_figure_size(
            row, col,
            base_width = 3, base_height = 2
        )
        width = max(width, 5)
        height = max(height, 5)
    else:
        pass
    fig, axes = plt.subplots(
        row, col, figsize = (width, height)
    )
    for ai in range(row * col):
        if ai in list(range(len(mismatch_list))):
            fig.axes[ai].hist(
                mismatch_list[ai].get_locations(),
                color = "#7a0177",
                label = mismatch_list[ai].label
            )
            fig.axes[ai].set_title(
                mismatch_list[ai].label
            )
        else:
            fig.axes[ai].set_frame_on(False)
            fig.axes[ai].set_axis_off()

    fig.supxlabel("Normalized read location")
    fig.supylabel("Mismatch count")
    plt.tight_layout()
    return fig


def plot_splice_type_count(splice_list,
                           width = None,
                           height = None):
    splicecolors = {
        "gt-ag": "#7a0177", "gc-ag": "#dd3497",
        "at-ac": "#f768a1", "other": "#969696"
    }
    splicetypes = ["gt-ag", "gc-ag",
                   "at-ac", "other"]
    row, col = get_facet_row_col(
        len(splice_list) + 1
    )
    if width is None and height is None:
        width, height = determine_figure_size(
            row, col,
            base_width = 2.5, base_height = 2
        )
        width = max(width, 5)
        height = max(height, 4)
    else:
        pass

    splist = list()

    for a in splice_list:
        alist = list()
        count_dict = a.get_splice_pair_count_dict()
        for s in splicetypes[0:3]:
            alist.append(count_dict[s])
        other = 0
        for s in count_dict:
            if s not in splicetypes[0:3]:
                other += count_dict[s]
            else:
                pass
        alist.append(other)
        alist.append(a.label)
        splist.append(alist)

    fig, axes = plt.subplots(
        row, col, figsize = (width, height)
    )
    for ai in range(row * col):
        if ai in list(range(len(splice_list))):
            for j in range(len(splicetypes)):
                st = splicetypes[j]
                fig.axes[ai].bar(
                    j, splist[ai][j],
                    color = splicecolors[st],
                    label = splicetypes[j]
                )
                if st == 'other':
                    stidx = len(splicetypes) - 1
                    ratio = splist[ai][stidx] / \
                        sum(splist[ai][0:stidx])
                    fig.axes[ai].text(
                        j, splist[ai][j],
                        '{:.2}'.format(ratio),
                        ha = "center", va = "bottom"
                    )
            fig.axes[ai].set_xticks(
                list(range(len(splicetypes)))
            )
            fig.axes[ai].set_xticklabels(splicetypes)
            fig.axes[ai].set_title(
                splist[ai][len(splicetypes)]
            )
        else:
            if ai == len(splice_list):
                Bar, Label = fig.axes[0].get_legend_handles_labels()
                fig.axes[ai].legend(
                    Bar, Label,
                    ncol = 1,
                    loc = 'center',
                    title = "Splice types"
                )
            else:
                pass
            fig.axes[ai].set_frame_on(False)
            fig.axes[ai].set_axis_off()
    fig.supylabel("Count of splice sites")
    plt.tight_layout()
    return fig


########################################
