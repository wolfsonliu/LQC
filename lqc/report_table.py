import pandas as pd


def create_readstat_table(readstat_list, readstat_sum):
    colnames = [
        'label', 'read_count', 'total_base',
        'read_length_mean',
        'read_length_median',
        'read_length_N50',
        'read_length_L50',
        'mean_insertion_per_read',
        'mean_insertion_per_read_per_kb',
        'insertion_per_kb',
        'mean_deletion_per_read',
        'mean_deletion_per_read_per_kb',
        'deletion_per_kb',
        'mean_mismatch_per_read',
        'mean_mismatch_per_read_per_kb',
        'mismatch_per_kb',
        'mean_intron_per_read',
        'mean_intron_per_read_per_kb'
    ]

    row_list = list()

    for a in readstat_list:
        N50, L50 = a.get_length_NL(50)
        row_list.append([
            a.label,
            a.get_read_count(),
            a.get_total_base(),
            a.get_mean_length(),
            a.get_median_length(),
            N50, L50,
            a.get_mean_insertions(),
            a.get_mean_length_normalized_insertions() * 1000,
            a.insertions_per_base() * 1000,
            a.get_mean_deletions(),
            a.get_mean_length_normalized_deletions() * 1000,
            a.deletions_per_base() * 1000,
            a.get_mean_mismatches(),
            a.get_mean_length_normalized_mismatches() * 1000,
            a.mismatches_per_base() * 1000,
            a.get_mean_introns(),
            a.get_mean_length_normalized_introns() * 1000
        ])
    # add Total
    N50, L50 = readstat_sum.get_length_NL(50)
    row_list.append([
        readstat_sum.label,
        readstat_sum.get_read_count(),
        readstat_sum.get_total_base(),
        readstat_sum.get_mean_length(),
        readstat_sum.get_median_length(),
        N50, L50,
        readstat_sum.get_mean_insertions(),
        readstat_sum.get_mean_length_normalized_insertions() * 1000,
        readstat_sum.insertions_per_base() * 1000,
        readstat_sum.get_mean_deletions(),
        readstat_sum.get_mean_length_normalized_deletions() * 1000,
        readstat_sum.deletions_per_base() * 1000,
        readstat_sum.get_mean_mismatches(),
        readstat_sum.get_mean_length_normalized_mismatches() * 1000,
        readstat_sum.mismatches_per_base() * 1000,
        readstat_sum.get_mean_introns(),
        readstat_sum.get_mean_length_normalized_introns() * 1000
    ])
    rstable = pd.DataFrame(
        row_list, columns = colnames
    )
    return rstable


def create_mismatch_normalized_read_location_table(mismatch_list,
                                                   mismatch_sum):
    cuts = [0, 0.1, 0.2, 0.3, 0.4, 0.5,
            0.6, 0.7, 0.8, 0.9, 1]
    mis_type_bin_counts = [
        mismatch_list[i].get_location_bin_count_by_type(cuts = cuts)
        for i in range(len(mismatch_list))
    ]
    sum_type_bin_count = mismatch_sum.get_location_bin_count_by_type(cuts = cuts)
    mistypes = list(sum_type_bin_count.keys())
    bins = list(
        sum_type_bin_count[mistypes[0]].keys()
    )

    data_list = list()
    for i in range(len(mismatch_list)):
        for ibin in bins:
            data_list.append(
                [mismatch_list[i].label, ibin] +
                [mis_type_bin_counts[i][c][ibin]
                 for c in mistypes]
            )
    for ibin in bins:
        data_list.append(
            [mismatch_sum.label, ibin] +
            [sum_type_bin_count[c][ibin]
             for c in mistypes]
        )
    mistable = pd.DataFrame(
        data_list,
        columns = ['label', 'bin'] + mistypes
    )
    return mistable


def create_indel_summary_table(indel_list, indel_sum):
    table = pd.DataFrame(
        [
            [indel_list[i].label,
             indel_list[i].get_total_count(),
             indel_list[i].get_total_length(),
             indel_list[i].get_mean_length(),
             indel_list[i].get_median_length()]
            for i in range(len(indel_list))
        ] + [
            [indel_sum.label,
             indel_sum.get_total_count(),
             indel_sum.get_total_length(),
             indel_sum.get_mean_length(),
             indel_sum.get_median_length()]
        ],
        columns = ['label', 'total_count',
                   'total_length',
                   'mean_length',
                   'median_length']
    )
    return table


def create_splice_table(splice_list, splice_sum):

    sptypes = list(
        splice_sum.get_splice_pair_count_dict().keys()
    )

    sptable = pd.DataFrame(
        [
            [splice_list[i].label] +
            [splice_list[i].get_splice_pair_count_dict()[a]
             for a in sptypes]
            for i in range(len(splice_list))
        ] + [
            ['Total'] +
            [splice_sum.get_splice_pair_count_dict()[a]
             for a in sptypes]
        ],
        columns = ['label'] + sptypes
    )

    return sptable


########################################
