import pandas as pd


def create_readstat_table(readstat_list):
    rssum = sum(readstat_list)
    rssum.label = 'Total'

    colnames = [
        'label', 'read_count', 'total_base',
        'read_length_mean',
        'read_length_median',
        'read_length_N50',
        'read_length_L50',
        'mean_insertion_per_read',
        'mean_deletion_per_read',
        'mean_mismatch_per_read',
        'mean_intron_per_read'
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
            a.get_mean_deletions(),
            a.get_mean_mismatches(),
            a.get_mean_introns()
        ])
    # add Total
    N50, L50 = rssum.get_length_NL(50)
    row_list.append([
        rssum.label,
        rssum.get_read_count(),
        rssum.get_total_base(),
        rssum.get_mean_length(),
        rssum.get_median_length(),
        N50, L50,
        rssum.get_mean_insertions(),
        rssum.get_mean_deletions(),
        rssum.get_mean_mismatches(),
        rssum.get_mean_introns()
    ])
    rstable = pd.DataFrame(
        row_list, columns = colnames
    )
    return rstable


def create_mismatch_table(mismatch_list):
    missum = sum(mismatch_list)
    missum.label = 'Total'

    mistypes = list(
        missum.get_count_dict().keys()
    )

    mistable = pd.DataFrame(
        [
            [mismatch_list[i].label] +
            [mismatch_list[i].get_count_dict()[a]
             for a in mistypes]
            for i in range(len(mismatch_list))
        ] + [
            ['Total'] +
            [missum.get_count_dict()[a]
             for a in mistypes]
        ],
        columns = ['label'] + mistypes
    )

    return mistable


def create_insertion_table(insertion_list):
    insum = sum(insertion_list)
    insum.label = 'Total'
    intable = pd.DataFrame(
        [
            [insertion_list[i].label,
             insertion_list[i].get_total_count(),
             insertion_list[i].get_total_insertion_length(),
             insertion_list[i].get_mean_insertion_length(),
             insertion_list[i].get_median_insertion_length()]
            for i in range(len(insertion_list))
        ] + [
            [insum.label,
             insum.get_total_count(),
             insum.get_total_insertion_length(),
             insum.get_mean_insertion_length(),
             insum.get_median_insertion_length()]
        ],
        columns = ['label', 'insertion_count',
                   'insertion_length',
                   'mean_insertion_length',
                   'median_insertion_length']
    )
    return intable


def create_deletion_table(deletion_list):
    desum = sum(deletion_list)
    desum.label = 'Total'
    detable = pd.DataFrame(
        [
            [deletion_list[i].label,
             deletion_list[i].get_total_count(),
             deletion_list[i].get_total_deletion_length(),
             deletion_list[i].get_mean_deletion_length(),
             deletion_list[i].get_median_deletion_length()]
            for i in range(len(deletion_list))
        ] + [
            [desum.label,
             desum.get_total_count(),
             desum.get_total_deletion_length(),
             desum.get_mean_deletion_length(),
             desum.get_median_deletion_length()]
        ],
        columns = ['label', 'deletion_count',
                   'deletion_length',
                   'mean_deletion_length',
                   'median_deletion_length']
    )
    return detable


def create_splice_table(splice_list):
    spsum = sum(splice_list)
    spsum.label = 'Total'

    sptypes = list(
        spsum.get_splice_pair_count_dict().keys()
    )

    sptable = pd.DataFrame(
        [
            [splice_list[i].label] +
            [splice_list[i].get_splice_pair_count_dict()[a]
             for a in sptypes]
            for i in range(len(splice_list))
        ] + [
            ['Total'] +
            [spsum.get_splice_pair_count_dict()[a]
             for a in sptypes]
        ],
        columns = ['label'] + sptypes
    )

    return sptable


########################################
