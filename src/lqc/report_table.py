import pandas as pd

from lqc.constants import MISMATCH_TYPES, TOTAL_LABEL

COL_LABEL = 'label'
COL_READ_COUNT = 'read_count'


def total_row(table, column):
    """Return the ``column`` value from the ``Total`` row of a summary table."""
    return table.loc[table[COL_LABEL] == TOTAL_LABEL, column].iloc[0]


def create_readstat_table(readstat_list, readstat_sum):
    colnames = [
        COL_LABEL, COL_READ_COUNT, 'total_base', 'aligned_base',
        'read_length_mean', 'read_length_median',
        'read_length_N50', 'read_length_L50',
        'mean_insertion_per_read', 'mean_insertion_per_read_per_kb',
        'insertion_per_query_kb', 'insertion_per_aligned_kb',
        'mean_deletion_per_read', 'mean_deletion_per_read_per_kb',
        'deletion_per_query_kb', 'deletion_per_aligned_kb',
        'mean_mismatch_per_read', 'mean_mismatch_per_read_per_kb',
        'mismatch_per_query_kb', 'mismatch_per_aligned_kb',
        'mean_intron_per_read', 'mean_intron_per_read_per_kb'
    ]

    def _row(a):
        N50, L50 = a.get_length_NL(50)
        return [
            a.label,
            a.get_read_count(),
            a.get_total_base(),
            a.get_total_aligned_base(),
            a.get_mean_length(),
            a.get_median_length(),
            N50, L50,
            a.get_mean_insertions(),
            a.get_mean_length_normalized_insertions() * 1000,
            a.insertions_per_base() * 1000,
            a.insertions_per_aligned_base() * 1000,
            a.get_mean_deletions(),
            a.get_mean_length_normalized_deletions() * 1000,
            a.deletions_per_base() * 1000,
            a.deletions_per_aligned_base() * 1000,
            a.get_mean_mismatches(),
            a.get_mean_length_normalized_mismatches() * 1000,
            a.mismatches_per_base() * 1000,
            a.mismatches_per_aligned_base() * 1000,
            a.get_mean_introns(),
            a.get_mean_length_normalized_introns() * 1000
        ]

    rows = [_row(a) for a in readstat_list] + [_row(readstat_sum)]
    return pd.DataFrame(rows, columns = colnames)


def create_mapping_table(readstat_list, readstat_sum):
    colnames = [
        COL_LABEL, COL_READ_COUNT, 'query_base', 'aligned_base',
        'aligned_fraction_mean', 'aligned_fraction_median',
        'mapq_mean', 'mapq_median',
        'reads_aligned_fraction_lt_0.9', 'reads_fully_aligned'
    ]

    def _row(a):
        return [
            a.label,
            a.get_read_count(),
            a.get_total_base(),
            a.get_total_aligned_base(),
            float(a.get_mean_aligned_fraction()),
            float(a.get_median_aligned_fraction()),
            float(a.get_mean_mapping_quality()),
            float(a.get_median_mapping_quality()),
            a.get_read_count_with_aligned_fraction_below(0.9),
            a.get_read_count_fully_aligned()
        ]

    rows = [_row(a) for a in readstat_list] + [_row(readstat_sum)]
    return pd.DataFrame(rows, columns = colnames)


MISMATCH_TYPES_IN_ORDER = list(MISMATCH_TYPES)


def create_mismatch_normalized_read_location_table(mismatch_list,
                                                   mismatch_sum):
    cuts = [0, 0.1, 0.2, 0.3, 0.4, 0.5,
            0.6, 0.7, 0.8, 0.9, 1]
    mis_type_bin_counts = [
        mismatch_list[i].get_location_bin_count_by_type(cuts = cuts)
        for i in range(len(mismatch_list))
    ]
    sum_type_bin_count = mismatch_sum.get_location_bin_count_by_type(cuts = cuts)
    mistypes = [
        t for t in MISMATCH_TYPES_IN_ORDER
        if t in sum_type_bin_count
    ]
    if not mistypes:
        return pd.DataFrame(columns = [COL_LABEL, 'bin'])
    bins = list(
        sum_type_bin_count[mistypes[0]].keys()
    )

    data_list = [
        [mismatch_list[i].label, ibin] +
        [mis_type_bin_counts[i][c][ibin]
         for c in mistypes]
        for i in range(len(mismatch_list))
        for ibin in bins
    ]
    data_list.extend(
        [mismatch_sum.label, ibin] +
        [sum_type_bin_count[c][ibin]
         for c in mistypes]
        for ibin in bins
    )
    mistable = pd.DataFrame(
        data_list,
        columns = [COL_LABEL, 'bin', *mistypes]
    )
    mistable = mistable.copy()
    mistable['bin_total'] = mistable[mistypes].sum(axis = 1)
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
        columns = [COL_LABEL, 'total_count',
                   'total_length',
                   'mean_length',
                   'median_length']
    )
    return table


def create_splice_table(splice_list, splice_sum):
    """Four-category splice summary (gt-ag / gc-ag / at-ac / other, count+pct)."""

    def _row(a):
        count_dict = a.get_splice_pair_count_dict()
        gtag = count_dict.get('gt-ag', 0)
        gcag = count_dict.get('gc-ag', 0)
        atac = count_dict.get('at-ac', 0)
        other = sum(
            v for k, v in count_dict.items()
            if k not in ('gt-ag', 'gc-ag', 'at-ac')
        )
        total = gtag + gcag + atac + other
        if total == 0:
            gtagp = gcagp = atacp = otherp = 0.0
        else:
            gtagp = gtag / total * 100
            gcagp = gcag / total * 100
            atacp = atac / total * 100
            otherp = other / total * 100
        return [a.label, gtag, gtagp, gcag, gcagp,
                atac, atacp, other, otherp]

    rows = [_row(a) for a in splice_list] + [_row(splice_sum)]
    return pd.DataFrame(
        rows,
        columns = [COL_LABEL, 'gt-ag', 'gt-ag_pct',
                   'gc-ag', 'gc-ag_pct',
                   'at-ac', 'at-ac_pct',
                   'other', 'other_pct']
    )


def create_splice_all_table(splice_list, splice_sum):
    """Full splice-pair matrix (one column per observed pair)."""

    sptypes = list(splice_sum.get_splice_pair_count_dict().keys())
    rows = [
        [splice_list[i].label] +
        [splice_list[i].get_splice_pair_count_dict().get(a, 0)
         for a in sptypes]
        for i in range(len(splice_list))
    ] + [
        [splice_sum.label] +
        [splice_sum.get_splice_pair_count_dict().get(a, 0)
         for a in sptypes]
    ]
    return pd.DataFrame(rows, columns = [COL_LABEL, *sptypes])


########################################
