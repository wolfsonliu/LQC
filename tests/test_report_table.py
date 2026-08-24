import pytest

from lqc import Indel, Mismatch, Splice
from lqc.readstat import ReadStat
from lqc.report_table import (
    create_indel_summary_table,
    create_mapping_table,
    create_mismatch_normalized_read_location_table,
    create_readstat_table,
    create_splice_all_table,
    create_splice_table,
)

READSTAT_COLUMNS = [
    'label', 'read_count', 'total_base', 'aligned_base',
    'read_length_mean', 'read_length_median',
    'read_length_N50', 'read_length_L50',
    'mean_insertion_per_read', 'mean_insertion_per_read_per_kb',
    'insertion_per_query_kb', 'insertion_per_aligned_kb',
    'mean_deletion_per_read', 'mean_deletion_per_read_per_kb',
    'deletion_per_query_kb', 'deletion_per_aligned_kb',
    'mean_mismatch_per_read', 'mean_mismatch_per_read_per_kb',
    'mismatch_per_query_kb', 'mismatch_per_aligned_kb',
    'mean_intron_per_read', 'mean_intron_per_read_per_kb',
]


def _two_readstat():
    r1 = ReadStat('chr1')
    r1.add_read(200, insertion=2, deletion=4, mismatch=6, intron=8,
                mapping_quality=60, aligned_length=190)
    r2 = ReadStat('chr2')
    r2.add_read(200, insertion=4, deletion=8, mismatch=12, intron=16,
                mapping_quality=60, aligned_length=190)
    total = sum([r1, r2])
    total.label = 'Total'
    return r1, r2, total


def test_readstat_table_column_contract():
    r1, r2, total = _two_readstat()
    df = create_readstat_table([r1, r2], total)
    assert list(df.columns) == READSTAT_COLUMNS


def test_readstat_table_values():
    r1, r2, total = _two_readstat()
    df = create_readstat_table([r1, r2], total)
    assert list(df['label']) == ['chr1', 'chr2', 'Total']
    total_row = df[df['label'] == 'Total'].iloc[0]
    assert total_row['read_count'] == 2
    assert total_row['total_base'] == 400
    assert total_row['mean_mismatch_per_read'] == pytest.approx(9)
    assert total_row['mean_mismatch_per_read_per_kb'] == pytest.approx(45)
    assert total_row['mean_intron_per_read'] == pytest.approx(12)


def test_readstat_table_aligned_columns():
    r1, r2, total = _two_readstat()
    df = create_readstat_table([r1, r2], total)
    total_row = df[df['label'] == 'Total'].iloc[0]
    # 2 reads x aligned 190 = 380 aligned bases
    assert total_row['aligned_base'] == 380
    assert total_row['insertion_per_query_kb'] == pytest.approx(6 / 400 * 1000)
    assert total_row['insertion_per_aligned_kb'] == pytest.approx(6 / 380 * 1000)
    assert total_row['mismatch_per_aligned_kb'] == pytest.approx(18 / 380 * 1000)


def test_indel_summary_table():
    i1 = Indel('chr1')
    i1.add_indel('a', 0.1)
    i2 = Indel('chr2')
    i2.add_indel('tt', 0.5)
    total = i1 + i2
    total.label = 'Total'
    df = create_indel_summary_table([i1, i2], total)
    assert list(df.columns) == [
        'label', 'total_count', 'total_length', 'mean_length', 'median_length'
    ]
    assert list(df['label']) == ['chr1', 'chr2', 'Total']
    assert total.get_total_count() == 2


def test_splice_table():
    s1 = Splice('chr1')
    s1.add_splice_pair('gt-ag')
    s2 = Splice('chr2')
    s2.add_splice_pair('gt-ag')
    total = s1 + s2
    total.label = 'Total'

    df = create_splice_table([s1, s2], total)
    assert list(df.columns) == [
        'label', 'gt-ag', 'gt-ag_pct', 'gc-ag', 'gc-ag_pct',
        'at-ac', 'at-ac_pct', 'other', 'other_pct',
    ]
    total_row = df[df['label'] == 'Total'].iloc[0]
    assert total_row['gt-ag'] == 2
    assert total_row['gt-ag_pct'] == 100.0
    assert total_row['other'] == 0
    assert total_row['other_pct'] == 0.0


def test_splice_all_table_preserves_full_matrix():
    s = Splice('chr1')
    s.add_splice_pair_count_dict({'gt-ag': 3, 'gc-ag': 1, 'tt-tt': 4})
    total = Splice('Total')
    total.add_splice_pair_count_dict({'gt-ag': 3, 'gc-ag': 1, 'tt-tt': 4})

    df = create_splice_all_table([s], total)
    assert list(df.columns) == ['label', 'gt-ag', 'gc-ag', 'tt-tt']
    total_row = df[df['label'] == 'Total'].iloc[0]
    assert total_row['gt-ag'] == 3
    assert total_row['tt-tt'] == 4


def test_splice_table_other_bucket_and_percent():
    s1 = Splice('chr1')
    s1.add_splice_pair_count_dict({'gt-ag': 3, 'gc-ag': 1, 'tt-tt': 4})
    s2 = Splice('chr2')
    s2.add_splice_pair_count_dict({'gt-ag': 1})
    total = s1 + s2
    total.label = 'Total'

    df = create_splice_table([s1, s2], total)
    total_row = df[df['label'] == 'Total'].iloc[0]
    assert total_row['gt-ag'] == 4
    assert total_row['gc-ag'] == 1
    assert total_row['at-ac'] == 0
    assert total_row['other'] == 4
    assert total_row['other_pct'] == pytest.approx(4 / 9 * 100)
    assert str(df['other_pct'].dtype) == 'float64'


def test_splice_all_table_fills_missing_pairs_with_zero():
    s = Splice('chr1')
    s.add_splice_pair_count_dict({'gt-ag': 3})
    total = Splice('Total')
    total.add_splice_pair_count_dict({'gt-ag': 3, 'gc-ag': 1})

    df = create_splice_all_table([s], total)
    assert list(df.columns) == ['label', 'gt-ag', 'gc-ag']
    assert df[df['label'] == 'chr1'].iloc[0]['gc-ag'] == 0
    assert df[df['label'] == 'Total'].iloc[0]['gc-ag'] == 1


def test_mapping_table_columns_and_values():
    r1 = ReadStat('chr1')
    r1.add_read(200, insertion=2, deletion=4, mismatch=6, intron=8,
                mapping_quality=60, aligned_length=180)
    r2 = ReadStat('chr2')
    r2.add_read(100, insertion=0, deletion=0, mismatch=0, intron=0,
                mapping_quality=30, aligned_length=100)
    total = sum([r1, r2])
    total.label = 'Total'
    df = create_mapping_table([r1, r2], total)

    assert list(df.columns) == [
        'label', 'read_count', 'query_base', 'aligned_base',
        'aligned_fraction_mean', 'aligned_fraction_median',
        'mapq_mean', 'mapq_median',
        'reads_aligned_fraction_lt_0.9', 'reads_fully_aligned'
    ]
    assert list(df['label']) == ['chr1', 'chr2', 'Total']
    t = df[df['label'] == 'Total'].iloc[0]
    assert t['read_count'] == 2
    assert t['query_base'] == 300
    assert t['aligned_base'] == 280
    assert t['aligned_fraction_mean'] == pytest.approx(0.95)
    assert t['aligned_fraction_median'] == pytest.approx(0.95)
    assert t['mapq_mean'] == pytest.approx(45)
    assert t['reads_aligned_fraction_lt_0.9'] == 0
    assert t['reads_fully_aligned'] == 1


def test_mismatch_table_fixed_type_order_and_bin_total():
    a = Mismatch('chr1')
    a.add_mismatch('ct', 0.05)
    a.add_mismatch('ag', 0.05)
    total = Mismatch('Total')
    total.add_mismatch('ct', 0.05)
    total.add_mismatch('ag', 0.05)

    df = create_mismatch_normalized_read_location_table([a], total)
    # fixed canonical order, present types only: ag before ct
    assert list(df.columns) == ['label', 'bin', 'ag', 'ct', 'bin_total']
    row = df[df['label'] == 'chr1'].iloc[0]
    assert row['ag'] == 1
    assert row['ct'] == 1
    assert row['bin_total'] == 2