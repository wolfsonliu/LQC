import pytest

from lqc import Indel, Splice
from lqc.readstat import ReadStat
from lqc.report_table import (
    create_indel_summary_table,
    create_mapping_table,
    create_readstat_table,
    create_splice_table,
)

READSTAT_COLUMNS = [
    'label', 'read_count', 'total_base',
    'read_length_mean', 'read_length_median',
    'read_length_N50', 'read_length_L50',
    'mean_insertion_per_read', 'mean_insertion_per_read_per_kb',
    'insertion_per_kb',
    'mean_deletion_per_read', 'mean_deletion_per_read_per_kb',
    'deletion_per_kb',
    'mean_mismatch_per_read', 'mean_mismatch_per_read_per_kb',
    'mismatch_per_kb',
    'mean_intron_per_read', 'mean_intron_per_read_per_kb',
]


def _two_readstat():
    r1 = ReadStat('chr1')
    r1.add_read(200, insertion=2, deletion=4, mismatch=6, intron=8)
    r2 = ReadStat('chr2')
    r2.add_read(200, insertion=4, deletion=8, mismatch=12, intron=16)
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
    assert list(df.columns) == ['label', 'gt-ag']
    assert df[df['label'] == 'Total']['gt-ag'].values[0] == 2


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