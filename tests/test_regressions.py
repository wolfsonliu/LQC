"""Regression tests for bugs found while smoke-testing lqc against a real BAM."""

from lqc import CS, Indel, Mismatch, Splice
from lqc.readstat import ReadStat
from lqc.report_figure import plot_splice_type_count
from lqc.report_html import html_add_splice_table
from lqc.report_table import (
    create_mismatch_normalized_read_location_table,
    create_splice_table,
)
from lqc.utils import list_bam_contigs


def test_cs_str_and_repr_use_mismatch_count_int():
    cs = CS.from_cs_tag_string(
        ':10*ag:5+tt:3~gt100ag:4-ccc:2',
        contig='chr1', start_pos=1000, strand='+',
    )
    assert str(cs) == 'CS tag: 1 introns, 1 mismatches, 2 indels'
    assert '1 mismatches' in repr(cs)
    assert '2 indels' in repr(cs)


def test_mismatch_type_count_unique_types():
    m = Mismatch('chr1')
    m.add_mismatch('ag', 0.1)
    m.add_mismatch('ag', 0.2)
    m.add_mismatch('ct', 0.1)
    assert m.get_type_count() == {'ag': 2, 'ct': 1}
    assert m.get_total_count() == 3
    assert m.get_locations() == [0.1, 0.2, 0.1]


def test_mismatch_table_with_no_mismatches_returns_empty_frame():
    m = Mismatch('chr1')
    total = Mismatch('Total')
    df = create_mismatch_normalized_read_location_table([m], total)
    assert list(df.columns) == ['label', 'bin']
    assert len(df) == 0


def test_indel_empty_gives_zero_median_and_no_crash():
    i = Indel('chr1')
    assert i.get_median_length() == 0
    assert i.get_mean_length() == 0
    assert i.get_longest_indel() == []
    assert i.get_most_abundant_indel() == []


def test_splice_html_table_with_no_introns_does_not_divide_by_zero():
    s = Splice('chr1')
    total = Splice('Total')
    table = create_splice_table([s], total)
    assert list(table.columns) == ['label']
    html = html_add_splice_table('<html>{%splice_table%}</html>', table, 0.0)
    assert '{%splice_table%}' not in html


def test_plot_splice_type_other_only_does_not_divide_by_zero():
    s = Splice('chr1')
    s.add_splice_pair_count_dict({'aa-aa': 5})
    fig = plot_splice_type_count([s], width=4, height=3)
    assert hasattr(fig, 'axes')


def test_list_bam_contigs(cs_bam):
    assert list_bam_contigs(cs_bam) == ['chr1']


def test_readstat_single_sum_keeps_identity_and_content():
    """sum([x]) == x (identity); relabeling safety lives at the call site."""
    rs = ReadStat('chr1')
    rs.add_read(100, insertion=1, deletion=0, mismatch=1, intron=0)
    total = sum([rs])
    assert total is rs
    assert total.get_read_count() == 1