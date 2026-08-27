from lqc import get_html_template
from lqc.readstat import ReadStat
from lqc.report_html import (
    _replace_tokens,
    html_add_deletion_table,
    html_add_insertion_table,
    html_add_mapping,
    html_add_readstat_table,
)
from lqc.report_table import create_mapping_table, create_readstat_table


def _readstat_table_with_total():
    r1 = ReadStat('chr1')
    r1.add_read(200, insertion=2, deletion=4, mismatch=6, intron=8)
    r2 = ReadStat('chr2')
    r2.add_read(200, insertion=4, deletion=8, mismatch=12, intron=16)
    total = sum([r1, r2])
    total.label = 'Total'
    return create_readstat_table([r1, r2], total)


def test_html_add_readstat_table_substitutes_placeholders():
    table = _readstat_table_with_total()
    html = html_add_readstat_table(get_html_template(), table)
    assert '{%readstat_total_read_count%}' not in html
    assert '{%readstat_total_base%}' not in html
    assert '{%readstat_N50%}' not in html
    assert '{%readstat_L50%}' not in html
    assert '{%readstat_table%}' not in html
    assert '<tr class="table-secondary">' in html


def test_html_add_readstat_table_contains_total_values():
    table = _readstat_table_with_total()
    html = html_add_readstat_table(get_html_template(), table)
    # Total read count 2, total base 400 (2 reads x 200 bp)
    assert '400' in html


def test_html_add_mapping_substitutes_placeholders():
    r1 = ReadStat('chr1')
    r1.add_read(200, insertion=2, deletion=4, mismatch=6, intron=8,
                mapping_quality=60, aligned_length=180)
    r2 = ReadStat('chr2')
    r2.add_read(100, insertion=0, deletion=0, mismatch=0, intron=0,
                mapping_quality=30, aligned_length=100)
    total = sum([r1, r2])
    total.label = 'Total'
    table = create_mapping_table([r1, r2], total)

    html = html_add_mapping(get_html_template(), table)
    assert '{%mapping_aligned_fraction_mean%}' not in html
    assert '{%mapping_aligned_fraction_median%}' not in html
    assert '{%mapping_mapq_mean%}' not in html
    assert '{%mapping_mapq_median%}' not in html
    assert '{%mapping_table%}' not in html
    assert '<tr class="table-secondary">' in html


def test_replace_tokens_substitutes_all_in_order():
    html = _replace_tokens(
        '<p>{%a%}-{%b%}</p>', {'{%a%}': '1', '{%b%}': '2'}
    )
    assert html == '<p>1-2</p>'


def test_indel_tables_share_implementation():
    import pandas as pd
    table = pd.DataFrame([{
        'label': 'Total', 'total_count': 1, 'total_length': 2,
        'mean_length': 2.0, 'median_length': 2.0,
    }])
    ins = html_add_insertion_table('<h>{%insertion_table%}</h>', table, 0.0)
    dele = html_add_deletion_table('<h>{%deletion_table%}</h>', table, 0.0)
    assert '{%insertion_table%}' not in ins
    assert '{%deletion_table%}' not in dele