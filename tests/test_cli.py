import pytest

from lqc.cli import (
    ELEMENT_BAR_SPECS,
    MULTI_FIG_SPECS,
    READSTAT_BAR_SPECS,
    main,
)


def test_version_flag(capsys):
    with pytest.raises(SystemExit) as excinfo:
        main(['--version'])
    assert excinfo.value.code == 0
    assert capsys.readouterr().out.strip() == 'lqc 0.0.8'


def test_missing_bam_file_exits_with_error(capsys):
    with pytest.raises(SystemExit) as excinfo:
        main([])
    assert excinfo.value.code == 2
    assert 'bam-file' in capsys.readouterr().err


def _read_table_labels(path):
    with open(path, 'r') as fh:
        lines = [line.rstrip('\n') for line in fh]
    return [line.split('\t')[0] for line in lines[1:]]


def test_main_end_to_end_single_contig(cs_bam, tmp_path):
    out = tmp_path / 'out'
    # chr20 and chrY are not in the tiny fixture; they must be skipped, not crash.
    rc = main([
        '-b', cs_bam, '-o', str(out),
        '-c', 'chr1', 'chr20', 'chrY',
        '-t', '1',
        '--output-pickle', '--output-cs',
    ])
    assert rc == 0

    labels = _read_table_labels(out / 'table' / 'read_stat.txt')
    assert labels == ['chr1', 'Total']

    html = (out / 'LQC_report.html').read_text()
    assert '{%' not in html
    assert (out / 'result.pickle').exists()
    assert (out / 'read.cs').exists()


def test_main_errors_when_all_contigs_absent(cs_bam, tmp_path):
    with pytest.raises(ValueError):
        main(['-b', cs_bam, '-o', str(tmp_path / 'out'), '-c', 'chrZ', '-t', '1'])


def test_main_report_is_self_contained(cs_bam, tmp_path):
    out = tmp_path / 'out'
    rc = main(['-b', cs_bam, '-o', str(out), '-c', 'chr1', '-t', '1'])
    assert rc == 0
    html = (out / 'LQC_report.html').read_text()
    assert 'src="fig/' not in html
    assert 'cdn.jsdelivr.net' not in html
    assert 'img.shields.io' not in html
    assert 'id="lqc-data"' in html
    assert 'data:image/png;base64,' in html
    assert 'data:image/svg+xml;base64,' in html
    assert '{%' not in html


def test_tables_identical_across_thread_counts(cs_bam, tmp_path):
    out1 = tmp_path / 'o1'
    out2 = tmp_path / 'o2'
    assert main(['-b', cs_bam, '-o', str(out1), '-c', 'chr1', '-t', '1']) == 0
    assert main(['-b', cs_bam, '-o', str(out2), '-c', 'chr1', '-t', '2']) == 0
    for name in (
        'read_stat.txt', 'insertion.txt', 'deletion.txt',
        'mismatch.txt', 'splice.txt', 'mapping.txt', 'splice_all.txt',
    ):
        assert (out1 / 'table' / name).read_bytes() == \
               (out2 / 'table' / name).read_bytes()


def test_output_cs_identical_across_thread_counts(cs_bam, tmp_path):
    out1 = tmp_path / 'o1'
    out2 = tmp_path / 'o2'
    main(['-b', cs_bam, '-o', str(out1), '-c', 'chr1', '-t', '1', '--output-cs'])
    main(['-b', cs_bam, '-o', str(out2), '-c', 'chr1', '-t', '2', '--output-cs'])
    assert (out1 / 'read.cs').read_bytes() == (out2 / 'read.cs').read_bytes()
    assert not list(out1.glob('.readcs-*.tmp'))
    assert not list(out2.glob('.readcs-*.tmp'))


def test_main_outputs_expected_figure_files(cs_bam, tmp_path):
    out = tmp_path / 'out'
    rc = main(['-b', cs_bam, '-o', str(out), '-t', '1',
               '-c', 'chr1', '--log-level', 'DEBUG'])
    assert rc == 0
    produced = {p.name for p in (out / 'fig').glob('*')}
    for expected in ['readstat_bar_Read_count.png',
                     'readstat_bar_mean_element_per_read.Total.png',
                     'insertion_bar_count.png',
                     'splice_type.png']:
        assert expected in produced


def test_figure_spec_lists_match_contract():
    assert READSTAT_BAR_SPECS[0]._fields == ('feature', 'stem')
    assert MULTI_FIG_SPECS[0]._fields == ('stem', 'plot_func', 'data_key')
    assert ELEMENT_BAR_SPECS[0]._fields == ('stem', 'data_key', 'kind')

    assert [(s.feature, s.stem) for s in READSTAT_BAR_SPECS] == [
        ('Read count', 'readstat_bar_Read_count'),
        ('Median read length', 'readstat_bar_Median_read_length'),
        ('Mean read length', 'readstat_bar_Mean_read_length'),
        ('Insertions per read', 'readstat_bar_insertions_per_read'),
        ('Insertions per read per kb',
         'readstat_bar_insertions_per_read_per_kb'),
        ('Deletions per read', 'readstat_bar_deletions_per_read'),
        ('Deletions per read per kb', 'readstat_bar_deletions_per_read_per_kb'),
        ('Mismatches per read', 'readstat_bar_mismatches_per_read'),
        ('Mismatches per read per kb', 'readstat_bar_mismatches_per_read_per_kb'),
        ('Mean intron number', 'readstat_bar_Mean_intron_number'),
        ('N50', 'readstat_bar_N50'),
        ('L50', 'readstat_bar_L50'),
    ]
    assert [(s.stem, s.data_key, s.kind) for s in ELEMENT_BAR_SPECS] == [
        ('insertion_bar_count', 'insertion', 'Insertion'),
        ('deletion_bar_count', 'deletion', 'Deletion'),
        ('mismatch_bar_count', 'mismatch', 'Mismatch'),
        ('intron_bar_count', 'splice', 'Intron'),
    ]
    assert [(s.stem, s.plot_func.__name__, s.data_key)
            for s in MULTI_FIG_SPECS] == [
        ('readstat_bar_mean_element_per_read',
         'plot_readstat_bar_mean_element_per_read', 'readstat'),
        ('readstat_bar_mean_element_per_read_per_kb',
         'plot_readstat_bar_mean_element_per_read_per_kb', 'readstat'),
        ('readstat_line_cumulative_length',
         'plot_readstat_cumulative_length', 'readstat'),
        ('readstat_bar_ratio_with_element',
         'plot_readstat_bar_ratio_with_element', 'readstat'),
        ('readstat_hist_length', 'plot_readstat_length_hist', 'readstat'),
        ('insertion_hist_length', 'plot_indel_hist_length', 'insertion'),
        ('deletion_hist_length', 'plot_indel_hist_length', 'deletion'),
        ('insertion_hist_location', 'plot_indel_hist_location', 'insertion'),
        ('deletion_hist_location', 'plot_indel_hist_location', 'deletion'),
        ('mismatch_type', 'plot_mismatch_type_count', 'mismatch'),
        ('mismatch_hist_location', 'plot_mismatch_hist_location', 'mismatch'),
        ('splice_type', 'plot_splice_type_count', 'splice'),
        ('mapping_hist_mapq', 'plot_mapping_mapq_hist', 'readstat'),
        ('mapping_hist_aligned_fraction',
         'plot_mapping_aligned_fraction_hist', 'readstat'),
        ('mapping_scatter_aligned_vs_query',
         'plot_mapping_aligned_vs_query', 'readstat'),
    ]
    data_keys = ({s.data_key for s in MULTI_FIG_SPECS}
                 | {s.data_key for s in ELEMENT_BAR_SPECS})
    assert data_keys <= {'readstat', 'insertion', 'deletion',
                         'mismatch', 'splice'}
