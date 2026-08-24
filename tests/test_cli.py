import pytest

from lqc.cli import main


def test_version_flag(capsys):
    with pytest.raises(SystemExit) as excinfo:
        main(['--version'])
    assert excinfo.value.code == 0
    assert capsys.readouterr().out.strip() == 'lqc 0.0.6'


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