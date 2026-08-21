import pytest

from lqc.cli import main


def test_version_flag(capsys):
    with pytest.raises(SystemExit) as excinfo:
        main(['--version'])
    assert excinfo.value.code == 0
    assert capsys.readouterr().out.strip() == 'lqc 0.0.5'


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