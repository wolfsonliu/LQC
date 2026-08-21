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