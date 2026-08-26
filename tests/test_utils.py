import pytest

from lqc import utils


def test_convert_complement():
    assert utils.convert_complement('acgtn-') == 'tgcan-'


def test_convert_reverse_complement():
    assert utils.convert_reverse_complement('aacg') == 'cgtt'
    assert utils.convert_reverse_complement('gt-ag') == 'ct-ac'


def test_convert_reverse_complement_uppercase_unchanged():
    # uppercase is not in _COMPLEMENT_TABLE, so it passes through unchanged
    assert utils.convert_reverse_complement('AG') == 'GA'


def test_check_cs_md_tag():
    tags = [('cs', ':5'), ('MD', '5'), ('NM', 0)]
    assert utils.check_cs_md_tag(tags) == ['cs', 'MD']


def test_bam_or_sam():
    assert utils.bam_or_sam('x.bam') == 'BAM'
    assert utils.bam_or_sam('x.sam') == 'SAM'


def test_bam_or_sam_rejects_unknown_extension():
    with pytest.raises(AssertionError):
        utils.bam_or_sam('x.txt')