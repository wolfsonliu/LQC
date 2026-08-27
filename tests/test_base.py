import pytest

from lqc._base import _LabelledStat, concat_stats
from lqc.mismatch import Mismatch
from lqc.readstat import ReadStat


def test_labelled_stat_rejects_nonstring_label():
    with pytest.raises(TypeError):
        _LabelledStat(123)


def test_labelled_stat_radd_zero_returns_self():
    stat = _LabelledStat('x')
    assert (0 + stat) is stat


def test_readstat_still_rejects_nonstring_label():
    with pytest.raises(TypeError):
        ReadStat(123)


def test_readstat_add_wrong_type_raises_typeerror():
    with pytest.raises(TypeError):
        ReadStat('a') + Mismatch('b')


def test_concat_stats_sums_and_relabels():
    a = ReadStat('chr1')
    a.add_read(100, insertion=1, deletion=0, mismatch=1, intron=0)
    b = ReadStat('chr2')
    b.add_read(200, insertion=0, deletion=1, mismatch=0, intron=1)
    total = concat_stats([a, b], label='Total')
    assert total.get_read_count() == 2
    assert total.label == 'Total'
    assert a.label == 'chr1'
    assert b.label == 'chr2'


def test_concat_stats_single_element_does_not_alias():
    a = ReadStat('chr1')
    a.add_read(100, insertion=1, deletion=0, mismatch=1, intron=0)
    total = concat_stats([a], label='Total')
    assert total is not a
    assert total.label == 'Total'
    assert a.label == 'chr1'
    assert total.get_read_count() == 1


def test_concat_stats_requires_nonempty():
    with pytest.raises(ValueError):
        concat_stats([], label='Total')