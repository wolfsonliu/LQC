import pytest

from lqc._base import _LabelledStat
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