import pytest

from lqc import Splice


def test_label_must_be_str():
    with pytest.raises(TypeError):
        Splice([])


def test_add_and_count():
    s = Splice('sp')
    s.add_splice_pair('gt-ag')
    s.add_splice_pair('gt-ag')
    s.add_splice_pair('gc-ag')
    assert s.get_splice_pair_count_dict() == {'gt-ag': 2, 'gc-ag': 1}
    assert s.get_total_splice_pair_count() == 3


def test_most_abundant_pair():
    s = Splice('sp')
    s.add_splice_pair_list(['gt-ag', 'gt-ag', 'gc-ag'])
    assert s.get_most_abundant_splice_pair() == [('gt-ag', 2)]


def test_reverse_complement():
    s = Splice('sp')
    s.add_splice_pair('gt-ag')
    rc = s.convert_reverse_complement()
    assert rc.get_splice_pair_count_dict() == {'ct-ac': 1}


def test_add_two_splices():
    a = Splice('x')
    a.add_splice_pair('gt-ag')
    b = Splice('y')
    b.add_splice_pair('gt-ag')
    b.add_splice_pair('gc-ag')
    c = a + b
    assert c.get_splice_pair_count_dict() == {'gt-ag': 2, 'gc-ag': 1}


def test_add_counts_accumulate_without_aliasing():
    a = Splice('a')
    b = Splice('b')
    a.add_splice_pair('gt-ag')
    b.add_splice_pair('ct-ac')
    b.add_splice_pair('gt-ag')
    total = a + b
    assert total.get_splice_pair_count_dict() == {'gt-ag': 2, 'ct-ac': 1}
    # the operands are untouched (no aliasing)
    assert a.get_splice_pair_count_dict() == {'gt-ag': 1}
    assert b.get_splice_pair_count_dict() == {'ct-ac': 1, 'gt-ag': 1}