import pytest

from lqc import Indel


@pytest.fixture()
def indel():
    i = Indel('ins')
    i.add_indel('a', 0.1)
    i.add_indel('a', 0.2)
    i.add_indel('tt', 0.5)
    return i


def test_label_must_be_str():
    with pytest.raises(TypeError):
        Indel(0)


def test_counts_and_lengths(indel):
    assert indel.get_total_count() == 3
    assert indel.get_total_length() == 4
    assert indel.get_lengths() == [1, 1, 2]
    assert indel.get_indel_count() == {'a': 2, 'tt': 1}


def test_mean_and_median_length(indel):
    assert indel.get_mean_length() == pytest.approx(4 / 3)
    assert indel.get_median_length() == pytest.approx(1.0)


def test_longest_and_most_abundant(indel):
    assert indel.get_longest_indel() == ['tt']
    assert indel.get_most_abundant_indel() == ['a']


def test_location_bin_count(indel):
    bins = indel.get_location_bin_count()
    assert bins['[0.0,0.25)'] == 2
    assert bins['[0.5,0.75)'] == 1
    assert sum(bins.values()) == 3


def test_reverse_complement(indel):
    rc = indel.convert_reverse_complement()
    assert rc.get_indel_count() == {'t': 2, 'aa': 1}


def test_add_two_indels():
    a = Indel('x')
    a.add_indel('a', 0.1)
    b = Indel('y')
    b.add_indel('a', 0.2)
    b.add_indel('tt', 0.5)
    c = a + b
    assert c.get_indel_count() == {'a': 2, 'tt': 1}
    assert c.get_total_count() == 3


def test_indel_empty_plus_nonempty_merge():
    empty = Indel('e')
    full = Indel('f')
    full.add_indel('a', 0.1)
    full.add_indel('b', 0.5)
    total = sum([empty, full])
    assert total.get_total_count() == 2
    assert total.get_lengths() == [1, 1]
    assert total.get_locations() == [0.1, 0.5]
    assert total.get_indel_count() == {'a': 1, 'b': 1}


def test_indel_locations_order_preserved():
    a = Indel('a')
    b = Indel('b')
    a.add_indel('x', 0.25)
    a.add_indel('y', 0.9)
    b.add_indel('z', 0.33)
    total = a + b
    assert total.get_locations() == [0.25, 0.9, 0.33]