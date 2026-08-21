import pytest

from lqc import Mismatch


@pytest.fixture()
def mismatch():
    m = Mismatch('mis')
    m.add_mismatch('ag', 0.1)
    m.add_mismatch('ag', 0.2)
    m.add_mismatch('ct', 0.1)
    return m


def test_label_must_be_str():
    with pytest.raises(TypeError):
        Mismatch(None)


def test_counts(mismatch):
    assert mismatch.get_total_count() == 3
    assert mismatch.get_type_count() == {'ag': 2, 'ct': 1}
    assert mismatch.get_locations() == [0.1, 0.2, 0.1]


def test_location_bin_count(mismatch):
    bins = mismatch.get_location_bin_count()
    assert bins['[0.0,0.25)'] == 3
    assert sum(bins.values()) == 3


def test_location_bin_count_by_type(mismatch):
    by_type = mismatch.get_location_bin_count_by_type()
    assert by_type['ag']['[0.0,0.25)'] == 2
    assert by_type['ct']['[0.0,0.25)'] == 1


def test_complement(mismatch):
    comp = mismatch.convert_complement()
    assert comp.get_type_count() == {'tc': 2, 'ga': 1}


def test_add_two_mismatches():
    a = Mismatch('x')
    a.add_mismatch('ag', 0.1)
    b = Mismatch('y')
    b.add_mismatch('ag', 0.2)
    b.add_mismatch('ct', 0.1)
    c = a + b
    assert c.get_type_count() == {'ag': 2, 'ct': 1}
    assert c.get_total_count() == 3