import pytest

from lqc.readstat import ReadStat


@pytest.fixture()
def rs():
    r = ReadStat('sample')
    r.add_read(200, insertion=2, deletion=4, mismatch=6, intron=8)
    r.add_read(200, insertion=4, deletion=8, mismatch=12, intron=16)
    return r


def test_label_must_be_str():
    with pytest.raises(TypeError):
        ReadStat(3)


def test_read_count_and_total_base(rs):
    assert rs.get_read_count() == 2
    assert rs.get_total_base() == 400
    assert len(rs) == 2


def test_length_stats(rs):
    assert rs.get_lengths() == [200, 200]
    assert rs.get_min_length() == 200
    assert rs.get_max_length() == 200
    assert rs.get_mean_length() == 200
    assert rs.get_median_length() == 200


def test_n50_l50(rs):
    assert rs.get_length_NL(50) == (200, 1)


def test_element_means(rs):
    assert rs.get_mean_insertions() == 3
    assert rs.get_mean_deletions() == 6
    assert rs.get_mean_mismatches() == 9
    assert rs.get_mean_introns() == 12


def test_per_base_rates(rs):
    assert rs.insertions_per_base() == pytest.approx(6 / 400)
    assert rs.deletions_per_base() == pytest.approx(12 / 400)
    assert rs.mismatches_per_base() == pytest.approx(18 / 400)


def test_length_normalized_means(rs):
    assert rs.get_mean_length_normalized_insertions() == pytest.approx(0.015)
    assert rs.get_mean_length_normalized_deletions() == pytest.approx(0.03)
    assert rs.get_mean_length_normalized_mismatches() == pytest.approx(0.045)
    assert rs.get_mean_length_normalized_introns() == pytest.approx(0.06)


def test_read_count_with_n_elements(rs):
    assert rs.get_read_count_with_n_insertions(2) == 1
    assert rs.get_read_count_with_n_deletions(4) == 1
    assert rs.get_read_count_with_n_mismatches(0) == 0


def test_min_max_medians(rs):
    assert rs.get_min_insertions() == 2
    assert rs.get_max_insertions() == 4
    assert rs.get_median_insertion() == 3
    assert rs.get_median_deletion() == 6
    assert rs.get_median_mismatch() == 9
    assert rs.get_median_intron() == 12


def test_add_two_readstats():
    a = ReadStat('a')
    a.add_read(100, insertion=1, deletion=0, mismatch=1, intron=0)
    b = ReadStat('b')
    b.add_read(200, insertion=0, deletion=1, mismatch=0, intron=1)
    c = a + b
    assert c.label == 'a b'
    assert c.get_read_count() == 2
    assert c.get_total_base() == 300
    assert c.get_lengths() == [100, 200]


def test_sum_readstats():
    a = ReadStat('a')
    a.add_read(100, insertion=1, deletion=0, mismatch=1, intron=0)
    b = ReadStat('b')
    b.add_read(200, insertion=0, deletion=1, mismatch=0, intron=1)
    total = sum([a, b])
    assert total.get_read_count() == 2
    assert total.get_total_base() == 300


def test_mapping_fields_default_to_zero(rs):
    assert rs.get_mapping_qualities() == [0, 0]
    assert rs.get_aligned_lengths() == [0, 0]
    assert rs.get_aligned_fractions() == [0.0, 0.0]
    assert rs.get_total_aligned_base() == 0


def test_mapping_getters():
    r = ReadStat('chr1')
    r.add_read(200, insertion=2, deletion=4, mismatch=6, intron=8,
               mapping_quality=60, aligned_length=180)
    r.add_read(100, insertion=0, deletion=0, mismatch=0, intron=0,
               mapping_quality=30, aligned_length=100)
    assert r.get_mapping_qualities() == [60, 30]
    assert r.get_aligned_lengths() == [180, 100]
    assert r.get_aligned_fractions() == [0.9, 1.0]
    assert r.get_total_aligned_base() == 280
    assert r.get_mean_mapping_quality() == 45
    assert r.get_median_mapping_quality() == 45
    assert r.get_mean_aligned_fraction() == pytest.approx(0.95)
    assert r.get_median_aligned_fraction() == pytest.approx(0.95)
    assert r.get_read_count_with_aligned_fraction_below(0.9) == 0
    assert r.get_read_count_with_aligned_fraction_below(0.95) == 1
    assert r.get_read_count_fully_aligned() == 1


def test_aligned_base_rates():
    r = ReadStat('chr1')
    r.add_read(200, insertion=2, deletion=4, mismatch=6, intron=0,
               mapping_quality=60, aligned_length=100)
    r.add_read(200, insertion=2, deletion=4, mismatch=6, intron=0,
               mapping_quality=60, aligned_length=100)
    assert r.get_total_aligned_base() == 200
    assert r.insertions_per_aligned_base() == pytest.approx(4 / 200)
    assert r.deletions_per_aligned_base() == pytest.approx(8 / 200)
    assert r.mismatches_per_aligned_base() == pytest.approx(12 / 200)


def test_aligned_base_rate_zero_when_no_aligned():
    r = ReadStat('chr1')
    r.add_read(200, insertion=2, deletion=4, mismatch=6, intron=0)
    assert r.insertions_per_aligned_base() == 0
    assert r.deletions_per_aligned_base() == 0
    assert r.mismatches_per_aligned_base() == 0


def test_add_two_readstats_sums_aligned_base():
    a = ReadStat('a')
    a.add_read(100, insertion=1, deletion=0, mismatch=1, intron=0,
               mapping_quality=60, aligned_length=90)
    b = ReadStat('b')
    b.add_read(200, insertion=0, deletion=1, mismatch=0, intron=1,
               mapping_quality=30, aligned_length=200)
    c = a + b
    assert c.get_total_aligned_base() == 290