from lqc.stat import stat_element_from_bam_by_contig
from lqc.utils import check_bam_with_cs_or_md


def test_check_bam_with_cs_or_md(cs_bam):
    assert check_bam_with_cs_or_md(cs_bam) == 'cs'


def test_stat_element_from_bam_by_contig(cs_bam):
    (
        readstat,
        insertion,
        deletion,
        mismatch,
        splice,
    ) = stat_element_from_bam_by_contig(
        bam_file=cs_bam,
        genome_file=None,
        contig='chr1',
        method='cs',
    )

    assert readstat.get_read_count() == 2
    assert readstat.get_total_base() == 20
    assert readstat.get_lengths() == [10, 10]
    assert readstat.get_mismatches() == [1, 0]
    assert readstat.get_introns() == [0, 0]
    assert readstat.get_insertions() == [0, 0]
    assert readstat.get_deletions() == [0, 0]

    assert mismatch.get_total_count() == 1
    assert mismatch.get_type_count() == {'ag': 1}

    assert insertion.get_total_count() == 0
    assert deletion.get_total_count() == 0
    assert splice.get_total_splice_pair_count() == 0


def test_stat_element_from_bam_records_indels(cs_bam_with_indel):
    (
        readstat,
        insertion,
        deletion,
        mismatch,
        splice,
    ) = stat_element_from_bam_by_contig(
        bam_file=cs_bam_with_indel,
        genome_file=None,
        contig='chr1',
        method='cs',
    )

    assert readstat.get_read_count() == 1
    assert readstat.get_insertions() == [1]
    assert readstat.get_deletions() == [1]

    assert insertion.get_total_count() == 1
    assert insertion.get_indel_count() == {'tt': 1}

    assert deletion.get_total_count() == 1
    assert deletion.get_indel_count() == {'cc': 1}
    assert deletion.get_total_length() == 2

    assert mismatch.get_total_count() == 0
    assert splice.get_total_splice_pair_count() == 0