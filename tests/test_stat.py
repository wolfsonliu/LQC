from lqc.stat import (
    prefetch_records,
    reduce_blocks_to_contigs,
    stat_element_from_bam_by_contig,
    stat_records,
)
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


def _assert_stat_tuples_equal(a, b):
    rs1, ins1, dele1, mis1, spl1 = a
    rs2, ins2, dele2, mis2, spl2 = b
    assert rs1.get_read_count() == rs2.get_read_count()
    assert rs1.label == rs2.label
    assert ins1.label == ins2.label
    assert dele1.label == dele2.label
    assert mis1.label == mis2.label
    assert spl1.label == spl2.label
    assert rs1.get_total_base() == rs2.get_total_base()
    assert rs1.get_lengths() == rs2.get_lengths()
    assert rs1.get_insertions() == rs2.get_insertions()
    assert rs1.get_deletions() == rs2.get_deletions()
    assert rs1.get_mismatches() == rs2.get_mismatches()
    assert rs1.get_introns() == rs2.get_introns()
    assert ins1.get_indel_count() == ins2.get_indel_count()
    assert dele1.get_indel_count() == dele2.get_indel_count()
    assert mis1.get_type_count() == mis2.get_type_count()
    assert spl1.get_splice_pair_count_dict() == spl2.get_splice_pair_count_dict()


def test_prefetch_and_reduce_matches_serial(cs_bam, cs_bam_with_indel):
    for bam in (cs_bam, cs_bam_with_indel):
        serial = stat_element_from_bam_by_contig(
            bam_file=bam, genome_file=None, contig='chr1', method='cs'
        )
        contig, records = prefetch_records(bam, ['chr1'], 'cs')[0]
        half = (len(records) + 1) // 2
        chunks = [records[:half], records[half:]]
        blocks = [
            stat_records((i, contig, chunk), None, 'cs')
            for i, chunk in enumerate(chunks)
            if chunk
        ]
        reduced = reduce_blocks_to_contigs(blocks, ['chr1'])
        _assert_stat_tuples_equal(reduced[0], serial)


def test_stat_records_writes_cs_lines(cs_bam, tmp_path):
    contig, records = prefetch_records(cs_bam, ['chr1'], 'cs')[0]
    block = stat_records((0, contig, records), None, 'cs', cs_dir=str(tmp_path))
    assert block.cs_path == str(tmp_path / '.readcs-00000000.tmp')
    assert (tmp_path / '.readcs-00000000.tmp').read_text().splitlines() == [
        'read1\tchr1\t0\t5\t:\t5',
        'read1\tchr1\t5\t6\t*\tag',
        'read1\tchr1\t6\t10\t:\t4',
        'read2\tchr1\t100\t110\t:\t10',
    ]