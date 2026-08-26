from lqc.stat import (
    plan_tasks,
    prefetch_records,
    reduce_blocks_to_contigs,
    stat_element_from_bam_by_contig,
    stat_records,
    stat_region,
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


def test_stat_records_captures_mapping_fields(cs_bam):
    contig, records = prefetch_records(cs_bam, ['chr1'], 'cs')[0]
    read0 = records[0]
    assert read0.mapping_quality == 60
    assert read0.aligned_length == 10
    block = stat_records((0, contig, records), None, 'cs')
    assert block.readstat.get_mapping_qualities() == [60, 60]
    assert block.readstat.get_aligned_lengths() == [10, 10]


def test_record_from_read_md_captures_mapping_fields():
    from lqc.stat import record_from_read

    class FakeRead:
        is_reverse = False
        reference_start = 5
        query_sequence = 'ACGTACGTAC'
        query_length = 10
        mapping_quality = 42
        query_alignment_length = 10
        query_name = 'read1'
        cigarstring = '10M'
        reference_end = 15

        def get_tag(self, tag):
            assert tag == 'MD'
            return '10'

    record = record_from_read(FakeRead(), 'chr1', 'md')
    assert record.mapping_quality == 42
    assert record.aligned_length == 10
    assert record.query_length == 10


def test_record_from_read_cs_uses_query_length_not_sequence():
    from lqc.stat import record_from_read

    class FakeRead:
        is_reverse = False
        reference_start = 0
        mapping_quality = 60
        query_alignment_length = 10
        query_name = 'read1'
        query_length = 10

        @property
        def query_sequence(self):
            raise AssertionError(
                'cs path must use query_length, not materialize query_sequence'
            )

        def get_tag(self, tag):
            assert tag == 'cs'
            return ':10'

    record = record_from_read(FakeRead(), 'chr1', 'cs')
    assert record.query_length == 10


def test_plan_tasks_splits_by_target_reads(cs_bam):
    # cs_bam: chr1 (LN=1,000,000) with 2 mapped reads.
    # target_reads=1  -> ceil(2/1)=2 windows of 500,000 bp each.
    assert plan_tasks(cs_bam, ['chr1'], target_reads=1) == [
        (0, 'chr1', 0, 500000),
        (1, 'chr1', 500000, 1000000),
    ]
    # target_reads >= mapped -> a single full-contig window.
    assert plan_tasks(cs_bam, ['chr1'], target_reads=100) == [
        (0, 'chr1', 0, 1000000),
    ]
    # A contig with zero mapped reads is skipped.
    assert plan_tasks(cs_bam, ['chr1', 'chr999'], target_reads=1) == [
        (0, 'chr1', 0, 500000),
        (1, 'chr1', 500000, 1000000),
    ]


def test_stat_region_matches_serial(cs_bam):
    # Two windows over chr1 (both reads fall in the first window; the second
    # is empty). The reduced result must equal the serial single-contig path.
    tasks = plan_tasks(cs_bam, ['chr1'], target_reads=1)
    blocks = [
        stat_region(task, cs_bam, None, 'cs')
        for task in tasks
    ]
    reduced = reduce_blocks_to_contigs(blocks, ['chr1'])
    serial = stat_element_from_bam_by_contig(cs_bam, None, 'chr1', 'cs')
    _assert_stat_tuples_equal(reduced[0], serial)
    assert reduced[0][0].get_read_count() == 2


def test_stat_region_assigns_boundary_read_once(cs_bam_boundary):
    tasks = plan_tasks(cs_bam_boundary, ['chr1'], target_reads=1)
    blocks = [
        stat_region(task, cs_bam_boundary, None, 'cs')
        for task in tasks
    ]
    # 'span' (reference_start=499990, 20M) crosses the 500000 boundary; it must
    # be counted in the first window and skipped in the second.
    assert [b.readstat.get_read_count() for b in blocks] == [1, 1]
    reduced = reduce_blocks_to_contigs(blocks, ['chr1'])
    serial = stat_element_from_bam_by_contig(cs_bam_boundary, None, 'chr1', 'cs')
    _assert_stat_tuples_equal(reduced[0], serial)
