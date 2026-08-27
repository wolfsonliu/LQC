import os
from collections import defaultdict
from contextlib import contextmanager
from typing import NamedTuple, Optional

import pysam

from lqc.cs import CS
from lqc.indel import Indel
from lqc.mismatch import Mismatch
from lqc.readstat import ReadStat
from lqc.splice import Splice
from lqc.utils import bam_or_sam, convert_complement, convert_reverse_complement


class ReadRecord(NamedTuple):
    """Lightweight per-read fields needed to build a CS and accumulate stats.

    ``method`` is ``'cs'`` or ``'md'``. ``cs_string`` is set for the cs path;
    the remaining optional fields are set for the MD+CIGAR path.
    """
    method: str
    contig: str
    start_pos: int
    strand: str
    query_length: int
    mapping_quality: int
    aligned_length: int
    query_name: str
    cs_string: Optional[str] = None
    cigar: Optional[str] = None
    md_string: Optional[str] = None
    query_sequence: Optional[str] = None
    reference_end: Optional[int] = None


class StatBlock(NamedTuple):
    """Result of one parallel stat chunk (see ``stat_records``)."""
    contig: str
    readstat: ReadStat
    insertion: Indel
    deletion: Indel
    mismatch: Mismatch
    splice: Splice
    cs_path: Optional[str] = None


class ContigStats(NamedTuple):
    """The five per-contig accumulators, in canonical order."""
    readstat: ReadStat
    insertion: Indel
    deletion: Indel
    mismatch: Mismatch
    splice: Splice


def record_from_read(read, contig, method):
    """Extract the minimal fields from a pysam read for stat accumulation."""
    strand = '-' if read.is_reverse else '+'
    if method == 'cs':
        return ReadRecord(
            method='cs',
            contig=contig,
            start_pos=read.reference_start,
            strand=strand,
            query_length=read.query_length,
            mapping_quality=read.mapping_quality,
            aligned_length=read.query_alignment_length,
            query_name=read.query_name,
            cs_string=read.get_tag('cs'),
        )
    return ReadRecord(
        method='md',
        contig=contig,
        start_pos=read.reference_start,
        strand=strand,
        query_length=read.query_length,
        mapping_quality=read.mapping_quality,
        aligned_length=read.query_alignment_length,
        query_name=read.query_name,
        cigar=read.cigarstring,
        md_string=read.get_tag('MD'),
        query_sequence=read.query_sequence,
        reference_end=read.reference_end,
    )


def _build_cs(record, genome):
    if record.method == 'cs':
        return CS.from_cs_tag_string(
            record.cs_string,
            contig=record.contig,
            start_pos=record.start_pos,
            strand=record.strand,
        )
    ref_seq = genome.fetch(
        record.contig, record.start_pos, record.reference_end
    )
    return CS.from_cigar_string(
        record.cigar,
        record.md_string,
        record.query_sequence,
        ref_seq,
        contig=record.contig,
        start_pos=record.start_pos,
        strand=record.strand,
    )


def process_record(record, genome,
                   readstat, insertion, deletion, mismatch, splice):
    """Accumulate one read into the five stat objects and return its CS."""
    cs = _build_cs(record, genome)
    strand = record.strand
    readstat.add_read(
        length=record.query_length,
        insertion=cs.get_insertion_count(),
        deletion=cs.get_deletion_count(),
        mismatch=cs.get_mismatch_count(),
        intron=cs.get_intron_count(),
        mapping_quality=record.mapping_quality,
        aligned_length=record.aligned_length,
    )
    insertions, deletions, mismatches = cs.get_indel_mismatches(
        coordinate='normalized_read'
    )
    for a, _, _, d in insertions:
        insertion.add_indel(
            d if strand == '+' else convert_reverse_complement(d), a
        )
    for a, _, _, d in deletions:
        deletion.add_indel(
            d if strand == '+' else convert_reverse_complement(d), a
        )
    for a, _, _, d in mismatches:
        mismatch.add_mismatch(
            d if strand == '+' else convert_complement(d), a
        )
    for s_str, s_ct in cs.get_splice_pair_count().items():
        splice.add_splice_pair_count_dict({
            s_str if strand == '+' else convert_reverse_complement(s_str): s_ct
        })
    return cs


def prefetch_records(bam_file, contigs, method):
    """Return ``[(contig, [ReadRecord, ...]), ...]`` in BAM traversal order.

    Materializes every read's ``ReadRecord`` in memory before pooling. For the
    cs path this is small (a cs string plus a few ints); for the MD path each
    record additionally carries the full ``query_sequence``.
    """
    file_type = bam_or_sam(bam_file)
    file_read = "rb" if file_type == "BAM" else "r"
    bam = pysam.AlignmentFile(bam_file, file_read)
    per_contig = []
    for contig in contigs:
        records = [
            record_from_read(read, contig, method)
            for read in bam.fetch(contig)
        ]
        per_contig.append((contig, records))
    bam.close()
    return per_contig


def plan_tasks(bam_file, contigs, target_reads=10000):
    """Return ``[(index, contig, start, end), ...]`` coordinate windows.

    Splits each contig's coordinate span into windows of roughly
    ``target_reads`` mapped reads so a large contig can be processed by
    several workers. Windows are emitted in requested-contig order, then
    ascending coordinate order, so concatenating per-window results preserves
    the current single-contig traversal order.
    """
    file_type = bam_or_sam(bam_file)
    file_read = "rb" if file_type == "BAM" else "r"
    bam = pysam.AlignmentFile(bam_file, file_read)
    try:
        stats = {
            stat.contig: stat.mapped
            for stat in bam.get_index_statistics()
        }
        tasks = []
        index = 0
        for contig in contigs:
            mapped = stats.get(contig, 0)
            if mapped == 0:
                continue
            length = bam.get_reference_length(contig)
            n_windows = max(1, (mapped + target_reads - 1) // target_reads)
            for w in range(n_windows):
                start = (w * length) // n_windows
                end = ((w + 1) * length) // n_windows
                tasks.append((index, contig, start, end))
                index += 1
    finally:
        bam.close()
    return tasks


@contextmanager
def _open_optional_cs(path):
    """Yield an open write handle for ``path``, or ``None`` if path is ``None``."""
    if path is None:
        yield None
    else:
        with open(path, 'w') as fh:
            yield fh


def stat_records(task, genome_file, method, cs_dir=None):
    """Process one ``(index, contig, records)`` task; return per-chunk stats.

    When ``cs_dir`` is set, also write this chunk's read.cs lines to a
    per-chunk temp file and return its path for ordered merging.
    """
    index, contig, records = task
    readstat = ReadStat(contig)
    insertion = Indel(contig)
    deletion = Indel(contig)
    mismatch = Mismatch(contig)
    splice = Splice(contig)

    genome = None
    if method not in ['cs', 'both']:
        genome = pysam.FastaFile(genome_file)

    cs_path = None
    if cs_dir is not None:
        cs_path = os.path.join(cs_dir, f'.readcs-{index:08d}.tmp')

    try:
        with _open_optional_cs(cs_path) as cs_fh:
            for record in records:
                cs = process_record(
                    record, genome,
                    readstat, insertion, deletion, mismatch, splice
                )
                if cs_fh is not None:
                    for low, high, mark, value in cs.get_contig_position():
                        cs_fh.write(
                            f'{record.query_name}\t{record.contig}\t'
                            f'{low}\t{high}\t{mark}\t{value}\n'
                        )
    finally:
        if genome is not None:
            genome.close()

    return StatBlock(
        contig, readstat, insertion, deletion, mismatch, splice, cs_path
    )


def stat_region(task, bam_file, genome_file, method, cs_dir=None):
    """Process one ``(index, contig, start, end)`` coordinate-window task.

    Opens the BAM itself, fetches the window, and assigns each read to the
    window containing its ``reference_start`` (a read that only overlaps the
    window's left edge belongs to the previous window and is skipped). The
    selected records are then accumulated by ``stat_records``.
    """
    index, contig, start, end = task
    file_type = bam_or_sam(bam_file)
    file_read = "rb" if file_type == "BAM" else "r"
    bam = pysam.AlignmentFile(bam_file, file_read)
    try:
        records = []
        for read in bam.fetch(contig, start, end):
            if read.reference_start < start:
                continue
            records.append(record_from_read(read, contig, method))
    finally:
        bam.close()
    return stat_records(
        (index, contig, records), genome_file, method, cs_dir
    )


def reduce_blocks_to_contigs(blocks, contigs):
    """Fold ``StatBlock`` chunks into per-contig ``ContigStats`` objects.

    ``sum`` concatenates each object's ``label``, so restore the label to the
    contig name afterwards — every block in a group shares that contig.
    """
    by_contig = defaultdict(list)
    for block in blocks:
        by_contig[block.contig].append(block)
    result = []
    for contig in contigs:
        contig_blocks = by_contig.get(contig)
        if not contig_blocks:
            continue
        readstat = sum(b.readstat for b in contig_blocks)
        insertion = sum(b.insertion for b in contig_blocks)
        deletion = sum(b.deletion for b in contig_blocks)
        mismatch = sum(b.mismatch for b in contig_blocks)
        splice = sum(b.splice for b in contig_blocks)
        for obj in (readstat, insertion, deletion, mismatch, splice):
            obj.label = contig
        result.append(ContigStats(
            readstat, insertion, deletion, mismatch, splice
        ))
    return result


def stat_element_from_bam_by_contig(bam_file,
                                    genome_file,
                                    contig,
                                    method='cs'):
    assert method in ['cs', 'MD', 'both'],\
        "method should be either: cs, MD, both."
    readstat = ReadStat(contig)
    insertion = Indel(contig)
    deletion = Indel(contig)
    mismatch = Mismatch(contig)
    splice = Splice(contig)

    file_type = bam_or_sam(bam_file)
    file_read = "rb" if file_type == "BAM" else "r"
    bam = pysam.AlignmentFile(bam_file, file_read)

    genome = None
    if method not in ['cs', 'both']:
        genome = pysam.FastaFile(genome_file)

    for read in bam.fetch(contig):
        record = record_from_read(read, contig, method)
        process_record(
            record, genome,
            readstat, insertion, deletion, mismatch, splice
        )

    bam.close()
    if genome is not None:
        genome.close()

    return ContigStats(readstat, insertion, deletion, mismatch, splice)


########################################
