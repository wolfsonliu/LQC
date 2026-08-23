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
    query_name: str
    cs_string: Optional[str] = None
    cigar: Optional[str] = None
    md_string: Optional[str] = None
    query_sequence: Optional[str] = None
    reference_end: Optional[int] = None


def record_from_read(read, contig, method):
    """Extract the minimal fields from a pysam read for stat accumulation."""
    strand = '-' if read.is_reverse else '+'
    if method == 'cs':
        return ReadRecord(
            method='cs',
            contig=contig,
            start_pos=read.reference_start,
            strand=strand,
            query_length=len(read.query_sequence),
            query_name=read.query_name,
            cs_string=read.get_tag('cs'),
        )
    return ReadRecord(
        method='md',
        contig=contig,
        start_pos=read.reference_start,
        strand=strand,
        query_length=len(read.query_sequence),
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
    """Return ``[(contig, [ReadRecord, ...]), ...]`` in BAM traversal order."""
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

    if genome is not None:
        genome.close()

    return contig, readstat, insertion, deletion, mismatch, splice, cs_path


def reduce_blocks_to_contigs(block_results, contigs):
    """Fold per-chunk ``(contig, rs, ins, dele, mis, spl)`` into per-contig tuples."""
    by_contig = defaultdict(list)
    for contig, rs, ins, dele, mis, spl in block_results:
        by_contig[contig].append((rs, ins, dele, mis, spl))
    result = []
    for contig in contigs:
        blocks = by_contig.get(contig)
        if not blocks:
            continue
        rss, inss, deles, miss, spls = zip(*blocks, strict=True)
        result.append((
            sum(rss), sum(inss), sum(deles), sum(miss), sum(spls)
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

    return readstat, insertion, deletion, mismatch, splice


########################################