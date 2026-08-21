"""Optional smoke tests against the real, gitignored BAM.

These run only when ``tmp/data/ENCFF417VHJ.chr22.sorted.bam`` is present, so CI
(which does not ship the ~200 MB file) skips them automatically while a local
checkout exercises the parser against real minimap2 `--cs` output.
"""

import pathlib

import pysam
import pytest

from lqc.cs import CS
from lqc.utils import check_bam_with_cs_or_md, list_bam_contigs

REPO_ROOT = pathlib.Path(__file__).resolve().parents[1]
REAL_BAM = REPO_ROOT / 'tmp' / 'data' / 'ENCFF417VHJ.chr22.sorted.bam'

pytestmark = pytest.mark.skipif(
    not REAL_BAM.exists(),
    reason='real BAM not present (gitignored tmp/data)',
)


def _sample_reads(limit=2000):
    """Yield up to ``limit`` reads spread across chr22 without reading it all."""
    bam = pysam.AlignmentFile(str(REAL_BAM), 'rb')
    try:
        mapped = {
            stat.contig: stat.mapped for stat in bam.get_index_statistics()
        }
        step = max(1, mapped.get('chr22', 0) // limit)
        for i, read in enumerate(bam.fetch('chr22')):
            if i % step == 0:
                yield read
    finally:
        bam.close()


def test_real_bam_has_cs_tags():
    assert check_bam_with_cs_or_md(str(REAL_BAM)) == 'cs'


def test_real_bam_contigs():
    assert list_bam_contigs(str(REAL_BAM)) == ['chr22']


def test_real_bam_cs_read_length_matches_cigar():
    """cs read length must equal the CIGAR query-consuming length (M + I)."""
    for read in _sample_reads():
        cs = CS.from_cs_tag_string(
            next(a[1] for a in read.tags if a[0] == 'cs'),
            contig='chr22',
            start_pos=read.reference_start,
            strand='-' if read.is_reverse else '+',
        )
        cigar_query_len = sum(
            length for op, length in read.cigartuples if op in (0, 1)
        )
        assert cs.get_read_length() == cigar_query_len


def test_real_bam_cs_reference_span_matches_cigar():
    """cs reference span must equal the CIGAR reference length (M + D + N)."""
    for read in _sample_reads():
        cs = CS.from_cs_tag_string(
            next(a[1] for a in read.tags if a[0] == 'cs'),
            contig='chr22',
            start_pos=read.reference_start,
            strand='-' if read.is_reverse else '+',
        )
        cigar_ref_len = sum(
            length for op, length in read.cigartuples if op in (0, 2, 3)
        )
        positions = cs.get_contig_position()
        assert positions[-1][1] - positions[0][0] == cigar_ref_len