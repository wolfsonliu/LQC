import pathlib

import pysam
import pytest

DATA_DIR = pathlib.Path(__file__).parent / 'data'


@pytest.fixture(scope='session')
def data_dir():
    return DATA_DIR


@pytest.fixture(scope='session')
def cs_test_data_records(data_dir):
    """Load the committed fixture: one tuple (qname, strand, cs, cigar, md, read_seq, ref_seq) per read."""
    path = data_dir / 'cs_test.test_data'
    with open(path, 'r') as fh:
        records = [tuple(line.strip().split()) for line in fh]
    return records


@pytest.fixture()
def cs_bam(tmp_path):
    """Write a tiny coordinate-sorted, indexed BAM carrying cs tags for integration tests."""
    bam_path = tmp_path / 'tiny.cs.bam'
    header = {
        'HD': {'VN': '1.6', 'SO': 'coordinate'},
        'SQ': [{'SN': 'chr1', 'LN': 1000000}],
    }
    with pysam.AlignmentFile(bam_path, 'wb', header=header) as out:
        a = pysam.AlignedSegment()
        a.query_name = 'read1'
        a.flag = 0
        a.reference_id = 0
        a.reference_start = 0
        a.mapping_quality = 60
        a.cigar = [(0, 10)]
        a.query_sequence = 'ACGTACGTAC'
        a.query_qualities = pysam.qualitystring_to_array('I' * 10)
        a.set_tag('cs', ':5*ag:4')
        out.write(a)

        b = pysam.AlignedSegment()
        b.query_name = 'read2'
        b.flag = 0
        b.reference_id = 0
        b.reference_start = 100
        b.mapping_quality = 60
        b.cigar = [(0, 10)]
        b.query_sequence = 'ACGTACGTAC'
        b.query_qualities = pysam.qualitystring_to_array('I' * 10)
        b.set_tag('cs', ':10')
        out.write(b)
    pysam.index(str(bam_path))
    return str(bam_path)


@pytest.fixture()
def cs_bam_with_indel(tmp_path):
    """Tiny indexed BAM with one read carrying one deletion (-cc) and one insertion (+tt)."""
    bam_path = tmp_path / 'indel.cs.bam'
    header = {
        'HD': {'VN': '1.6', 'SO': 'coordinate'},
        'SQ': [{'SN': 'chr1', 'LN': 1000000}],
    }
    with pysam.AlignmentFile(bam_path, 'wb', header=header) as out:
        a = pysam.AlignedSegment()
        a.query_name = 'read1'
        a.flag = 0
        a.reference_id = 0
        a.reference_start = 0
        a.mapping_quality = 60
        a.cigar = [(0, 3), (2, 2), (0, 4), (1, 2), (0, 2)]
        a.query_sequence = 'AAAAABBBBCC'
        a.query_qualities = pysam.qualitystring_to_array('I' * 11)
        a.set_tag('cs', ':3-cc:4+tt:2')
        out.write(a)
    pysam.index(str(bam_path))
    return str(bam_path)