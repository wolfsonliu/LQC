import sys
import pysam
from collections import defaultdict

genome = pysam.FastaFile(sys.argv[1])
cssam = pysam.AlignmentFile(sys.argv[2], 'rb')
mdsam = pysam.AlignmentFile(sys.argv[3], 'rb')

outfile = open(sys.argv[4] + '.test_data', 'w')

reads = defaultdict(dict)

for read in mdsam.fetch('chr1', 0, 200000):
    reads[read.qname]['read_strand'] = '-' if read.is_reverse else '+'
    reads[read.qname]['read_seq'] = read.query_sequence
    reads[read.qname]['ref_seq'] = genome.fetch(
        'chr1', read.reference_start, read.reference_end
    )
    reads[read.qname]['cigarstring'] = read.cigarstring
    reads[read.qname]['mdstring'] = [
        a[1] for a in read.tags if a[0] == 'MD'
    ][0]

for read in cssam.fetch('chr1', 0, 200000):
    if reads[read.qname]['cigarstring'] == read.cigarstring:
        reads[read.qname]['csstring'] = [
            a[1] for a in read.tags if a[0] == 'cs'
        ][0]
    else:
        reads.pop(read.qname)

for qname in reads:
    outfile.write(
        '{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
            qname,
            reads[qname]['read_strand'],
            reads[qname]['csstring'],
            reads[qname]['cigarstring'],
            reads[qname]['mdstring'],
            reads[qname]['read_seq'],
            reads[qname]['ref_seq']
        )
    )


genome.close()
outfile.close()
cssam.close()
mdsam.close()
