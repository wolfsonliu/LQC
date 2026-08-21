import sys
from collections import defaultdict

import pysam

genome = pysam.FastaFile(sys.argv[1])
cssam = pysam.AlignmentFile(sys.argv[2], 'rb')
mdsam = pysam.AlignmentFile(sys.argv[3], 'rb')

reads = defaultdict(dict)

for read in mdsam.fetch('chr1', 0, 200000):
    reads[read.qname]['read_strand'] = '-' if read.is_reverse else '+'
    reads[read.qname]['read_seq'] = read.query_sequence
    reads[read.qname]['ref_seq'] = genome.fetch(
        'chr1', read.reference_start, read.reference_end
    )
    reads[read.qname]['cigarstring'] = read.cigarstring
    reads[read.qname]['mdstring'] = next(
        a[1] for a in read.tags if a[0] == 'MD'
    )

for read in cssam.fetch('chr1', 0, 200000):
    if reads[read.qname]['cigarstring'] == read.cigarstring:
        reads[read.qname]['csstring'] = next(
            a[1] for a in read.tags if a[0] == 'cs'
        )
    else:
        reads.pop(read.qname)

with open(sys.argv[4] + '.test_data', 'w') as outfile:
    outfile.writelines(
        '{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(
            qname,
            reads[qname]['read_strand'],
            reads[qname]['csstring'],
            reads[qname]['cigarstring'],
            reads[qname]['mdstring'],
            reads[qname]['read_seq'],
            reads[qname]['ref_seq']
        )
        for qname in reads
    )

genome.close()
cssam.close()
mdsam.close()