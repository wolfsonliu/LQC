import pysam
from lqc.cs import CS
from lqc.readstat import ReadStat
from lqc.indel import Indel
from lqc.mismatch import Mismatch
from lqc.splice import Splice
from lqc.utils import bam_or_sam
from lqc.utils import convert_reverse_complement
from lqc.utils import convert_complement


def stat_element_from_bam_by_contig(bam_file,
                                    genome_file,
                                    contig,
                                    method = 'cs'):
    assert method in ['cs', 'MD', 'both'],\
        "method should be either: cs, MD, both."
    # make counters
    readstat = ReadStat(contig)
    insertion = Indel(contig)
    deletion = Indel(contig)
    mismatch = Mismatch(contig)
    splice = Splice(contig)

    # open bam file
    file_type = bam_or_sam(bam_file)
    file_read = "rb" if file_type == "BAM" else "r"
    bam = pysam.AlignmentFile(bam_file, file_read)

    if method not in ['cs', 'both']:
        genome = pysam.FastaFile(genome_file)
    else:
        pass

    for read in bam.fetch(contig):
        strand = '-' if read.is_reverse else '+'
        if method == 'cs':
            # there're cs tags in the bam file
            cs_string = [
                a[1] for a in read.tags
                if a[0] == 'cs'
            ][0]
            cs = CS.from_cs_tag_string(
                cs_tag_string = cs_string,
                contig = contig,
                start_pos = read.reference_start,
                strand = strand
            )
        else:
            # there's no cs tag in the bam file, and there're MD tags.
            read_seq = read.query_sequence
            ref_seq = genome.fetch(
                contig, read.reference_start, read.reference_end
            )
            md_string = [
                a[1] for a in read.tags
                if a[0] == 'MD'
            ][0]
            cs = CS.from_cigar_string(
                cigar_string = read.cigarstring,
                md_string = md_string,
                read_seq = read_seq,
                ref_seq = ref_seq,
                contig = contig,
                start_pos = read.reference_start,
                strand = strand
            )
        # read stat
        readstat.add_read(
            length = len(read.query_sequence),
            insertion = cs.get_insertion_count(),
            deletion = cs.get_deletion_count(),
            mismatch = cs.get_mismatch_count(),
            intron = cs.get_intron_count()
        )
        # insertion
        for a, b, c, d in cs.get_insertions(
                coordinate='normalized_read'
        ):
            if strand == '+':
                indel_string = d
            else:
                indel_string = convert_reverse_complement(d)
            insertion.add_indel(
                indel = indel_string,
                normalized_read_location = a
            )
        # deletion
        for a, b, c, d in cs.get_insertions(
                coordinate = 'normalized_read'
        ):
            if strand == '+':
                indel_string = d
            else:
                indel_string = convert_reverse_complement(d)
            deletion.add_indel(
                indel = indel_string,
                normalized_read_location = a
            )
        # mismatch
        for a, b, c, d in cs.get_mismatches(
                coordinate = 'normalized_read'
        ):
            if strand == '+':
                mis_string = d
            else:
                mis_string = convert_complement(d)
            mismatch.add_mismatch(
                mismatch = mis_string,
                normalized_read_location = a
            )
        # splice
        splice_pair_count = cs.get_splice_pair_count()
        for s_str, s_ct in splice_pair_count.items():
            if strand == '+':
                splice_string = s_str
            else:
                splice_string = convert_reverse_complement(s_str)
            splice.add_splice_pair_count_dict(
                {splice_string: s_ct}
            )

    bam.close()
    if method not in ['cs', 'both']:
        genome.close()
    else:
        pass

    return readstat, insertion, deletion, mismatch, splice


########################################
