import pysam
from mes.cs import CS
from mes.readstat import ReadStat
from mes.indel import Insertion, Deletion
from mes.mismatch import Mismatch
from mes.splice import Splice
from mes.utils import bam_or_sam


def stat_element_from_bam_by_contig(bam_file,
                                    genome_file,
                                    contig,
                                    method = 'cs',
                                    reverse_complement = True):
    assert method in ['cs', 'MD', 'both'],\
        "method should be either: cs, MD, both."
    # make counters
    readstat = {'+': ReadStat(contig), '-': ReadStat(contig)}
    insertion = {'+': Insertion(contig), '-': Insertion(contig)}
    deletion = {'+': Deletion(contig), '-': Deletion(contig)}
    mismatch = {'+': Mismatch(contig), '-': Mismatch(contig)}
    splice = {'+': Splice(contig), '-': Splice(contig)}

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
        readstat[strand].add_read(
            length = len(read.query_sequence),
            insertion = cs.get_insertion_count(),
            deletion = cs.get_deletion_count(),
            mismatch = cs.get_mismatch_count(),
            intron = cs.get_intron_count()
        )
        # insertion
        insertion[strand].add_insertion_list(
            [a[3] for a in cs.get_insertions()]
        )
        # deletion
        deletion[strand].add_deletion_list(
            [a[3] for a in cs.get_deletions()]
        )
        # mismatch
        mismatch[strand].add_count_dict(
            cs.get_mismatch_type_count()
        )
        # splice
        splice[strand].add_splice_pair_count_dict(
            cs.get_splice_pair_count()
        )

    bam.close()
    if method not in ['cs', 'both']:
        genome.close()
    else:
        pass

    if reverse_complement:
        # reverse complement
        insertion['-'] = insertion['-'].convert_reverse_complement()
        deletion['-'] = deletion['-'].convert_reverse_complement()
        mismatch['-'] = mismatch['-'].convert_complement()
        splice['-'] = splice['-'].convert_reverse_complement()
    else:
        pass

    return readstat, insertion, deletion, mismatch, splice


########################################
