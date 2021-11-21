import pysam


def convert_complement(string):
    ntpair = {'a': 't', 'c': 'g',
              'g': 'c', 't': 'a',
              'n': 'n', '-': '-'}
    return ''.join([ntpair[c] for c in string])


def convert_reverse_complement(string):
    return convert_complement(string)[::-1]


def check_cs_md_tag(tag_list):
    '''
    The tag_list should be a list of tuples (tag, value).
    '''
    return [a for a, b in tag_list
            if a == 'cs' or a == 'MD']


def bam_or_sam(file_path):
    expansion = file_path.split('.')[-1].upper()
    assert expansion in ['BAM', 'SAM'], 'Not a bam or sam file expansion.'
    return expansion


def check_bam_with_cs_or_md(bam_file):
    file_type = bam_or_sam(bam_file)
    file_read = "rb" if file_type == "BAM" else "r"

    bam = pysam.AlignmentFile(bam_file, file_read)
    i = 0
    read_cs_md = list()
    bam_type = None
    for read in bam:
        i += 1
        if i >= 10:
            break
        else:
            rcsmd = check_cs_md_tag(read.tags)
            if len(rcsmd) > 0:
                if "cs" in rcsmd and "MD" in rcsmd:
                    bam_type = "both"
                    return bam_type
                elif "cs" in rcsmd:
                    if bam_type == "MD":
                        bam_type = "both"
                        return bam_type
                    else:
                        bam_type = "cs"
                elif "MD" in rcsmd:
                    if bam_type == "cs":
                        bam_type = "both"
                        return bam_type
                    else:
                        bam_type = "MD"
                else:
                    pass
            else:
                pass
    bam.close()
    return bam_type


def write_readcs(bam_file,
                 genome_file,
                 output_file,
                 method = 'cs'):
    assert method in ['cs', 'MD', 'both'],\
        "method should be either: cs, MD, both."
    file_type = bam_or_sam(bam_file)
    file_read = "rb" if file_type == "BAM" else "r"
    bam = pysam.AlignmentFile(bam_file, file_read)

    if method not in ['cs', 'both']:
        genome = pysam.FastaFile(genome_file)
    else:
        pass

    output = open(output_file, 'w')
    output.write(
        '\t'.join([
            'read_name', 'contig',
            'low', 'high',
            'cs_mark', 'cs_value'
        ]) + '\n'
    )
    for read in bam:
        strand = '-' if read.is_reverse else '+'
        if method == 'cs':
            # there're cs tags in the bam file
            cs_string = [
                a[1] for a in read.tags
                if a[0] == 'cs'
            ][0]
            cs = CS.from_cs_tag_string(
                cs_tag_string = cs_string,
                contig = read.reference_name,
                start_pos = read.reference_start,
                strand = strand
            )
        else:
            # there's no cs tag in the bam file, and there're MD tags.
            read_seq = read.query_sequence
            ref_seq = genome.fetch(
                read.reference_name,
                read.reference_start,
                read.reference_end
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
                contig = read.reference_name,
                start_pos = read.reference_start,
                strand = strand
            )
        # read cs
        cs_list = cs.get_contig_position()
        for line in cs_list:
            output.write(
                '\t'.join(
                    [read.query_name,
                     read.reference_name] +
                    ['{}'.format(a) for a in line]
                ) + '\n'
            )
    output.close()
    bam.close()
    if method not in ['cs', 'both']:
        genome.close()
    else:
        pass
    return

########################################
