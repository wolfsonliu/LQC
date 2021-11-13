import re
import pysam


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


########################################
