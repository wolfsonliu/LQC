import os
import unittest
from lqc import CS

file_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(file_dir)


class TestCS(unittest.TestCase):

    def test_init_methods(self):
        test_file = os.path.join(
            file_dir, 'cs_test.test_data'
        )
        with open(test_file, 'r') as f:
            for line in f:
                linesplit = line.strip().split()
                cs_string = linesplit[2]
                cigar_string = linesplit[3]
                md_string = linesplit[4]
                read_seq = linesplit[5]
                ref_seq = linesplit[6]
                bycs = CS.from_cs_tag_string(
                    cs_tag_string = cs_string,
                    contig = 'chr1',
                    start_pos = 0,
                    strand = '+'
                )
                bycigar = CS.from_cigar_string(
                    cigar_string = cigar_string,
                    md_string = md_string,
                    read_seq = read_seq,
                    ref_seq = ref_seq,
                    contig = 'chr1',
                    start_pos = 0,
                    strand = '+'
                )
                self.assertEqual(bycs._cs, bycigar._cs)

    def test_element_count(self):
        cs_string = ':29*ga:19*at:61*ag:8*ga:19*ag:31*ag:15*ga:8*ta:4*ag:5*tc:42*ct:45*gc:130~ct140ac:45*gc:23~ct757ac:48*tc:104~ct659ac:24*tc:56*ga:11*gc:5+cccc:3*gc*at*gc:3+a:4*ta:23+tggtggtgc:23~ct88ac:2*ag:43*tc:48*ac:27*ag:18*tc:59~ct177ac:36+g:96*cg:3~ct239ac:110*ga:13*ca:12~ct172ac:47+g:5*ga:61*ct:30*ct~ct206ac:100*tc:11~ct6360ac:154~ct4420ac:29'

        cs = CS.from_cs_tag_string(
            cs_tag_string = cs_string,
            contig = 'chr1',
            start_pos = 1000,
            strand = '+'
        )
        mismatch_type_count = cs.get_mismatch_type_count()
        self.assertEqual(
            sum(b for a, b in mismatch_type_count.items()),
            33
        )
        self.assertEqual(
            sum(b for a, b in mismatch_type_count.items()),
            len(cs.get_mismatches())
        )
        self.assertEqual(
            sum(b for a, b in mismatch_type_count.items()),
            cs.get_mismatch_count()
        )
        splice_p = cs.get_splice_pair_count()
        splice_l, splice_r = cs.get_splice_site_count()
        self.assertEqual(
            sum(b for a, b in splice_p.items()),
            10
        )
        self.assertEqual(
            sum(b for a, b in splice_p.items()),
            len(cs.get_introns())
        )
        insertion_count = cs.get_insertion_count()
        insertion_length = cs.get_insertion_length()
        self.assertEqual(insertion_count, 5)
        self.assertEqual(insertion_length, 16)
        deletion_count = cs.get_deletion_count()
        deletion_length = cs.get_deletion_length()
        self.assertEqual(deletion_count, 0)
        self.assertEqual(deletion_length, 0)


if __name__ == '__main__':
    unittest.main()

########################################
