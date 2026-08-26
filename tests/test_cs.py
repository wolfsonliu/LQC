import re

import pytest

from lqc import CS
from lqc.cs import cigar_to_list, cs_to_list, md_to_list

# ---- pure parser helpers ----

def test_cs_to_list():
    assert cs_to_list(':29*ga:19') == [
        [0, 29, ':', '29'],
        [29, 30, '*', 'ga'],
        [30, 49, ':', '19'],
    ]


def test_cigar_to_list():
    assert cigar_to_list('10M2D1I') == [
        [0, 10, 'M', '10'],
        [10, 12, 'D', '2'],
        [12, 12, 'I', '1'],
    ]


def test_md_to_list():
    assert md_to_list('5A2^CT3') == [
        ['5', 'M', 5],
        ['A', 'X', 1],
        ['2', 'M', 2],
        ['^CT', 'D', 2],
        ['3', 'M', 3],
    ]


# ---- CIGAR+MD -> cs conversion (covers mismatch / insertion / deletion) ----

@pytest.mark.parametrize(
    'cigar, md, read_seq, ref_seq, expected',
    [
        ('10M', '4C5', 'AAAAG' + 'A' * 5, 'AAAAC' + 'A' * 5, ':4*cg:5'),
        ('9M1I10M', '19', 'A' * 9 + 'C' + 'A' * 10, 'A' * 19, ':9+c:10'),
        ('5M2D5M', '5^CC5', 'A' * 10, 'AAAAA' + 'CC' + 'AAAAA', ':5-cc:5'),
    ],
)
def test_convert_cigar_md_to_cs_list(cigar, md, read_seq, ref_seq, expected):
    cs_obj = CS.from_cigar_string(
        cigar, md, read_seq, ref_seq,
        contig='chr1', start_pos=0, strand='+'
    )
    assert cs_obj.get_cs_tag_string() == expected


# ---- CS object: counts, positions, mutators ----

CS_STRING = ':10*ag:5+tt:3~gt100ag:4-ccc:2'


# NOTE: CS_STRING is already defined at module scope; do NOT redefine it here.


def _ref_cs_to_list(cs_string):
    """Reference copy of the pre-M1 regex tokenizer, for byte-parity only."""
    pos = 0
    cslenfuncs = {
        ':': int,
        '*': lambda x: 1,
        '+': lambda x: 0,
        '-': len,
        '~': lambda x: int(re.sub('[a-z]', '', x)),
    }
    cs_mark = re.sub('[0-9a-z]', ' ', cs_string).strip().split()
    cs_value = re.sub(r'[:*\-+~]', ' ', cs_string).strip().split()
    cslist = []
    for a, b in zip(cs_mark, cs_value, strict=True):
        low = pos
        pos += cslenfuncs[a](b)
        high = pos
        cslist.append([low, high, a, b])
    return cslist


def test_cs_to_list_all_marks():
    assert cs_to_list(CS_STRING) == [
        [0, 10, ':', '10'],
        [10, 11, '*', 'ag'],
        [11, 16, ':', '5'],
        [16, 16, '+', 'tt'],
        [16, 19, ':', '3'],
        [19, 119, '~', 'gt100ag'],
        [119, 123, ':', '4'],
        [123, 126, '-', 'ccc'],
        [126, 128, ':', '2'],
    ]


PARITY_STRINGS = [
    ':29*ga:19',
    CS_STRING,
    '*ag',
    ':10',
    '+atcg',
    '-c',
    '~gt100ag',
    ':4~ct636ac:5',
    ':29*ga:19*at:61~ct140ac:45*gc',
]


@pytest.mark.parametrize('cs_string', PARITY_STRINGS)
def test_cs_to_list_matches_reference(cs_string):
    assert cs_to_list(cs_string) == _ref_cs_to_list(cs_string)


def test_cs_to_list_matches_reference_on_fixtures(cs_test_data_records):
    for record in cs_test_data_records:
        cs_string = record[2]
        assert cs_to_list(cs_string) == _ref_cs_to_list(cs_string), record[0]


@pytest.fixture()
def cs_obj():
    return CS.from_cs_tag_string(
        CS_STRING, contig='chr1', start_pos=1000, strand='+'
    )


def test_cs_tag_string_round_trips(cs_obj):
    assert cs_obj.get_cs_tag_string() == CS_STRING


def test_read_length(cs_obj):
    assert cs_obj.get_read_length() == 27


def test_mismatch_counts(cs_obj):
    assert cs_obj.get_mismatch_count() == 1
    assert cs_obj.get_mismatch_type_count() == {'ag': 1}
    assert len(cs_obj.get_mismatches()) == 1


def test_insertion_counts(cs_obj):
    assert cs_obj.get_insertion_count() == 1
    assert cs_obj.get_insertion_length() == 2


def test_deletion_counts(cs_obj):
    assert cs_obj.get_deletion_count() == 1
    assert cs_obj.get_deletion_length() == 3


def test_get_indel_mismatches_groups(cs_obj):
    insertions, deletions, mismatches = cs_obj.get_indel_mismatches(
        coordinate='normalized_read'
    )
    # the grouped getter must be exactly equivalent to the three single-mark getters
    assert insertions == cs_obj.get_insertions(coordinate='normalized_read')
    assert deletions == cs_obj.get_deletions(coordinate='normalized_read')
    assert mismatches == cs_obj.get_mismatches(coordinate='normalized_read')
    assert [a[2] for a in insertions] == ['+']
    assert [a[2] for a in deletions] == ['-']
    assert [a[2] for a in mismatches] == ['*']
    assert insertions[0][3] == 'tt'
    assert deletions[0][3] == 'ccc'
    assert mismatches[0][3] == 'ag'


def test_indel_mismatch_normalized_eager_matches_lazy(cs_obj):
    # cached attribute is in cs-tag order (for CS_STRING: *, +, -)
    cached = cs_obj._indel_mismatch_normalized
    assert [a[2] for a in cached] == ['*', '+', '-']
    ins, dels, mism = cs_obj.get_indel_mismatches(coordinate='normalized_read')
    # the grouped getter must still match the single-mark lazy getters
    assert ins == cs_obj.get_insertions(coordinate='normalized_read')
    assert dels == cs_obj.get_deletions(coordinate='normalized_read')
    assert mism == cs_obj.get_mismatches(coordinate='normalized_read')
    # regroup the cached items by mark; must reproduce the grouped getter exactly
    by_mark = {'+': [], '-': [], '*': []}
    for a in cached:
        by_mark[a[2]].append(a)
    assert by_mark['+'] == ins
    assert by_mark['-'] == dels
    assert by_mark['*'] == mism


def test_get_indel_mismatches_reverse_strand_mirrors():
    rev = CS.from_cs_tag_string(
        CS_STRING, contig='chr1', start_pos=1000, strand='-'
    )
    ins, dels, mism = rev.get_indel_mismatches(coordinate='normalized_read')
    # read_len == 27; *ag at read_low=10, read_high=11 -> mirror to (16,17)
    assert len(mism) == 1
    assert mism[0][0] == (27 - 11) / 27
    assert mism[0][1] == (27 - 10) / 27
    assert mism[0][2] == '*'
    assert mism[0][3] == 'ag'
    # +tt at read_low=16, read_high=18 -> mirror to (9,11)
    assert len(ins) == 1
    assert ins[0][0] == (27 - 18) / 27
    assert ins[0][1] == (27 - 16) / 27
    # -ccc at read_low=25, read_high=25 -> mirror to (2,2)
    assert len(dels) == 1
    assert dels[0][0] == (27 - 25) / 27
    assert dels[0][1] == (27 - 25) / 27


def test_indel_mismatch_other_coordinates_unchanged(cs_obj):
    # non-normalized coordinates still go through the lazy path
    ins_c, dels_c, mism_c = cs_obj.get_indel_mismatches(coordinate='contig')
    assert [a[2] for a in ins_c] == ['+']
    assert dels_c[0][2] == '-'
    assert mism_c[0][2] == '*'
    assert mism_c[0][0] == 1010  # start_pos 1000 + ref_low 10


def test_intron_and_splice_counts(cs_obj):
    assert cs_obj.get_intron_count() == 1
    assert cs_obj.get_splice_pair_count() == {'gt-ag': 1}
    assert cs_obj.get_splice_site_count() == ({'gt': 1}, {'ag': 1})


def test_contig_position_adds_start_pos(cs_obj):
    assert cs_obj.get_contig_position()[0] == [1000, 1010, ':', '10']


def test_modify_contig_and_start_pos(cs_obj):
    assert cs_obj.modify_contig('chr2') == 'chr1'
    assert cs_obj.get_contig() == 'chr2'
    assert cs_obj.modify_start_pos(50) == 1000
    assert cs_obj.get_start_pos() == 50


# ---- committed fixture: cs path == CIGAR+MD path on every read ----

def test_cs_and_cigar_md_produce_identical_lists(cs_test_data_records):
    assert len(cs_test_data_records) > 0
    for record in cs_test_data_records:
        qname = record[0]
        cs_string, cigar, md, read_seq, ref_seq = record[2:7]
        bycs = CS.from_cs_tag_string(
            cs_string, contig='chr1', start_pos=0, strand='+'
        )
        bycigar = CS.from_cigar_string(
            cigar, md, read_seq, ref_seq,
            contig='chr1', start_pos=0, strand='+'
        )
        assert bycs._cs == bycigar._cs, qname


# ---- original long-read element-count regression ----

def test_element_count():
    cs_string = ':29*ga:19*at:61*ag:8*ga:19*ag:31*ag:15*ga:8*ta:4*ag:5*tc:42*ct:45*gc:130~ct140ac:45*gc:23~ct757ac:48*tc:104~ct659ac:24*tc:56*ga:11*gc:5+cccc:3*gc*at*gc:3+a:4*ta:23+tggtggtgc:23~ct88ac:2*ag:43*tc:48*ac:27*ag:18*tc:59~ct177ac:36+g:96*cg:3~ct239ac:110*ga:13*ca:12~ct172ac:47+g:5*ga:61*ct:30*ct~ct206ac:100*tc:11~ct6360ac:154~ct4420ac:29'

    cs_obj = CS.from_cs_tag_string(
        cs_string, contig='chr1', start_pos=1000, strand='+'
    )
    mismatch_type_count = cs_obj.get_mismatch_type_count()
    assert sum(mismatch_type_count.values()) == 33
    assert sum(mismatch_type_count.values()) == len(cs_obj.get_mismatches())
    assert sum(mismatch_type_count.values()) == cs_obj.get_mismatch_count()
    splice_p = cs_obj.get_splice_pair_count()
    assert sum(splice_p.values()) == 10
    assert splice_p == {'ct-ac': 10}
    assert sum(splice_p.values()) == len(cs_obj.get_introns())
    assert cs_obj.get_insertion_count() == 5
    assert cs_obj.get_insertion_length() == 16
    assert cs_obj.get_deletion_count() == 0
    assert cs_obj.get_deletion_length() == 0