from collections import Counter

from lqc._base import _LabelledStat
from lqc.utils import convert_reverse_complement


class Splice(_LabelledStat):
    """
    A class to store splice pair counts.
    """
    def __init__(self, label = ''):
        super().__init__(label)
        self._pair_count_dict = Counter()

    def add_splice_pair(self, splice_pair):
        """
        splice_pair string should be like: "gt-ag"
        """
        self._pair_count_dict[splice_pair] += 1

    def add_splice_pair_list(self, splice_pair_list):
        for pair in splice_pair_list:
            self.add_splice_pair(pair)

    def add_splice_pair_count_dict(self, count_dict):
        for string, count in count_dict.items():
            self._pair_count_dict[string] += count

    def get_splice_pair_count_dict(self):
        return self._pair_count_dict

    def get_total_splice_pair_count(self):
        return sum(self._pair_count_dict.values())

    def convert_reverse_complement(self):
        p_dict = self.get_splice_pair_count_dict()
        new_dict = {}
        for a, b in p_dict.items():
            new_string = convert_reverse_complement(a)
            new_dict[new_string] = b
        new_splice = type(self)(self.label)
        new_splice.add_splice_pair_count_dict(
            new_dict
        )
        return new_splice

    def get_most_abundant_splice_pair(self):
        p_dict = self.get_splice_pair_count_dict()
        max_count = max(p_dict.values())
        aim_pair = [
            (a, b)
            for a, b in p_dict.items()
            if b == max_count
        ]
        return aim_pair

    def __repr__(self):
        outstring = '\n'.join([
            f"Splice {self.label}:",
            f"  {self.get_total_splice_pair_count()} splice site pairs",
            "  the most abundant splice pair: {}".format(
                ', '.join([
                    f'{a} (count: {b})'
                    for a, b in self.get_most_abundant_splice_pair()
                ])
            )
        ])
        return outstring

    def __str__(self):
        outstring = f"Splice {self.label}: {self.get_total_splice_pair_count()} splice site pairs"
        return outstring

    def __add__(self, other):
        other = self._require_same_type(other)
        new_dict = (
            self.get_splice_pair_count_dict() +
            other.get_splice_pair_count_dict()
        )
        newSp = type(self)(
            f'{self.label} {other.label}'
        )
        newSp.add_splice_pair_count_dict(new_dict)
        return newSp

########################################
