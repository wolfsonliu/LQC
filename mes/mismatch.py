from collections import Counter

class Mismatch(object):
    """
    A class to store mismatch count information.
    """
    def __init__(self, label = ''):
        super(Mismatch, self).__init__()
        if isinstance(label, str):
            pass
        else:
            raise TypeError(
                "label should be string."
            )
        self.label = label
        self._mismatch_count = Counter()

    def add_count(self, mismatch_type, mismatch_count):
        self._mismatch_count[mismatch_type] += mismatch_count

    def add_count_dict(self, count_dict):
        for a, b in count_dict.items():
            self._mismatch_count[a] += b

    def get_count_dict(self):
        return self._mismatch_count

    def get_total_count(self):
        return sum(self._mismatch_count.values())

    def convert_complement(self):
        ntpair = {'a': 't', 'c': 'g',
                  'g': 'c', 't': 'a',
                  'n': 'n', '=': '=', '>': '>'}
        new_dict = dict()
        for a, b in self._mismatch_count.items():
            new_string = ''.join([
                ntpair[c] for c in a
            ])
            new_dict[new_string] = b
        new_mismatch = type(self)(self.label)
        new_mismatch.add_count_dict(
            new_dict
        )
        return new_mismatch

    def __repr__(self):
        outstring = '\n'.join([
            "Mismatch {}:".format(self.label),
            '\n'.join([
                "  {}: {}".format(a, b)
                for a, b in self._mismatch_count.items()
            ])
        ])
        return outstring

    def __str__(self):
        total_counts = self.get_total_count()
        outstring = "Mismatch {}: {} mismatches".format(
            self.label, total_counts
        )
        return outstring

    def __add__(self, other):
        assert other.__class__ == self.__class__, 'wrong object to add'
        count_dict = self.get_count_dict()
        other_dict = other.get_count_dict()
        sumMis = type(self)(' '.join([self.label, other.label]))
        sumMis.add_count_dict(count_dict + other_dict)
        return sumMis

    def __radd__(self, other):
        if other == 0:
            return self
        else:
            return self.__add__(other)

########################################
