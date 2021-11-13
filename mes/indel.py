from collections import Counter
from itertools import accumulate


class Indel(object):
    """
    A class to store insertion or deletion count information.
    """
    def __init__(self, label = ''):
        super(Indel, self).__init__()
        if isinstance(label, str):
            pass
        else:
            raise TypeError(
                "label should be string."
            )
        self.label = label
        self._count_dict = Counter()
        self._total_count = 0
        self._total_length = 0

    def _add_string(self, string):
        self._count_dict[string] += 1
        self._total_length += len(string)
        self._total_count += 1

    def _add_string_list(self, string_list):
        for string in string_list:
            self._add_string(string)

    def _add_string_count_dict(self, count_dict):
        for string, count in count_dict.items():
            self._count_dict[string] += count
            self._total_length += (
                len(string) * count
            )
        self._total_count = sum(self._count_dict.values())

    def _get_count_dict(self):
        return self._count_dict

    def _get_total_count(self):
        return self._total_count

    def _get_total_length(self):
        return self._total_length

    def _get_mean_length(self):
        if self._total_count != 0:
            return self._total_length / self._total_count
        else:
            return 0

    def _get_median_length(self):
        count_dict = self._get_count_dict()
        count_list = sorted(
            list(count_dict.items()),
            key = lambda a: len(a[0])
        )
        count_list_cumsum = list(
            accumulate([b for a, b in count_list])
        )
        total_count = max(count_list_cumsum)
        median = 0
        if total_count % 2 == 0:
            median_idx1 = int(total_count / 2) - 1
            median_idx2 = int(total_count / 2)
            list_idx1 = 0
            list_idx2 = 0
            for i in range(1, len(count_list)):
                if (count_list_cumsum[i - 1] - 1 < median_idx1):
                    if (count_list_cumsum[i] - 1 >= median_idx1):
                        list_idx1 = i
                        break
                    else:
                        pass
                else:
                    pass
            for i in range(1, len(count_list)):
                if (count_list_cumsum[i - 1] - 1 < median_idx2):
                    if (count_list_cumsum[i] - 1 >= median_idx2):
                        list_idx2 = i
                        break
                    else:
                        pass
                else:
                    pass
            median = (
                len(count_list[list_idx1][0]) +
                len(count_list[list_idx2][0])
            ) / 2
        else:
            median_idx = int(total_count / 2)
            list_idx = 0
            for i in range(1, len(count_list)):
                if (count_list_cumsum[i - 1] - 1 < median_idx):
                    if (count_list_cumsum[i] - 1 >= median_idx):
                        list_idx = i
                        break
                    else:
                        pass
                else:
                    pass
            median = len(count_list[list_idx][0])
        return median

    def _get_longest_indel(self):
        count_dict = self._get_count_dict()
        max_length = max([
            len(a) for a in count_dict.keys()
        ])
        aims = [
            (a, b) for a, b in count_dict.items()
            if len(a) == max_length
        ]
        return aims

    def _get_most_abundant_indel(self):
        count_dict = self._get_count_dict()
        max_count = max(count_dict.values())
        aims = [
            (a, b) for a, b in count_dict.items()
            if b == max_count
        ]
        return aims

    def _convert_reverse_complement(self):
        ntpair = {'a': 't', 'c': 'g',
                  'g': 'c', 't': 'a',
                  'n': 'n'}
        new_dict = dict()
        for a, b in self._count_dict.items():
            new_string = ''.join([
                ntpair[c] for c in a[::-1]
            ])
            new_dict[new_string] = b
        new_item = type(self)(self.label)
        new_item._add_string_count_dict(
            new_dict
        )
        return new_item

    def __add__(self, other):
        assert other.__class__ == self.__class__, 'wrong object to add'
        sumIndel = type(self)(' '.join([self.label, other.label]))
        sumIndel._count_dict = self._count_dict + other._count_dict
        sumIndel._total_count = self._total_count + other._total_count
        sumIndel._total_length = self._total_length + other._total_length

        return sumIndel

    def __radd__(self, other):
        if other == 0:
            return self
        else:
            return self.__add__(other)


class Insertion(Indel):
    """
    A class to store insertion count information.
    """
    def __init__(self, label):
        super(Insertion, self).__init__(label)
        self.category = 'Insertion'

    def add_insertion(self, insertion):
        self._add_string(insertion)

    def add_insertion_list(self, insertion_list):
        self._add_string_list(insertion_list)

    def add_insertion_count_dict(self, count_dict):
        self._add_string_count_dict(count_dict)

    def get_count_dict(self):
        return self._get_count_dict()

    def get_total_count(self):
        return self._get_total_count()

    def get_total_insertion_length(self):
        return self._get_total_length()

    def get_mean_insertion_length(self):
        return self._get_mean_length()

    def get_median_insertion_length(self):
        return self._get_median_length()

    def get_longest_insertion(self):
        return self._get_longest_indel()

    def get_most_abundant_insertion(self):
        return self._get_most_abundant_indel()

    def convert_reverse_complement(self):
        return self._convert_reverse_complement()

    def __str__(self):
        outstring = "Insertion {}: {} insertions with total length {}".format(
            self.label,
            self.get_total_count(),
            self.get_total_insertion_length()
        )
        return outstring

    def __repr__(self):
        outstring = '\n'.join([
            "Insertion {}:".format(self.label),
            "  {} insertions".format(
                self.get_total_count()
            ),
            "  {} bp length in total".format(
                self.get_total_insertion_length()
            ),
            "  the longest insertion: {}".format(
                ', '.join([
                    '{} (count: {})'.format(a, b)
                    for a, b in self.get_longest_insertion()
                ])
            ),
            "  the most abundant insertion: {}".format(
                ', '.join([
                    '{} (count: {})'.format(a, b)
                    for a, b in self.get_most_abundant_insertion()
                ])
            )
        ])
        return outstring


class Deletion(Indel):
    """
    A class to store deletion count information.
    """
    def __init__(self, label):
        super(Deletion, self).__init__(label)
        self.category = 'Deletion'

    def add_deletion(self, deletion):
        self._add_string(deletion)

    def add_deletion_list(self, deletion_list):
        self._add_string_list(deletion_list)

    def add_deletion_count_dict(self, count_dict):
        self._add_string_count_dict(count_dict)

    def get_count_dict(self):
        return self._get_count_dict()

    def get_total_count(self):
        return self._get_total_count()

    def get_total_deletion_length(self):
        return self._get_total_length()

    def get_mean_deletion_length(self):
        return self._get_mean_length()

    def get_median_deletion_length(self):
        return self._get_median_length()

    def get_longest_deletion(self):
        return self._get_longest_indel()

    def get_most_abundant_deletion(self):
        return self._get_most_abundant_indel()

    def convert_reverse_complement(self):
        return self._convert_reverse_complement()

    def __str__(self):
        outstring = "Deletion {}: {} deletions with total length {}".format(
            self.label,
            self.get_total_count(),
            self.get_total_deletion_length()
        )
        return outstring

    def __repr__(self):
        outstring = '\n'.join([
            "Deletion {}:".format(self.label),
            "  {} deletions".format(
                self.get_total_count()
            ),
            "  {} bp length in total".format(
                self.get_total_deletion_length()
            ),
            "  the longest deletion: {}".format(
                ', '.join(
                    [
                        '{} (count: {})'.format(a, b)
                        for a, b in self.get_longest_deletion()
                    ]
                )
            ),
            "  the most abundant deletion: {}".format(
                ', '.join(
                    [
                        '{} (count: {})'.format(a, b)
                        for a, b in self.get_most_abundant_deletion()
                    ]
                )
            )
        ])
        return outstring

########################################
