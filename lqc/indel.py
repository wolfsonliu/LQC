import numpy as np
from copy import deepcopy
from collections import Counter
from lqc.utils import convert_reverse_complement


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
        self._indel_strings = list()
        self._indels = list()

    def add_indel(self, indel,
                  normalized_read_location):
        iidx = -1
        if indel in self._indel_strings:
            iidx, = np.where([
                a == indel
                for a in self._indel_strings
            ])
            iidx = iidx[0]
        else:
            self._indel_strings.append(indel)
            iidx = len(self._indel_strings) - 1
        self._indels.append([
            iidx, len(indel),
            normalized_read_location
        ])

    def get_indel_count(self):
        indel_count = Counter()
        for iidx, ilen, loc in self._indels:
            indel = self._indel_strings[iidx]
            indel_count[indel] += 1
        return indel_count

    def get_total_count(self):
        return len(self._indels)

    def get_total_length(self):
        return sum([
            ilen
            for iidx, ilen, loc in self._indels
        ])

    def get_mean_length(self):
        total_count = self.get_total_count()
        total_length = self.get_total_length()
        if total_count != 0:
            return total_length / total_count
        else:
            return 0

    def get_median_length(self):
        len_list = sorted(self.get_lengths())
        total_count = len(len_list)
        median = 0
        if total_count % 2 == 0:
            median_idx1 = int(total_count / 2) - 1
            median_idx2 = int(total_count / 2)
            median = (
                len_list[median_idx1] +
                len_list[median_idx2]
            ) / 2
        else:
            median_idx = int(total_count / 2)
            median = float(len_list[median_idx])
        return median

    def get_longest_indel(self):
        indel_count = self.get_indel_count()
        len_list = [
            len(indel)
            for indel in indel_count
        ]
        max_length = max(len_list)
        aims = [
            indel
            for indel in indel_count
            if len(indel) == max_length
        ]
        return aims

    def get_most_abundant_indel(self):
        indel_count = self.get_indel_count()
        max_count = max(indel_count.values())
        aims = [
            indel
            for indel,count in indel_count.items()
            if count == max_count
        ]
        return aims

    def get_lengths(self):
        return [
            ilen
            for iidx, ilen, loc in self._indels
        ]

    def get_locations(self):
        return [
            loc
            for iidx, ilen, loc in self._indels
        ]

    def convert_reverse_complement(self):
        newIndel = type(self)(self.label)
        newIndel._indels = deepcopy(self._indels)
        new_strings = list()
        for indel in self._indel_strings:
            new_strings.append(
                convert_reverse_complement(indel)
            )

        newIndel._indel_strings = new_strings
        return newIndel

    def get_location_bin_count(self,
                               cuts = [0, 0.25, 0.5, 0.75, 1]):
        hist, edges = np.histogram(
            self.get_locations(),
            bins = cuts,
            density = False
        )
        bin_count = Counter()
        for i in range(len(hist)):
            if i < (len(hist) - 1):
                label = '[{},{})'.format(
                    edges[i], edges[i+1]
                )
            else:
                label = '[{},{}]'.format(
                    edges[i], edges[i+1]
                )
            bin_count[label] = hist[i]
        return bin_count

    def __add__(self, other):
        assert type(other) == type(self), 'wrong object to add'
        newIndel = type(self)(
            ' '.join([self.label, other.label])
        )
        new_strings = deepcopy(self._indel_strings)
        other_new_idx = list()
        for indel in other._indel_strings:
            if indel not in new_strings:
                new_strings.append(indel)
                other_new_idx.append(
                    len(new_strings) - 1
                )
            else:
                iidx, = np.where([
                    a == indel
                    for a in new_strings
                ])
                iidx = iidx[0]
                other_new_idx.append(iidx)
        new_indels = deepcopy(self._indels)
        for iidx, ilen, loc in other._indels:
            new_iidx = other_new_idx[iidx]
            new_indels.append([
                new_iidx, ilen, loc
            ])
        newIndel._indel_strings = new_strings
        newIndel._indels = new_indels
        return newIndel

    def __radd__(self, other):
        if other == 0:
            return self
        else:
            return self.__add__(other)

    def __str__(self):
        outstring = "Indel {}: {} indels with total length {}".format(
            self.label,
            self.get_total_count(),
            self.get_total_length()
        )
        return outstring

    def __repr__(self):
        outstring = '\n'.join([
            "Indel {}:".format(self.label),
            "  {} indels".format(
                self.get_total_count()
            ),
            "  {} bp length in total".format(
                self.get_total_length()
            ),
            "  the mean indel length: {}".format(
                self.get_mean_length()
            ),
            "  the median indel length: {}".format(
                self.get_median_length()
            )
        ])
        return outstring


########################################
