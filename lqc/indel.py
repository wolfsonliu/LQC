from collections import Counter
from copy import deepcopy

import numpy as np

from lqc.utils import convert_reverse_complement


class Indel:
    """
    A class to store insertion or deletion count information.
    """
    def __init__(self, label = ''):
        super().__init__()
        if isinstance(label, str):
            pass
        else:
            raise TypeError(
                "label should be string."
            )
        self.label = label
        self._indel_strings = list()
        self._indel_index = dict()
        self._indels = list()

    def add_indel(self, indel,
                  normalized_read_location):
        iidx = self._indel_index.get(indel)
        if iidx is None:
            iidx = len(self._indel_strings)
            self._indel_strings.append(indel)
            self._indel_index[indel] = iidx
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
        if total_count == 0:
            return 0
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
        if not indel_count:
            return []
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
        if not indel_count:
            return []
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
        new_strings = [
            convert_reverse_complement(indel)
            for indel in self._indel_strings
        ]
        newIndel._indel_strings = new_strings
        newIndel._indel_index = {
            indel: i for i, indel in enumerate(new_strings)
        }
        return newIndel

    def get_location_bin_count(self,
                               cuts = [0, 0.25, 0.5, 0.75, 1]):
        edges, hist = self.get_location_histogram(cuts = cuts)
        bin_count = Counter()
        for i in range(len(hist)):
            if i < (len(hist) - 1):
                label = f'[{edges[i]},{edges[i+1]})'
            else:
                label = f'[{edges[i]},{edges[i+1]}]'
            bin_count[label] = hist[i]
        return bin_count

    def get_location_histogram(self, cuts = None):
        """Return (edges, counts) for the normalized-location histogram
        without materializing the raw locations as a Python list."""
        if cuts is None:
            cuts = [i / 10 for i in range(11)]
        locs = np.fromiter(
            (loc for iidx, ilen, loc in self._indels),
            dtype = float,
            count = len(self._indels)
        )
        hist, edges = np.histogram(
            locs, bins = cuts, density = False
        )
        return edges, hist

    def __add__(self, other):
        assert type(other) == type(self), 'wrong object to add'
        newIndel = type(self)(
            ' '.join([self.label, other.label])
        )
        new_strings = deepcopy(self._indel_strings)
        string_index = {
            indel: i for i, indel in enumerate(new_strings)
        }
        other_new_idx = list()
        for indel in other._indel_strings:
            if indel in string_index:
                other_new_idx.append(string_index[indel])
            else:
                string_index[indel] = len(new_strings)
                new_strings.append(indel)
                other_new_idx.append(len(new_strings) - 1)
        new_indels = deepcopy(self._indels)
        for iidx, ilen, loc in other._indels:
            new_iidx = other_new_idx[iidx]
            new_indels.append([
                new_iidx, ilen, loc
            ])
        newIndel._indel_strings = new_strings
        newIndel._indels = new_indels
        newIndel._indel_index = string_index
        return newIndel

    def __radd__(self, other):
        if other == 0:
            return self
        else:
            return self.__add__(other)

    def __str__(self):
        outstring = f"Indel {self.label}: {self.get_total_count()} indels with total length {self.get_total_length()}"
        return outstring

    def __repr__(self):
        outstring = '\n'.join([
            f"Indel {self.label}:",
            f"  {self.get_total_count()} indels",
            f"  {self.get_total_length()} bp length in total",
            f"  the mean indel length: {self.get_mean_length()}",
            f"  the median indel length: {self.get_median_length()}"
        ])
        return outstring


########################################
