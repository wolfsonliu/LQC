from collections import Counter

import numpy as np

from lqc._base import _LabelledStat
from lqc.utils import convert_reverse_complement


class Indel(_LabelledStat):
    """
    A class to store insertion or deletion count information.
    """
    def __init__(self, label = ''):
        super().__init__(label)
        self._indel_strings = []
        self._indel_index = {}
        self._indels = []
        self._indel_arrays = None  # (iidx int32, ilen int32, loc float64)

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
        self._indel_arrays = None

    def _finalize(self):
        if self._indel_arrays is None:
            rows = np.asarray(self._indels, dtype = np.float64).reshape(-1, 3)
            self._indel_arrays = (
                rows[:, 0].astype(np.int32, copy = False),
                rows[:, 1].astype(np.int32, copy = False),
                rows[:, 2],
            )
            self._indels = None
        return self._indel_arrays

    def get_indel_count(self):
        indel_count = Counter()
        iidx, _, _ = self._finalize()
        for i in iidx.tolist():
            indel_count[self._indel_strings[i]] += 1
        return indel_count

    def get_total_count(self):
        return int(self._finalize()[0].size)

    def get_total_length(self):
        return int(self._finalize()[1].sum())

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
        return self._finalize()[1].tolist()

    def get_locations(self):
        return self._finalize()[2].tolist()

    def convert_reverse_complement(self):
        newIndel = type(self)(self.label)
        iidx, ilen, loc = self._finalize()
        newIndel._indel_arrays = (iidx.copy(), ilen.copy(), loc.copy())
        newIndel._indels = None
        new_strings = [
            convert_reverse_complement(indel)
            for indel in self._indel_strings
        ]
        newIndel._indel_strings = new_strings
        newIndel._indel_index = {
            indel: i for i, indel in enumerate(new_strings)
        }
        return newIndel

    def get_location_bin_count(self, cuts = None):
        if cuts is None:
            cuts = [0, 0.25, 0.5, 0.75, 1]
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
        if cuts is None:
            cuts = [i / 10 for i in range(11)]
        hist, edges = np.histogram(
            self._finalize()[2], bins = cuts, density = False
        )
        return edges, hist

    def __add__(self, other):
        other = self._require_same_type(other)
        newIndel = type(self)(
            f'{self.label} {other.label}'
        )
        siidx, silen, sloc = self._finalize()
        oiidx, oilen, oloc = other._finalize()
        # string table: keep self's ordering, append other's new strings once
        new_strings = list(self._indel_strings)
        string_index = {
            indel: i for i, indel in enumerate(new_strings)
        }
        other_new_idx = []
        for indel in other._indel_strings:
            if indel in string_index:
                other_new_idx.append(string_index[indel])
            else:
                string_index[indel] = len(new_strings)
                new_strings.append(indel)
                other_new_idx.append(len(new_strings) - 1)
        remapped_oiidx = np.asarray(
            [other_new_idx[i] for i in oiidx.tolist()],
            dtype = np.int32,
        )
        newIndel._indel_arrays = (
            np.concatenate([siidx, remapped_oiidx]),
            np.concatenate([silen, oilen]),
            np.concatenate([sloc, oloc]),
        )
        newIndel._indels = None
        newIndel._indel_strings = new_strings
        newIndel._indel_index = string_index
        return newIndel

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
