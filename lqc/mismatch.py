import numpy as np
from copy import deepcopy
from collections import defaultdict
from collections import Counter
from lqc.utils import convert_complement


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
        self._mismatch_types = list()
        self._mismatches = list()

    def add_mismatch(self, mismatch,
                     normalized_read_location):
        midx = -1
        if mismatch in self._mismatch_types:
            midx, = np.where([
                a == mismatch
                for a in self._mismatch_types
            ])
            midx = midx[0]
        else:
            self._mismatch_types.append(mismatch)
            midx = len(self._mismatch_types) - 1
        self._mismatches.append([
            midx, normalized_read_location
        ])

    def get_total_count(self):
        return len(self._mismatches)

    def convert_complement(self):
        newMis = type(self)(self.label)
        newMis._mismatches = deepcopy(
            self._mismatches
        )
        new_types = list()
        for a in self._mismatch_types:
            new_types.append(
                convert_complement(a)
            )
        newMis._mismatch_types = new_types
        return newMis

    def _get_bin_count(self, value_list, cuts):
        hist, edges = np.histogram(
            value_list,
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

    def get_location_bin_count(self,
                               cuts = [0, 0.25, 0.5, 0.75, 1]):
        bin_count = self._get_bin_count(
            value_list = [b for a, b in self._mismatches],
            cuts = cuts
        )
        return bin_count

    def get_location_bin_count_by_type(self,
                                       cuts = [0, 0.25, 0.5, 0.75, 1]):
        type_bin_count_dict = defaultdict(
            Counter
        )
        type_idxs = list(set(
            [a for a, b in self._mismatches]
        ))
        for ty in type_idxs:
            type_list = [
                b for a, b in self._mismatches
                if a == ty
            ]
            type_bin_count_dict[
                self._mismatch_types[ty]
            ] += self._get_bin_count(
                type_list, cuts = cuts
            )
        return type_bin_count_dict

    def get_locations(self):
        return [
            loc
            for iidx, loc in self._mismatches
        ]

    def get_type_count(self):
        type_count = Counter()
        for tidx, loc in self._mismatches:
            thetype = self._mismatch_types[tidx]
            type_count[thetype] += 1
        return type_count

    def __repr__(self):
        type_count = self.get_type_count()
        outstring = '\n'.join([
            "Mismatch {}:".format(self.label),
            '\n'.join(
                ['  {}=>{}: {}'.format(a[0], a[1], b)
                 for a, b in type_count.items()]
            )
        ])
        return outstring

    def __str__(self):
        total_count = self.get_total_count()
        outstring = "Mismatch {}: {} mismatches".format(
            self.label, total_count
        )
        return outstring

    def __add__(self, other):
        assert type(other) == type(self), 'wrong object to add'
        newMis = type(self)(' '.join([self.label, other.label]))
        new_types = deepcopy(self._mismatch_types)
        for ty in other._mismatch_types:
            if ty not in new_types:
                new_types.append(ty)
            else:
                pass
        new_mismatches = deepcopy(self._mismatches)
        for tidx, loc in other._mismatches:
            type_label = other._mismatch_types[tidx]
            new_tidx, = np.where([
                a == type_label
                for a in new_types
            ])
            new_tidx = new_tidx[0]
            new_mismatches.append([new_tidx, loc])
        newMis._mismatch_types = new_types
        newMis._mismatches = new_mismatches
        return newMis

    def __radd__(self, other):
        if other == 0:
            return self
        else:
            return self.__add__(other)

########################################
