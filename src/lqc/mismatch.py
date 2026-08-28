from collections import Counter, defaultdict

import numpy as np

from lqc._base import _LabelledStat
from lqc.utils import convert_complement


class Mismatch(_LabelledStat):
    """
    A class to store mismatch count information.
    """
    def __init__(self, label = ''):
        super().__init__(label)
        self._mismatch_types = []
        self._mismatch_index = {}
        self._type_count = Counter()
        self._type_locations = defaultdict(list)
        self._type_location_arrays = None  # {type: float64 array}

    def add_mismatch(self, mismatch,
                     normalized_read_location):
        if mismatch not in self._mismatch_index:
            self._mismatch_index[mismatch] = len(self._mismatch_types)
            self._mismatch_types.append(mismatch)
        self._type_count[mismatch] += 1
        self._type_locations[mismatch].append(
            normalized_read_location
        )
        self._type_location_arrays = None

    def _finalize(self):
        if self._type_location_arrays is None:
            self._type_location_arrays = {
                ty: np.asarray(locs, dtype = np.float64)
                for ty, locs in self._type_locations.items()
            }
            self._type_locations = None
        return self._type_location_arrays

    def get_total_count(self):
        return sum(self._type_count.values())

    def _get_bin_count(self, value_list, cuts):
        hist, edges = np.histogram(
            value_list, bins = cuts, density = False
        )
        bin_count = Counter()
        for i in range(len(hist)):
            if i < (len(hist) - 1):
                label = f'[{edges[i]},{edges[i+1]})'
            else:
                label = f'[{edges[i]},{edges[i+1]}]'
            bin_count[label] = hist[i]
        return bin_count

    def get_location_bin_count(self, cuts = None):
        if cuts is None:
            cuts = [0, 0.25, 0.5, 0.75, 1]
        bin_count = Counter()
        arrays = self._finalize()
        for ty in self._mismatch_types:
            bin_count += self._get_bin_count(
                arrays[ty], cuts = cuts
            )
        return bin_count

    def get_location_histogram(self, cuts = None):
        if cuts is None:
            cuts = [i / 10 for i in range(11)]
        arrays = self._finalize()
        total_hist = None
        edges = None
        for ty in self._mismatch_types:
            hist, edges = np.histogram(
                arrays[ty], bins = cuts, density = False
            )
            if total_hist is None:
                total_hist = hist.astype(np.float64)
            else:
                total_hist += hist
        if total_hist is None:
            total_hist, edges = np.histogram(
                [], bins = cuts, density = False
            )
        return edges, total_hist

    def get_location_bin_count_by_type(self, cuts = None):
        if cuts is None:
            cuts = [0, 0.25, 0.5, 0.75, 1]
        type_bin_count_dict = defaultdict(Counter)
        arrays = self._finalize()
        for ty in self._mismatch_types:
            type_bin_count_dict[ty] += self._get_bin_count(
                arrays[ty], cuts = cuts
            )
        return type_bin_count_dict

    def get_locations(self):
        arrays = self._finalize()
        locations = []
        for ty in self._mismatch_types:
            locations.extend(arrays[ty].tolist())
        return locations

    def convert_complement(self):
        newMis = type(self)(self.label)
        new_types = [
            convert_complement(a)
            for a in self._mismatch_types
        ]
        newMis._mismatch_types = new_types
        newMis._mismatch_index = {
            ty: i for i, ty in enumerate(new_types)
        }
        newMis._type_count = Counter({
            new_ty: self._type_count[old_ty]
            for old_ty, new_ty in zip(
                self._mismatch_types, new_types, strict = True
            )
        })
        arrays = self._finalize()
        newMis._type_location_arrays = {
            new_ty: arrays[old_ty].copy()
            for old_ty, new_ty in zip(
                self._mismatch_types, new_types, strict = True
            )
        }
        newMis._type_locations = None
        return newMis

    def get_type_count(self):
        return Counter(self._type_count)

    def __repr__(self):
        type_count = self.get_type_count()
        outstring = '\n'.join([
            f"Mismatch {self.label}:",
            '\n'.join(
                [f'  {a[0]}=>{a[1]}: {b}'
                 for a, b in type_count.items()]
            )
        ])
        return outstring

    def __str__(self):
        total_count = self.get_total_count()
        outstring = f"Mismatch {self.label}: {total_count} mismatches"
        return outstring

    def __add__(self, other):
        other = self._require_same_type(other)
        newMis = type(self)(f'{self.label} {other.label}')
        new_types = list(self._mismatch_types)
        type_index = {ty: i for i, ty in enumerate(new_types)}
        for ty in other._mismatch_types:
            if ty not in type_index:
                type_index[ty] = len(new_types)
                new_types.append(ty)
        s_arrs = self._finalize()
        o_arrs = other._finalize()
        empty = np.array([], dtype = np.float64)
        new_arrays = {}
        for ty in new_types:
            new_arrays[ty] = np.concatenate([
                s_arrs.get(ty, empty), o_arrs.get(ty, empty)
            ])
        newMis._mismatch_types = new_types
        newMis._mismatch_index = type_index
        newMis._type_count = self._type_count + other._type_count
        newMis._type_location_arrays = new_arrays
        newMis._type_locations = None
        return newMis

########################################
