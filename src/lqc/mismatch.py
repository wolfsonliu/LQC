from collections import Counter, defaultdict

import numpy as np

from lqc.utils import convert_complement


class Mismatch:
    """
    A class to store mismatch count information.
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
        self._mismatch_types = []
        self._mismatch_index = {}
        self._type_count = Counter()
        self._type_locations = defaultdict(list)

    def add_mismatch(self, mismatch,
                     normalized_read_location):
        if mismatch not in self._mismatch_index:
            self._mismatch_index[mismatch] = len(self._mismatch_types)
            self._mismatch_types.append(mismatch)
        self._type_count[mismatch] += 1
        self._type_locations[mismatch].append(
            normalized_read_location
        )

    def get_total_count(self):
        return sum(self._type_count.values())

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
        newMis._type_locations = defaultdict(list, {
            new_ty: list(self._type_locations[old_ty])
            for old_ty, new_ty in zip(
                self._mismatch_types, new_types, strict = True
            )
        })
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
                label = f'[{edges[i]},{edges[i+1]})'
            else:
                label = f'[{edges[i]},{edges[i+1]}]'
            bin_count[label] = hist[i]
        return bin_count

    def get_location_bin_count(self, cuts = None):
        if cuts is None:
            cuts = [0, 0.25, 0.5, 0.75, 1]
        bin_count = Counter()
        for ty in self._mismatch_types:
            bin_count += self._get_bin_count(
                self._type_locations[ty], cuts = cuts
            )
        return bin_count

    def get_location_histogram(self, cuts = None):
        """Return (edges, counts) for the overall normalized-location
        histogram, summed across mismatch types without materializing the
        full location list."""
        if cuts is None:
            cuts = [i / 10 for i in range(11)]
        total_hist = None
        edges = None
        for ty in self._mismatch_types:
            hist, edges = np.histogram(
                self._type_locations[ty],
                bins = cuts,
                density = False
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
        type_bin_count_dict = defaultdict(
            Counter
        )
        for ty in self._mismatch_types:
            type_bin_count_dict[ty] += self._get_bin_count(
                self._type_locations[ty], cuts = cuts
            )
        return type_bin_count_dict

    def get_locations(self):
        locations = []
        for ty in self._mismatch_types:
            locations.extend(self._type_locations[ty])
        return locations

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
        assert isinstance(other, type(self)), 'wrong object to add'
        newMis = type(self)(f'{self.label} {other.label}')
        new_types = list(self._mismatch_types)
        type_index = {ty: i for i, ty in enumerate(new_types)}
        for ty in other._mismatch_types:
            if ty not in type_index:
                type_index[ty] = len(new_types)
                new_types.append(ty)

        new_locs = defaultdict(list)
        for ty in new_types:
            merged = list(self._type_locations.get(ty, ()))
            merged.extend(other._type_locations.get(ty, ()))
            new_locs[ty] = merged

        newMis._mismatch_types = new_types
        newMis._mismatch_index = type_index
        newMis._type_count = self._type_count + other._type_count
        newMis._type_locations = new_locs
        return newMis

    def __radd__(self, other):
        if other == 0:
            return self
        else:
            return self.__add__(other)

########################################