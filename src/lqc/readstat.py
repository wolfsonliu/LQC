import numpy as np


class ReadStat:
    """
    A class to store read count and read statistic information.
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
        self._read_count = 0
        self._total_base = 0
        self._total_aligned_base = 0
        self._reads = []
        self._reads_array = None
        self._sorted_lengths = None

    def add_read(self,
                 length,
                 insertion,
                 deletion,
                 mismatch,
                 intron,
                 mapping_quality = 0,
                 aligned_length = 0):
        self._read_count += 1
        self._reads.append([
            length,
            insertion, deletion,
            mismatch, intron,
            mapping_quality, aligned_length
        ])
        self._total_base += length
        self._total_aligned_base += aligned_length
        self._reads_array = None
        self._sorted_lengths = None

    def _finalize(self):
        """Convert the boxed row list to an ``(n, 7)`` int32 array (lazily) and
        drop the boxed list. Columns: length, insertion, deletion, mismatch,
        intron, mapping_quality, aligned_length."""
        if self._reads_array is None:
            self._reads_array = np.asarray(
                self._reads, dtype = np.int32
            ).reshape(-1, 7)
            self._reads = None
        return self._reads_array

    def get_read_count(self):
        return self._read_count

    def get_read_count_with_n_insertions(self, n):
        elements = [
            a for a in self.get_insertions()
            if a == n
        ]
        return len(elements)

    def get_read_count_with_n_deletions(self, n):
        elements = [
            a for a in self.get_deletions()
            if a == n
        ]
        return len(elements)

    def get_read_count_with_n_mismatches(self, n):
        elements = [
            a for a in self.get_mismatches()
            if a == n
        ]
        return len(elements)

    def get_total_base(self):
        return self._total_base

    def get_lengths(self):
        return self._finalize()[:, 0].tolist()

    def get_insertions(self):
        return self._finalize()[:, 1].tolist()

    def get_deletions(self):
        return self._finalize()[:, 2].tolist()

    def get_mismatches(self):
        return self._finalize()[:, 3].tolist()

    def get_introns(self):
        return self._finalize()[:, 4].tolist()

    def get_length_normalized_insertions(self):
        arr = self._finalize()
        return (arr[:, 1] / arr[:, 0]).tolist()

    def get_length_normalized_deletions(self):
        arr = self._finalize()
        return (arr[:, 2] / arr[:, 0]).tolist()

    def get_length_normalized_mismatches(self):
        arr = self._finalize()
        return (arr[:, 3] / arr[:, 0]).tolist()

    def get_length_normalized_introns(self):
        arr = self._finalize()
        return (arr[:, 4] / arr[:, 0]).tolist()

    def _sorted_lengths_desc(self):
        if self._sorted_lengths is None:
            self._sorted_lengths = sorted(
                self.get_lengths(), reverse = True
            )
        return self._sorted_lengths

    def get_length_NL(self, percent):
        assert percent >= 0 and percent <= 100,\
            "percent value should be between 0 and 100."
        lengths = self._sorted_lengths_desc()
        basesum = 0
        previous_basesum = 0
        percentbase = self._total_base * percent / 100
        for i in range(self._read_count):
            length = lengths[i]
            previous_basesum = basesum
            basesum += length
            if (previous_basesum <= percentbase) \
               and (basesum >= percentbase):
                N = length
                L = i + 1
                break
            else:
                continue
        return N, L

    def get_mean_length(self):
        return self._total_base / self._read_count

    def get_min_length(self):
        return min(self.get_lengths())

    def get_max_length(self):
        return max(self.get_lengths())

    def insertions_per_base(self):
        return sum(
            self.get_insertions()
        ) / self.get_total_base()

    def get_mean_insertions(self):
        return sum(
            self.get_insertions()
        ) / self._read_count

    def get_mean_length_normalized_insertions(self):
        return sum(
            self.get_length_normalized_insertions()
        ) / self._read_count

    def get_min_insertions(self):
        return min(self.get_insertions())

    def get_max_insertions(self):
        return max(self.get_insertions())

    def deletions_per_base(self):
        return sum(
            self.get_deletions()
        ) / self.get_total_base()

    def get_mean_deletions(self):
        return sum(
            self.get_deletions()
        ) / self._read_count

    def get_mean_length_normalized_deletions(self):
        return sum(
            self.get_length_normalized_deletions()
        ) / self._read_count

    def get_min_deletions(self):
        return min(self.get_deletions())

    def get_max_deletions(self):
        return max(self.get_deletions())

    def mismatches_per_base(self):
        return sum(
            self.get_mismatches()
        ) / self.get_total_base()

    def get_mean_mismatches(self):
        return sum(
            self.get_mismatches()
        ) / self._read_count

    def get_mean_length_normalized_mismatches(self):
        return sum(
            self.get_length_normalized_mismatches()
        ) / self._read_count

    def get_min_mismatches(self):
        return min(self.get_mismatches())

    def get_max_mismatches(self):
        return max(self.get_mismatches())

    def get_mean_introns(self):
        return sum(
            self.get_introns()
        ) / self._read_count

    def get_mean_length_normalized_introns(self):
        return sum(
            self.get_length_normalized_introns()
        ) / self._read_count

    def get_min_introns(self):
        return min(self.get_introns())

    def get_max_introns(self):
        return max(self.get_introns())

    def get_mapping_qualities(self):
        return self._finalize()[:, 5].tolist()

    def get_aligned_lengths(self):
        return self._finalize()[:, 6].tolist()

    def get_aligned_fractions(self):
        arr = self._finalize()
        length = arr[:, 0]
        aligned = arr[:, 6]
        fraction = np.divide(
            aligned, length,
            out = np.zeros_like(length, dtype = np.float64),
            where = length > 0,
        )
        return fraction.tolist()

    def get_total_aligned_base(self):
        return self._total_aligned_base

    def get_mean_mapping_quality(self):
        qualities = self.get_mapping_qualities()
        if not qualities:
            return 0.0
        return sum(qualities) / len(qualities)

    def get_median_mapping_quality(self):
        return self._get_median(self.get_mapping_qualities())

    def get_mean_aligned_fraction(self):
        fractions = self.get_aligned_fractions()
        if not fractions:
            return 0.0
        return sum(fractions) / len(fractions)

    def get_median_aligned_fraction(self):
        return self._get_median(self.get_aligned_fractions())

    def get_read_count_with_aligned_fraction_below(self, threshold = 0.9):
        arr = self._finalize()
        length = arr[:, 0]
        aligned = arr[:, 6]
        fraction = np.divide(
            aligned, length,
            out = np.zeros_like(length, dtype = np.float64),
            where = length > 0,
        )
        return int(np.sum(fraction < threshold))

    def get_read_count_fully_aligned(self):
        arr = self._finalize()
        return int(np.sum(arr[:, 6] == arr[:, 0]))

    def insertions_per_aligned_base(self):
        if self._total_aligned_base == 0:
            return 0.0
        return sum(self.get_insertions()) / self._total_aligned_base

    def deletions_per_aligned_base(self):
        if self._total_aligned_base == 0:
            return 0.0
        return sum(self.get_deletions()) / self._total_aligned_base

    def mismatches_per_aligned_base(self):
        if self._total_aligned_base == 0:
            return 0.0
        return sum(self.get_mismatches()) / self._total_aligned_base

    def _get_median(self, item_list):
        sorted_item_list = sorted(item_list)
        median = 0
        if self._read_count % 2 == 0:
            median = (
                sorted_item_list[int(self._read_count / 2) - 1] +
                sorted_item_list[int(self._read_count / 2)]
            ) / 2
        else:
            median = sorted_item_list[int(self._read_count / 2)]
        return median

    def get_median_length(self):
        return self._get_median(self.get_lengths())

    def get_median_insertion(self):
        return self._get_median(self.get_insertions())

    def get_median_deletion(self):
        return self._get_median(self.get_deletions())

    def get_median_mismatch(self):
        return self._get_median(self.get_mismatches())

    def get_median_intron(self):
        return self._get_median(self.get_introns())

    def __repr__(self):
        N50, L50 = self.get_length_NL(50)
        outstring = '\n'.join([
            "ReadStat:",
            f"  read counts: {self.get_read_count()}",
            f"  total bases: {self.get_total_base()}",
            f"  mean of read length: {float(self.get_mean_length()):.4}",
            f"  median of read length: {float(self.get_median_length()):.4}",
            f"  [min, max] of read length: [{self.get_min_length()}, {self.get_max_length()}]",
            f"  N50 of read length: {N50}",
            f"  L50 of read length: {L50}",
            f"  mean of insertions per read: {float(self.get_mean_insertions()):.4}",
            f"  mean of deletions per read: {float(self.get_mean_deletions()):.4}",
            f"  mean of mismatches per read: {float(self.get_mean_mismatches()):.4}",
            f"  mean of introns per read: {float(self.get_mean_introns()):.4}",
            f"  mean of aligned fraction: {float(self.get_mean_aligned_fraction()):.4}",
            f"  median of mapping quality: {float(self.get_median_mapping_quality()):.4}"
        ])
        return outstring

    def __str__(self):
        outstring = f"ReadStat: {self.get_read_count()} reads"
        return outstring

    def __add__(self, other):
        assert isinstance(other, type(self)),\
            'wrong object to add'
        sumReadStat = type(self)(
            f'{self.label} {other.label}'
        )
        sumReadStat._reads_array = np.concatenate([
            self._finalize(), other._finalize()
        ])
        sumReadStat._reads = None
        sumReadStat._read_count = self._read_count + other._read_count
        sumReadStat._total_base = self._total_base + other._total_base
        sumReadStat._total_aligned_base = self._total_aligned_base + other._total_aligned_base
        return sumReadStat

    def __radd__(self, other):
        if other == 0:
            return self
        else:
            return self.__add__(other)

    def __len__(self):
        return self._read_count

########################################
