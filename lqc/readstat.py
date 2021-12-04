from copy import deepcopy


class ReadStat(object):
    """
    A class to store read count and read statistic information.
    """
    def __init__(self, label = ''):
        super(ReadStat, self).__init__()
        if isinstance(label, str):
            pass
        else:
            raise TypeError(
                "label should be string."
            )
        self.label = label
        self._read_count = 0
        self._total_base = 0
        self._reads = []

    def add_read(self,
                 length,
                 insertion,
                 deletion,
                 mismatch,
                 intron):
        self._read_count += 1
        self._reads.append([
            length,
            insertion, deletion,
            mismatch, intron
        ])
        self._total_base += length

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
        return [a[0] for a in self._reads]

    def get_insertions(self):
        return [a[1] for a in self._reads]

    def get_deletions(self):
        return [a[2] for a in self._reads]

    def get_mismatches(self):
        return [a[3] for a in self._reads]

    def get_introns(self):
        return [a[4] for a in self._reads]

    def get_length_normalized_insertions(self):
        return [a[1] / a[0] for a in self._reads]

    def get_length_normalized_deletions(self):
        return [a[2] / a[0] for a in self._reads]

    def get_length_normalized_mismatches(self):
        return [a[3] / a[0] for a in self._reads]

    def get_length_normalized_introns(self):
        return [a[4] / a[0] for a in self._reads]

    def get_length_NL(self, percent):
        assert percent >= 0 and percent <= 100,\
            "percent value should be between 0 and 100."
        lengths = self.get_lengths()
        lengths = sorted(lengths, reverse = True)
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
            "  read counts: {}".format(
                self.get_read_count()
            ),
            "  total bases: {}".format(
                self.get_total_base()
            ),
            "  mean of read length: {:.4}".format(
                float(self.get_mean_length())
            ),
            "  median of read length: {:.4}".format(
                float(self.get_median_length())
            ),
            "  [min, max] of read length: [{}, {}]".format(
                self.get_min_length(), self.get_max_length()
            ),
            "  N50 of read length: {}".format(N50),
            "  L50 of read length: {}".format(L50),
            "  mean of insertions per read: {:.4}".format(
                float(self.get_mean_insertions())
            ),
            "  mean of deletions per read: {:.4}".format(
                float(self.get_mean_deletions())
            ),
            "  mean of mismatches per read: {:.4}".format(
                float(self.get_mean_mismatches())
            ),
            "  mean of introns per read: {:.4}".format(
                float(self.get_mean_introns())
            )
        ])
        return outstring

    def __str__(self):
        outstring = "ReadStat: {} reads".format(
            self.get_read_count()
        )
        return outstring

    def __add__(self, other):
        assert type(other) == type(self),\
            'wrong object to add'
        sumReadStat = type(self)(
            ' '.join([self.label, other.label])
        )
        sumReadStat._reads = deepcopy(
            self._reads + other._reads
        )
        sumReadStat._read_count = deepcopy(
            self._read_count + other._read_count
        )
        sumReadStat._total_base = deepcopy(
            self._total_base + other._total_base
        )
        return sumReadStat

    def __radd__(self, other):
        if other == 0:
            return self
        else:
            return self.__add__(other)

    def __len__(self):
        return self._read_count

########################################
