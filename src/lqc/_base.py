"""Internal base for the four labelled accumulator classes."""

from copy import copy


class _LabelledStat:
    def __init__(self, label=''):
        if not isinstance(label, str):
            raise TypeError('label should be string.')
        self.label = label

    def __radd__(self, other):
        if other == 0:
            return self
        return self.__add__(other)

    def _require_same_type(self, other):
        if not isinstance(other, type(self)):
            raise TypeError('wrong object to add')
        return other


def concat_stats(iterable, label):
    """Fold accumulator objects (via ``__add__``) into one relabelled object.

    Always returns a new object so relabelling never aliases a caller's
    per-contig instance (the exact bug the old ``sum([x]) is x`` path guarded
    against). Only the first element is shallow-copied, sharing the read-only
    numpy arrays rather than duplicating them.
    """
    it = iter(iterable)
    try:
        acc = copy(next(it))
    except StopIteration:
        raise ValueError('concat_stats requires at least one object') from None
    for obj in it:
        acc = acc + obj
    acc.label = label
    return acc